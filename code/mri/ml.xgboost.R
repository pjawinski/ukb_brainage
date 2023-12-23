#!/usr/bin/env Rscript
# =======================================
# === Estimate brainAGE using xgboost === 
# =======================================

# get script location
initial.options = commandArgs(trailingOnly = FALSE)
file.arg.name = "--file="
script.name = sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename = dirname(script.name)

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=6) {
  stop(paste0('expected 6 arguments, but ', length(args), ' arguments provided.'), call.=FALSE)
}

# get arguments from command line
tissue = args[1] # tissue="gm"
targetDir = args[2] # targetDir="results/mri/ml.xgb/gm"
Rfile = args[3] # Rfile="results/mri/prepML.discovery.RData"
matFile = args[4] # matFile="results/mri/prepML.discovery.mat"
matlabpath = args[5] # matlabpath="/opt/matlab/bin/matlab"
threads = as.numeric(args[6]) # threads = 56

message(paste0('\n--- XGBoost settings ---',
               '\ntissue: ', tissue,
               '\ntargetDir: ', targetDir,
               '\nRfile: ', Rfile,
               '\nmatFile: ', matFile,
               '\nmatlabpath: ', matlabpath,
               '\nthreads: ', threads, '\n'))

# load required packages
for (pkg in c('data.table','doParallel','foreach','xgboost')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# create targetDir and clear logfile in case it already exists
system(sprintf('mkdir -p %s', targetDir))
system(sprintf('rm -f %s/logfile.txt', targetDir))

# load data (cat12-preprocessed T1, 'gm', 'wm', 'covs', and 'u')
Rfile = normalizePath(Rfile)
matFile = normalizePath(matFile)
load(Rfile)

# get features and target variable
x = get(tissue)
y = data.frame(covs$t1_age)

# set number of folds and repeats
nfolds = max(u)
nrepeats = ncol(u)

# create output array for age estimates
pred = array(NA, list(nrow(y),nrepeats,2))

# start pool for parallel computations
cl = makeCluster(nfolds)
registerDoParallel(cl)

# loop across folds and repeats
for (j in 1:nrepeats) {
    message(sprintf('Starting with repeat %03d.', j))
    start_time = Sys.time()

    # parallel matlab pca
    message(' - running matlab pca for each fold.')
    foreach(i=c(1:nfolds)) %dopar% {
      
      # run pca in matlab
      system(sprintf('mkdir -p %s/models/xgmodel_%03d_%02d', targetDir, j, i))
      system(sprintf('cp %s/ml_pcamatlab.m %s/models/xgmodel_%03d_%02d/ml_pcamatlab.m', script.basename, targetDir, j, i))
      system(sprintf('cd %s/models/xgmodel_%03d_%02d/; %s -nodesktop -nodisplay -r \'matFile = \"%s\"; tissue = \"%s\"; i = %d; j = %d; maxNumCompThreads(%d); run(\"ml_pcamatlab.m\");\' | tail +11', targetDir, j, i, matlabpath, matFile, tissue, i, j, max(floor(threads/nfolds),1)))
    
    }
      
    # train and apply models with xgboost
    for (i in 1:nfolds) {
      
      # y variable: select train, validation, and test sample
      message(sprintf('Starting with xgmodel_%03d_%02d.', j, i))
      Train_Targets = data.frame(y[u[,j]!=i & u[,j]!=v[i,j],1])
      Validation_Targets = data.frame(y[u[,j]==v[i,j],1])
      Test_Targets = data.frame(y[u[,j]==i,1])
      
      # x variable: import results from matlab pca
      Train_Samples_pca = data.frame(fread(sprintf('%s/models/xgmodel_%03d_%02d/Train_Samples_pca.txt', targetDir, j, i), sep = "\t", header = F))
      Validation_Samples_pca = data.frame(fread(sprintf('%s/models/xgmodel_%03d_%02d/Validation_Samples_pca.txt', targetDir, j, i), sep = "\t", header = F))
      Test_Samples_pca = data.frame(fread(sprintf('%s/models/xgmodel_%03d_%02d/Test_Samples_pca.txt', targetDir, j, i), sep = "\t", header = F))

      # convert do xgb dmatrix
      dtrain = xgb.DMatrix(data = as.matrix(Train_Samples_pca), label=as.matrix(Train_Targets))
      dvalid = xgb.DMatrix(data = as.matrix(Validation_Samples_pca), label=as.matrix(Validation_Targets))
      dtest = xgb.DMatrix(data = as.matrix(Test_Samples_pca), label=as.matrix(Test_Targets))
        
      # train models
      message(' - training tree-based model.')
      watchlist = list(train=dtrain, test=dvalid)
      model_tree = xgb.train(data = dtrain,
                              watchlist = watchlist,
                              params = list(booster = 'gbtree'),
                              eval_metric = "mae",
                              max.depth = 3,
                              eta = 0.02,
                              nthread = threads,
                              nrounds = 5000,
                              early_stopping_rounds = 50,
                              verbose = 0,
                              objective = "reg:linear")
        
      message(' - training linear model.')
      model_linear = xgb.train(data = dtrain,
                                watchlist = watchlist,
                                params = list(booster = 'gblinear'),
                                eval_metric = "mae",
                                eta = 0.02,
                                nthread = threads,
                                nrounds = 5000,
                                early_stopping_rounds = 50,
                                verbose = 0,
                                objective = "reg:linear")
        
        # make predictions for test sample
        pred[u[,j]==i,j,1] = predict(model_tree, dtest)
        pred[u[,j]==i,j,2] = predict(model_linear, dtest)
        
        # save models and predictions
        message(' - saving model predictions.')
        save(model_tree, model_linear, file = sprintf('%s/models/xgmodel_%03d_%02d/models.RData', targetDir, j, i))
        pred_fold = array(0L, list(length(y[,1]),2))
        pred_fold[u[,j]==i,1] = pred[u[,j]==i,j,1]
        pred_fold[u[,j]==i,2] = pred[u[,j]==i,j,2]
        write.table(pred_fold, sprintf('%s/models/xgmodel_%03d_%02d/pred_fold.txt', targetDir, j, i), sep = '\t', row.names = F, col.names = F)

        # delete .m and .txt files
        system(sprintf('rm -f %s/models/xgmodel_%03d_%02d/*.m', targetDir, j, i))
        system(sprintf('rm -f %s/models/xgmodel_%03d_%02d/*_pca.txt', targetDir, j, i))        
    }
  
    # to monitor progress type into terminal: watch -n 1 tail {yourtargetDir}/logfile.txt
    line = sprintf("Repeat: %d - Corr Tree: %.4f - Corr Linear: %.4f", j, cor(pred[,j,1],y[,1]), cor(pred[,j,2],y[,1]))
    message(line)
    write(line,file=paste0(targetDir, '/logfile.txt'), append=TRUE)
}

# stop cluster
stopCluster(cl)

# save predicted age
message('Saving predicted age.')
brainage = data.frame(IID = covs$IID, brainage_tree = rowMeans(pred[,1:100,1], na.rm = T), brainage_lin = rowMeans(pred[,1:100,2], na.rm = T))
write.table(pred, sprintf('%s/pred_brainage.txt', targetDir), sep="\t", col.names = F, row.names = F, quote = F)
write.table(brainage, sprintf('%s/pred_brainage_means.txt', targetDir), sep="\t", col.names = T, row.names = F, quote = F)
save(pred, brainage, file = sprintf('%s/pred_brainage.RData', targetDir))
message('--- Completed: XGBoost ---\n')
