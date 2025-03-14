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
targetDir = args[2] # targetDir="results/mri/ml.xgb/singlemodel/"
Rfile = args[3] # Rfile="results/mri/prepML.discovery.RData"
matFile = args[4] # matFile="results/mri/prepML.discovery.mat"
matlabpath = args[5] # matlabpath="/opt/matlab/bin/matlab"
threads = as.numeric(args[6]) # threads = 56

message(paste0('\n--- XGBoost settings (singlemodel) ---',
               '\ntissue: ', tissue,
               '\ntargetDir: ', targetDir,
               '\nRfile: ', Rfile,
               '\nmatFile: ', matFile,
               '\nmatlabpath: ', matlabpath,
               '\nthreads: ', threads, '\n'))

# load required packages
for (pkg in c('data.table','xgboost')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# load data (cat12-preprocessed T1, 'gm', 'wm', 'covs', and 'u')
Rfile = normalizePath(Rfile)
matFile = normalizePath(matFile)
load(Rfile)

# get features and target variable
x = get(tissue)
y = data.frame(covs$t1_age)

# set test and validation dataset as defined in first iteration (repeat 1, fold 1)
i = 1
j = 1
    
# run pca in matlab
system(sprintf('mkdir -p %s', targetDir))
message(' - running pca in matlab.')
system(sprintf('cp %s/ml_pcamatlab_singlemodel.m %s/ml_pcamatlab_singlemodel.m', script.basename, targetDir))
system(sprintf('cd %s/; %s -nodesktop -nodisplay -r \'matFile = \"%s\"; tissue = \"%s\"; i = %d; j = %d; maxNumCompThreads(%d); run(\"ml_pcamatlab_singlemodel.m\");\' | tail +11', targetDir, matlabpath, matFile, tissue, i, j, threads))

# y variable: select train and validation sample
message(' - selecting training and validation samples.')
Train_Targets = data.frame(y[u[,j]!=v[i,j],1])
Validation_Targets = data.frame(y[u[,j]==v[i,j],1])

# x variable: import results from matlab pca
Train_Samples_pca = data.frame(fread(sprintf('%s/xgb_%s_train_samples_pca.txt', targetDir, tissue), sep = "\t", header = F))
Validation_Samples_pca = data.frame(fread(sprintf('%s/xgb_%s_validation_samples_pca.txt', targetDir, tissue), sep = "\t", header = F))

# convert do xgb dmatrix
dtrain = xgb.DMatrix(data = as.matrix(Train_Samples_pca), label=as.matrix(Train_Targets))
dvalid = xgb.DMatrix(data = as.matrix(Validation_Samples_pca), label=as.matrix(Validation_Targets))
  
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
                        nrounds = 10000,
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
                          nrounds = 10000,
                          early_stopping_rounds = 50,
                          verbose = 0,
                          objective = "reg:linear")
# delete .m and gzip .txt files
system(sprintf('rm -f %s/*.m', targetDir))
system(sprintf('gzip -dc %s/xgb_%s_train_samples_pca.txt', targetDir, tissue))   
system(sprintf('gzip -dc %s/xgb_%s_validation_samples_pca.txt', targetDir, tissue))   

# save models
message(' - saving models.')
save(model_tree, model_linear, file = sprintf('%s/xgb_%s_models.RData', targetDir, tissue), version = 2)
message('--- Completed: XGBoost (singlemodel) ---\n')
