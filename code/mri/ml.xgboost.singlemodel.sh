#!/bin/bash

# ==========================================
# ======== xgboost machine learning ========
# ==========================================

# -----------------
# run pca in Matlab
# -----------------

# enter cluster
jawinskp@cluster1.psychologie.hu-berlin.de
cd /home/groups/markett/UK_Biobank/11b_brainage

# open matlab
tmux
/opt/matlab/bin/matlab -nodesktop -nodisplay

% load data
parent = '/home/groups/markett/UK_Biobank/11b_brainage';
load(strcat(parent,'/04_primary_dataset.mat'))

% make directory
system(char(strcat({'mkdir -p '}, parent, '/05_xgboost/whole_sample')));

% -----------
% Grey matter
% -----------

% define samples
i=1;
j = 1;
data = gm;

Train_Samples = data(u(:,j)~= v(i,j),:); % Train_Samples = gm(u(:,j)~=i & u(:,j)~= v(i,j),:);
Validation_Samples = data(u(:,j)==v(i,j),:);

% run pca
[Train_Samples_pca_coeff,Train_Samples_pca_score] = pca(Train_Samples);

% get pca scores of the first 500 principal components
Train_Samples_pca_score = Train_Samples_pca_score(:,1:500);
Train_Samples_pca_coeff = Train_Samples_pca_coeff(:,1:500);
Train_Samples_means = mean(Train_Samples);

% calculate pca scores for validation and test sample 
Validation_Samples_centered = Validation_Samples-repmat(Train_Samples_means,size(Validation_Samples,1),1);
Validation_Samples_pca_score = Validation_Samples_centered/Train_Samples_pca_coeff'; %'

% save data
dlmwrite(strcat(parent, '/05_xgboost/whole_sample/xgb_gm_Train_Samples_pca.txt'), Train_Samples_pca_score, 'delimiter', '\t', 'precision', 16);
dlmwrite(strcat(parent, '/05_xgboost/whole_sample/xgb_gm_Validation_Samples_pca.txt'), Validation_Samples_pca_score, 'delimiter', '\t', 'precision', 16);
save(strcat(parent, '/05_xgboost/whole_sample/xgb_gm_pca_training.mat'), 'Train_Samples_means', 'Train_Samples_pca_coeff')

% ------------
% White matter
% ------------

% define samples
i=1;
j = 1;
data = wm;

Train_Samples = data(u(:,j)~= v(i,j),:); % Train_Samples = gm(u(:,j)~=i & u(:,j)~= v(i,j),:);
Validation_Samples = data(u(:,j)==v(i,j),:);

% run pca
[Train_Samples_pca_coeff,Train_Samples_pca_score] = pca(Train_Samples);

% get pca scores of the first 500 principal components
Train_Samples_pca_score = Train_Samples_pca_score(:,1:500);
Train_Samples_pca_coeff = Train_Samples_pca_coeff(:,1:500);
Train_Samples_means = mean(Train_Samples);

% calculate pca scores for validation and test sample 
Validation_Samples_centered = Validation_Samples-repmat(Train_Samples_means,size(Validation_Samples,1),1);
Validation_Samples_pca_score = Validation_Samples_centered/Train_Samples_pca_coeff'; %'

% save data
dlmwrite(strcat(parent, '/05_xgboost/whole_sample/xgb_wm_Train_Samples_pca.txt'), Train_Samples_pca_score, 'delimiter', '\t', 'precision', 16);
dlmwrite(strcat(parent, '/05_xgboost/whole_sample/xgb_wm_Validation_Samples_pca.txt'), Validation_Samples_pca_score, 'delimiter', '\t', 'precision', 16);
save(strcat(parent, '/05_xgboost/whole_sample/xgb_wm_pca_training.mat'), 'Train_Samples_means', 'Train_Samples_pca_coeff')

% exit
exit

# open R
R

# set libPath, load required packages, and set parent
.libPaths("/home/jawinskp/R/lib")
library(data.table)
library(xgboost)
parent = '/home/groups/markett/UK_Biobank/11b_brainage'

# load data (cat12-preprocessed T1, 'gm', 'wm', 'covs', and 'u')
load(paste0(parent, '/05_primary_dataset.RData'))

# set Target
i = 1;
j = 1;
Target = data.frame(covs[,5])
Train_Targets = data.frame(Target[u[,j]!=v[i,j],1])
Validation_Targets = data.frame(Target[u[,j]==v[i,j],1])

# -----------
# Grey matter 
# -----------

# load pca results
Train_Samples_pca = data.frame(fread(paste0(parent, '/05_xgboost/whole_sample/xgb_gm_Train_Samples_pca.txt'), sep = "\t", header = F))
Validation_Samples_pca = data.frame(fread(paste0(parent, '/05_xgboost/whole_sample/xgb_gm_Validation_Samples_pca.txt'), sep = "\t", header = F))

# run machine learning
dtrain = xgb.DMatrix(data = as.matrix(Train_Samples_pca), label=as.matrix(Train_Targets))
dvalid = xgb.DMatrix(data = as.matrix(Validation_Samples_pca), label=as.matrix(Validation_Targets))
        

watchlist = list(train=dtrain, test=dvalid)
model_tree = xgb.train(data = dtrain,
    watchlist = watchlist,
    params = list(booster = 'gbtree'),
    eval_metric = "mae",
    max.depth = 3,
    eta = 0.02,
    nthread = 28,
    nrounds = 10000,
    early_stopping_rounds = 50,
    verbose = 1,
    objective = "reg:linear")
        
model_linear = xgb.train(data = dtrain,
    watchlist = watchlist,
    params = list(booster = 'gblinear'),
    eval_metric = "mae",
    eta = 0.02,
    nthread = 28,
    nrounds = 10000,
    early_stopping_rounds = 50,
    verbose = 1,
    objective = "reg:linear")
        
# save models and predictions
save(model_tree, model_linear, file = paste0(parent, '/05_xgboost/whole_sample/xgb_gm_models.RData'), version = 2)

# ------------
# White matter
# ------------

# load pca results
Train_Samples_pca = data.frame(fread(paste0(parent, '/05_xgboost/whole_sample/xgb_wm_Train_Samples_pca.txt'), sep = "\t", header = F))
Validation_Samples_pca = data.frame(fread(paste0(parent, '/05_xgboost/whole_sample/xgb_wm_Validation_Samples_pca.txt'), sep = "\t", header = F))

# run machine learning
dtrain = xgb.DMatrix(data = as.matrix(Train_Samples_pca), label=as.matrix(Train_Targets))
dvalid = xgb.DMatrix(data = as.matrix(Validation_Samples_pca), label=as.matrix(Validation_Targets))
        
watchlist = list(train=dtrain, test=dvalid)
model_tree = xgb.train(data = dtrain,
    watchlist = watchlist,
    params = list(booster = 'gbtree'),
    eval_metric = "mae",
    max.depth = 3,
    eta = 0.02,
    nthread = 28,
    nrounds = 10000,
    early_stopping_rounds = 50,
    verbose = 1,
    objective = "reg:linear")
        
model_linear = xgb.train(data = dtrain,
    watchlist = watchlist,
    params = list(booster = 'gblinear'),
    eval_metric = "mae",
    eta = 0.02,
    nthread = 28,
    nrounds = 10000,
    early_stopping_rounds = 50,
    verbose = 1,
    objective = "reg:linear")
        
# save models and predictions
save(model_tree, model_linear, file = paste0(parent, '/05_xgboost/whole_sample/xgb_wm_models.RData'), version = 2)

# --------
# Clean up
# --------

# delete .m and .txt files
system(paste0('rm -f ', parent, '/05_xgboost/whole_sample/*pca.txt'))
        
# quit R and delete "05b_" files
quit('no')

# ----
# Sync
# ----

rsync -avP jawinskp@cluster1.psychologie.hu-berlin.de:/home/groups/markett/UK_Biobank/11b_brainage/05_xgboost/whole_sample/ /Volumes/psymol-1/PsyMol/projects/UK_Biobank/11b_brainage/05_xgboost/whole_sample/

ssh /Volumes/psymol-1/PsyMol/projects/UK_Biobank/11b_brainage/05_xgboost/whole_sample/xgboost_pred.R jawinskp@cluster1.psychologie.hu-berlin.de:/home/groups/markett/UK_Biobank/11b_brainage/05_xgboost/whole_sample/xgboost_pred.R

