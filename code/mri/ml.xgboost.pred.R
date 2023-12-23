#!/usr/bin/env Rscript
# ============================
# === Apply xgboost models === 
# ============================

# get arguments from command line
args = commandArgs(trailingOnly=TRUE)
modelFile = args[1] # modelFile="results/mri/ml.xgb/gm/xgb_001_01/models.RData"
dataFile = args[2] # dataFile="results/mri/ml.xgb/gm/Test_samples_001_01_pca.txt"
outFile = args[3] # outFile="results/mri/ml.xgb/gm/Test_samples_001_01_pred.txt"

# load required package
library(xgboost)

# load model and data
load(modelFile)
Test_Samples_pca = read.table(dataFile, sep='\t')

# predict
dtest = xgb.DMatrix(data = as.matrix(Test_Samples_pca), label=as.matrix(c(1:length(Test_Samples_pca[,1]))))
pred = data.frame(cbind(predict(model_tree, dtest), predict(model_linear, dtest)))

# save models and predictions
write.table(pred, file = outFile, sep = '\t', row.names = F, col.names = F)
