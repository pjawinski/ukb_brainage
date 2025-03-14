#!/bin/bash

# =====================================================================================
# === cat12 data collection, sample filtering, and preparation for machine learning ===
# =====================================================================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage
conda activate envs/default

# r2020: create dataset of cat12-preprocessed MRI scans released until Feb 2020 (discovery and some replication individuals)
/opt/matlab/bin/matlab -nodesktop -nodisplay -r "scanDir = {'data/t1w/r2020/'}; outfile = 'results/mri/cat12.r2020'; spmPath = '/fast/software/matlab/spm12/'; ncores = 50; run code/mri/cat12collect.m"

# r2024: create dataset of cat12-preprocessed MRI scans released until Feb 2024 (discovery and replication)
/opt/matlab/bin/matlab -nodesktop -nodisplay -r "scanDir = {'data/t1w/r2020/','data/t1w/r2021/','data/t1w/r2024/'}; outfile = 'results/mri/cat12.r2024'; spmPath = '/fast/software/matlab/spm12/'; ncores = 50; run code/mri/cat12collect.m"

# r2024.retest: create dataset with all cat12-preprocessed retest MRI scans
/opt/matlab/bin/matlab -nodesktop -nodisplay -r "scanDir = {'data/t1w/r2022_retest/','data/t1w/r2024_retest/'}; outfile = 'results/mri/cat12.r2024.retest'; spmPath = '/fast/software/matlab/spm12/'; ncores = 50; run code/mri/cat12collect.m"

# run sample filtering (divide into discovery and replication, select unrelated individuals)
Rscript code/mri/sampleFiltering.R

# prepare datasets with filtered samples for machine learning in MATLAB and R
/opt/matlab/bin/matlab -nodesktop -nodisplay -r "workingDir = pwd; run code/mri/prepML.m"
Rscript code/mri/prepML.R

# run machine learning using rvm: kfold-cross validation + single model (whole sample)
matFile="results/mri/prepML.discovery.mat"
spmPath="/fast/software/matlab/spm12/"
rvmPath="/fast/software/matlab/RVM/"
sbPath="/fast/software/matlab/RVM/SB2_Release_200/"
threads=100
for tissue in gm wm; do
    targetDir="results/mri/ml.rvm/${tissue}/"
    /opt/matlab/bin/matlab -nodesktop -nodisplay -r "tissue = '${tissue}'; workingDir = pwd; targetDir = '${targetDir}'; matFile = '${matFile}';  maxNumCompThreads(${threads}); spmPath = '${spmPath}'; rvmPath = '${rvmPath}'; sbPath = '${sbPath}'; run code/mri/ml_rvm.m"
    targetDir="results/mri/ml.rvm/singlemodel/"
    /opt/matlab/bin/matlab -nodesktop -nodisplay -r "tissue = '${tissue}'; workingDir = pwd; targetDir = '${targetDir}'; matFile = '${matFile}';  maxNumCompThreads(${threads}); spmPath = '${spmPath}'; rvmPath = '${rvmPath}'; sbPath = '${sbPath}'; run code/mri/ml_rvm_singlemodel.m"
done

# run machine learning using xgboost: kfold-cross validation + single model (whole sample)
conda activate envs/xgb
Rfile="results/mri/prepML.discovery.RData"
matFile="results/mri/prepML.discovery.mat"
matlabpath="/opt/matlab/bin/matlab"
threads=50
for tissue in gm wm; do
    targetDir="results/mri/ml.xgb/${tissue}/"
    Rscript code/mri/ml.xgboost.R "${tissue}" "${targetDir}" "${Rfile}" "${matFile}" "${matlabpath}" "${threads}"
    targetDir="results/mri/ml.xgb/singlemodel/"
    Rscript code/mri/ml.xgboost.singlemodel.R "${tissue}" "${targetDir}" "${Rfile}" "${matFile}" "${matlabpath}" "${threads}"
done

# benchmark single-model (single model from full discovery) vs. multi-model approach (models from cross-validation)
/opt/matlab/bin/matlab -nodesktop -nodisplay -r "workingDir = pwd; run code/mri/singlemodel_benchmark.m"

# extract feature importances for XGBoost models with tree booster
conda activate envs/xgb
width=7.7
height=3.0
for tissue in gm wm; do
    modelFile="results/mri/ml.xgb/singlemodel/xgb_${tissue}_models.RData"
    outFile="results/mri/ml.xgb/singlemodel/xgb_${tissue}_importance.png"
    ./code/mri/singlemodel.importance.R ${modelFile} ${outFile} ${width} ${height}
done

# merge single-model plots
conda activate envs/default
benchmark="results/mri/singlemodel.benchmark.png"
gmImportance="results/mri/ml.xgb/singlemodel/xgb_gm_importance.png"
wmImportance="results/mri/ml.xgb/singlemodel/xgb_wm_importance.png"
outFile="results/combined/singlemodel.composite.png"
./code/mri/singlemodel.combine.R "${benchmark}" "${gmImportance}" "${wmImportance}" "${outFile}"

# discovery cohort: collect and stack brain age estimates
/opt/matlab/bin/matlab -nodesktop -nodisplay -r "workingDir = pwd; run code/mri/brainageDiscovery.m"

# replication cohort & retest measures: apply models and collect brain age estimates
conda activate envs/xgb
spmPath="/fast/software/matlab/spm12/"
rvmPath="/fast/software/matlab/RVM/"
sbPath="/fast/software/matlab/RVM/SB2_Release_200/"
threads=100
/opt/matlab/bin/matlab -nodesktop -nodisplay -r "workingDir = pwd; threads = ${threads}; spmPath = '${spmPath}'; rvmPath = '${rvmPath}'; sbPath = '${sbPath}'; run code/mri/brainageReplication.m"
/opt/matlab/bin/matlab -nodesktop -nodisplay -r "workingDir = pwd; threads = ${threads}; spmPath = '${spmPath}'; rvmPath = '${rvmPath}'; sbPath = '${sbPath}'; run code/mri/retestDiscovery.m"
/opt/matlab/bin/matlab -nodesktop -nodisplay -r "workingDir = pwd; threads = ${threads}; spmPath = '${spmPath}'; rvmPath = '${rvmPath}'; sbPath = '${sbPath}'; run code/mri/retestReplication.m"

# LIFE-Adult replication sample: collect brain age estimates (models applied on-site).
/opt/matlab/bin/matlab -nodesktop -nodisplay -r "workingDir = pwd; run code/mri/brainageLIFE.m"

# get accuracy metrics and draw plot
/opt/matlab/bin/matlab -nodesktop -nodisplay -r "workingDir = pwd; run code/mri/getAccuracy.m"

# get Freesurfer variables 
conda activate envs/default
surfacePath="data/surface"
outFile="results/mri/freesurfer.tsv"
./code/mri/getFreesurfer.R "${surfacePath}" "${outFile}"

