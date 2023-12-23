#!/bin/bash

# Contact: Philippe Jawinski (philippe.jawinski@hu-berlin.de)

# =============================================
# === create R library with package xgboost ===
# =============================================

# access hcp
ssh comps06h04
tmux
cd /data/p_02330

# create folder for Rlib
mkdir Rlib

# open R
R
.libPaths("/data/p_02330/Rlib")

# install xgboost
require(devtools)
install_version("xgboost", version = "0.82.1", repos = "http://cran.us.r-project.org")

# exit R
quit('no')

# start Matlab and load cat12 data required for brain age estimates
matlab -nodesktop -nodisplay 

% addpath
addpath /data/p_02330/functions/
addpath /data/p_02330/toolbox/spm12/
addpath /data/p_02330/toolbox/RVM/
addpath /data/p_02330/toolbox/RVM/SB2_Release_200

% set parent
parent = '/data/p_02330';

% load dataset
load(strcat(parent,'/01_cat12_data.mat'));

% -----------
% --- rvm ---
% -----------

% load models
load(strcat(parent, '/models/RVM_gm.mat'))
load(strcat(parent, '/models/RVM_wm.mat'))

% select features
gm = gm(:,RVM_gm.train_logical);
wm = wm(:,RVM_wm.train_logical);

% get rvm_gm estimates
RVM_gm.rv_index = 1:size(RVM_gm.rv_mu,1);
Test_Samples = gm;
Test_Samples_centered = Test_Samples-repmat(RVM_gm.train_means,size(Test_Samples,1),1);
Test_Samples_pca_score = Test_Samples_centered/RVM_gm.train_pca_coeff'; %'
[gm_y_mu, ~] = rvm_test(RVM_gm,Test_Samples_pca_score);

% get rvm_wm estimates
RVM_wm.rv_index = 1:size(RVM_wm.rv_mu,1);
Test_Samples = wm;
Test_Samples_centered = Test_Samples-repmat(RVM_wm.train_means,size(Test_Samples,1),1);
Test_Samples_pca_score = Test_Samples_centered/RVM_wm.train_pca_coeff'; %'
[wm_y_mu, ~] = rvm_test(RVM_wm,Test_Samples_pca_score);
pred_rvm = [gm_y_mu wm_y_mu];

% ---------------
% --- xgboost --- 
% ---------------

% gm: transform to pca scores
Test_Samples = gm;
x = load(strcat(parent, '/models/xgb_gm_pca_training.mat'));
Test_Samples_centered = Test_Samples-repmat(x.Train_Samples_means,size(Test_Samples,1),1);
Test_Samples_pca_score = Test_Samples_centered/transpose(x.Train_Samples_pca_coeff);
dlmwrite(strcat(parent, '/models/xgb_gm_Test_Samples_pca.txt'), Test_Samples_pca_score, 'delimiter', '\t', 'precision', 16)
 
% wm: transform to pca scores
Test_Samples = wm;
x = load(strcat(parent, '/models/xgb_wm_pca_training.mat'));
Test_Samples_centered = Test_Samples-repmat(x.Train_Samples_means,size(Test_Samples,1),1);
Test_Samples_pca_score = Test_Samples_centered/transpose(x.Train_Samples_pca_coeff);
dlmwrite(strcat(parent, '/models/xgb_wm_Test_Samples_pca.txt'), Test_Samples_pca_score, 'delimiter', '\t', 'precision', 16)

% predict in R
system(char(strcat({'Rscript '}, parent, {'/models/xgboost_pred.R '}, parent, '/models')));

% get estimates
pred_xgb = importdata(char(strcat(parent, '/models/xgb_pred.txt')));

% clean up
delete(strcat(parent, '/models/*Test_Samples*.txt'))
delete(strcat(parent, '/models/xgb_pred.txt'))

% ------------------------------
% --- put into one structure ---
% ------------------------------

% create structure
brainage = struct();
brainage.data = NaN(size(gm,1),36);
brainage.varnames = {'brainage_gm_xgtree','brainage_gm_xglin', 'brainage_gm_rvm', 'brainage_gm_stack', 'brainage_wm_xgtree', 'brainage_wm_xglin', 'brainage_wm_rvm', 'brainage_wm_stack', 'brainage_gwm_xgtree', 'brainage_gwm_xglin', 'brainage_gwm_rvm', 'brainage_gwm_stack', ...
    'brainage_gap_gm_xgtree','brainage_gap_gm_xglin', 'brainage_gap_gm_rvm', 'brainage_gap_gm_stack', 'brainage_gap_wm_xgtree', 'brainage_gap_wm_xglin', 'brainage_gap_wm_rvm', 'brainage_gap_wm_stack', 'brainage_gap_gwm_xgtree', 'brainage_gap_gwm_xglin', 'brainage_gap_gwm_rvm', 'brainage_gap_gwm_stack', ...
    'brainage_gap_adj_gm_xgtree','brainage_gap_adj_gm_xglin', 'brainage_gap_adj_gm_rvm', 'brainage_gap_adj_gm_stack', 'brainage_gap_adj_wm_xgtree', 'brainage_gap_adj_wm_xglin', 'brainage_gap_adj_wm_rvm', 'brainage_gap_adj_wm_stack', 'brainage_gap_adj_gwm_xgtree', 'brainage_gap_adj_gwm_xglin', 'brainage_gap_adj_gwm_rvm', 'brainage_gap_adj_gwm_stack'};

brainage.data(:,[1 2 5 6]) = pred_xgb;
brainage.data(:,[3 7]) = pred_rvm;

% stack estimates
load(strcat(parent, '/models/stack_coefficients.mat'))
brainage.data(:,4) = stack_coefficients.gm(1) + stack_coefficients.gm(2) * brainage.data(:,1) + stack_coefficients.gm(3) * brainage.data(:,2) + stack_coefficients.gm(4) * brainage.data(:,3);
brainage.data(:,8) = stack_coefficients.wm(1) + stack_coefficients.wm(2) * brainage.data(:,5) + stack_coefficients.wm(3) * brainage.data(:,6) + stack_coefficients.wm(4) * brainage.data(:,7);
brainage.data(:,9) = stack_coefficients.gwm_xgtree(1) + stack_coefficients.gwm_xgtree(2) * brainage.data(:,1) + stack_coefficients.gwm_xgtree(3) * brainage.data(:,5);
brainage.data(:,10) = stack_coefficients.gwm_xglin(1) + stack_coefficients.gwm_xglin(2) * brainage.data(:,2) + stack_coefficients.gwm_xglin(3) * brainage.data(:,6);
brainage.data(:,11) = stack_coefficients.gwm_rvm(1) + stack_coefficients.gwm_rvm(2) * brainage.data(:,3) + stack_coefficients.gwm_rvm(3) * brainage.data(:,7);
brainage.data(:,12) = stack_coefficients.gwm(1) + stack_coefficients.gwm(2) * brainage.data(:,1) + stack_coefficients.gwm(3) * brainage.data(:,2) + stack_coefficients.gwm(4) * brainage.data(:,3) + stack_coefficients.gwm(5) * brainage.data(:,5) + stack_coefficients.gwm(6) * brainage.data(:,6) + stack_coefficients.gwm(7) * brainage.data(:,7);

% ------------------
% --- make table ---
% ------------------

% get name of nifti and date
niftiname  = cell(size(meta.files));
niftidate  = cell(size(meta.files));
for i = 1:size(meta.files,1)
	niftiname{i,1} = strcat(meta.files(i).folder(25:end), '/', meta.files(i).name(1:(end-10)), '.nii.gz');
	niftidate{i,1} = meta.files(i).folder(47:54);
end

% merge variables
merged = [meta.IID, niftiname, niftidate, num2cell(meta.ratings), num2cell(brainage.data(:,1:12))];

% convert to table
colNames = [{'IID', 'NIFTI', 'DATE', 'TIV', 'CSF', 'GM', 'WM', 'VRES', 'RMS','NCR','ICR', 'IQR', 'IQR_poor'}, brainage.varnames(:,1:12)];
output = array2table(merged,'VariableNames',colNames);

% write table
writetable(output, strcat(parent, '/02_estimates_stacked.txt'), 'Delimiter', '\t', 'WriteRowNames', 0)

