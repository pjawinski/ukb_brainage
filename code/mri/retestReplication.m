% =======================================================================================================
% === collect and stack brain age estimates for repeat-imaging visit (cross-validation + singlemodel) ===
% =======================================================================================================
% /opt/matlab/bin/matlab -nodesktop -nodisplay -r "workingDir = pwd; threads = 100; spmPath = '/fast/software/matlab/spm12/'; rvmPath = '/fast/software/matlab/RVM/'; sbPath = '/fast/software/matlab/RVM/SB2_Release_200/'; run code/mri/brainageRetest.m"
fprintf('\n--- Collect & stack brain age estimates for repeat-imaging visit: discovery cohort ---\n')

% set working directory
cd(workingDir)

% addpath
addpath('code/functions/')
addpath(spmPath)
addpath(rvmPath)
addpath(sbPath)

% load dataset
fprintf('Loading dataset.\n')
replication = load('results/mri/brainage.replication.mat');
    gwm = load('results/mri/prepML.replication.mat','gm','wm');
    replication.gm = gwm.gm;
    replication.wm = gwm.wm;
clearvars gwm
load('results/mri/prepML.retest.mat');

% select replication dataset
member_logical = ismember(meta.IID, replication.meta.IID);
meta.IID = meta.IID(member_logical,:);
meta.ratings = meta.ratings(member_logical,:);
meta.files = meta.files(member_logical,:);
covs.table = covs.table(member_logical,:);
gm = gm(member_logical,:);
wm = wm(member_logical,:);

% select individualds with IQR < 3
idx = meta.ratings(:,9) < 3;
brainage.data = brainage.data(idx,:);
meta.ratings = meta.ratings(idx,:);
meta.IID = meta.IID(idx,:);
covs.table = covs.table(idx,:);

% ----------------------------------------------------------
% ------------------------ xgboost ------------------------
% transform features to pca scores and apply xgboost models
fprintf('Starting with predictions from cross-validation models.\n')

% start parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = threads;
saveProfile(myCluster);
parpool('local',threads);

% apply models
fprintf(' - applying xgb models.\n')
pred_xgb = zeros(size(gm,1),2,2,10,100);
parfor j = 1:100
    for i = 1:10
        %fprintf('Applying model %04d/1000.\n', j*10-10+i);

        % gm: define files 
        modelFile = sprintf('results/mri/ml.xgb/gm/models/xgmodel_%03d_%02d/models.RData',  j, i);
        pcaTrainingFile = sprintf('results/mri/ml.xgb/gm/models/xgmodel_%03d_%02d/pca_training.mat', j, i);
        dataFile = sprintf('results/mri/Test_Samples_%03d_%02d_pca.txt', j, i);
        outFile = sprintf('results/mri/Test_Samples_%03d_%02d_pred.txt', j, i);

        % gm: transform to pca scores
        Test_Samples = gm;
        x = parload(pcaTrainingFile);
        Test_Samples_centered = Test_Samples-repmat(x.Train_Samples_means,size(Test_Samples,1),1);
        Test_Samples_pca_score = Test_Samples_centered/transpose(x.Train_Samples_pca_coeff);
        dlmwrite(dataFile, Test_Samples_pca_score, 'delimiter', '\t', 'precision', 16);
        
        % gm: predict brain age
        system(sprintf('Rscript code/mri/ml.xgboost.pred.R %s %s %s', modelFile, dataFile, outFile));
        tmp_gm = importdata(outFile);

        % gm: clean up
        system(sprintf('rm -f %s', dataFile));
        system(sprintf('rm -f %s', outFile));

        % wm: define files 
        modelFile = sprintf('results/mri/ml.xgb/wm/models/xgmodel_%03d_%02d/models.RData',  j, i);
        pcaTrainingFile = sprintf('results/mri/ml.xgb/wm/models/xgmodel_%03d_%02d/pca_training.mat', j, i);
        dataFile = sprintf('results/mri/Test_Samples_%03d_%02d_pca.txt', j, i);
        outFile = sprintf('results/mri/Test_Samples_%03d_%02d_pred.txt', j, i);

        % wm: transform to pca scores
        Test_Samples = wm;
        x = parload(pcaTrainingFile);
        Test_Samples_centered = Test_Samples-repmat(x.Train_Samples_means,size(Test_Samples,1),1);
        Test_Samples_pca_score = Test_Samples_centered/transpose(x.Train_Samples_pca_coeff);
        dlmwrite(dataFile, Test_Samples_pca_score, 'delimiter', '\t', 'precision', 16);
        
        % wm: predict brain age
        system(sprintf('Rscript code/mri/ml.xgboost.pred.R %s %s %s', modelFile, dataFile, outFile));
        tmp_wm = importdata(outFile);

        % wm: clean up
        system(sprintf('rm -f %s', dataFile));
        system(sprintf('rm -f %s', outFile));

        % send temporary files to master
        tmp = [tmp_gm tmp_wm];
        tmp = reshape(tmp, [size(tmp,1),2,2]);
        pred_xgb(:,:,:,i,j) = tmp;
            
    end
end

if sum(pred_xgb(:)==0) == 0
    fprintf(' - xgb models successfully applied.\n')
else
    fprintf(' - missing predictions identified, exiting.\n')
    exit
end

% -------------------------------------------------------
% ------------------------- rvm -------------------------
% transform features to pca scores and apply rvm models

% apply models
fprintf(' - applying rvm models.\n')
pred_rvm = zeros(size(gm,1),2,10,100);
parfor j = 1:100
    for i = 1:10
        % fprintf('Applying model %04d/1000.\n', j*10-10+i);

        % gm
        Test_Samples = gm;
        x = parload(sprintf('results/mri/ml.rvm/gm/rvm_%03d_%02d.mat', j, i));
        Test_Samples_centered = Test_Samples-repmat(x.Train_Samples_means,size(Test_Samples,1),1);
        Test_Samples_pca_score = Test_Samples_centered/transpose(x.Train_Samples_pca_coeff);
        [tmp_gm, ~] = rvm_test(x.RVM,Test_Samples_pca_score);
        
        % wm
        Test_Samples = wm;
        x = parload(sprintf('results/mri/ml.rvm/wm/rvm_%03d_%02d.mat', j, i));
        Test_Samples_centered = Test_Samples-repmat(x.Train_Samples_means,size(Test_Samples,1),1);
        Test_Samples_pca_score = Test_Samples_centered/transpose(x.Train_Samples_pca_coeff);
        [tmp_wm, ~] = rvm_test(x.RVM,Test_Samples_pca_score);
        
        % send temporary files to master
        tmp = [tmp_gm tmp_wm];
        pred_rvm(:,:,i,j) = tmp;

    end
end

if sum(pred_rvm(:)==0) == 0
    fprintf(' - rvm models successfully applied.\n')
else
    fprintf(' - missing predictions identified, exiting.\n')
    exit
end

% --------------------
% --- stack models ---

% load stackmodels
fprintf(' - loading stackmodels.\n')
load('results/mri/brainage.stackmodels.mat')

% gm: get stacked brainage estimates
fprintf(' - stacking brain age estimates: grey matter.\n')
stack_predict_gm = NaN(size(pred_rvm,1),10,10,100);
for i = 1:100
    for j = 1:10
        predictions = NaN(size(pred_rvm,1),10);
        Test_Samples = [pred_xgb(:,:,1,j,i) pred_rvm(:,1,j,i)];
        for k = 1:10
            predictions(:,k) = predict(stackmodels.gm{k,j,i},Test_Samples);           
        end
        stack_predict_gm(:,:,j,i) = predictions;
    end
end

% wm: get stacked brainage estimates
fprintf(' - stacking brain age estimates: white matter.\n')
stack_predict_wm = NaN(size(pred_rvm,1),10,10,100);
for i = 1:100
    for j = 1:10
        predictions = NaN(size(pred_rvm,1),10);
        Test_Samples = [pred_xgb(:,:,2,j,i) pred_rvm(:,2,j,i)];
        for k = 1:10
            predictions(:,k) = predict(stackmodels.wm{k,j,i},Test_Samples);           
        end
        stack_predict_wm(:,:,j,i) = predictions;
    end
end

% gwm xgtree: get stacked brainage estimates
fprintf(' - stacking brain age estimates: grey and white matter (xgtree).\n')
stack_predict_gwm_xgtree = NaN(size(pred_rvm,1),10,10,100);
for i = 1:100
    for j = 1:10
        predictions = NaN(size(pred_rvm,1),10);
        Test_Samples = [pred_xgb(:,1,1,j,i) pred_xgb(:,1,2,j,i)];
        for k = 1:10
            predictions(:,k) = predict(stackmodels.gwm_xgtree{k,j,i},Test_Samples);           
        end
        stack_predict_gwm_xgtree(:,:,j,i) = predictions;
    end
end

% gwm xglin: get stacked brainage estimates
fprintf(' - stacking brain age estimates: grey and white matter (xglin).\n')
stack_predict_gwm_xglin = NaN(size(pred_rvm,1),10,10,100);
for i = 1:100
    for j = 1:10
        predictions = NaN(size(pred_rvm,1),10);
        Test_Samples = [pred_xgb(:,2,1,j,i) pred_xgb(:,2,2,j,i)];
        for k = 1:10
            predictions(:,k) = predict(stackmodels.gwm_xglin{k,j,i},Test_Samples);           
        end
        stack_predict_gwm_xglin(:,:,j,i) = predictions;
    end
end

% gwm rvm: get stacked brainage estimates
fprintf(' - stacking brain age estimates: grey and white matter (rvm).\n')
stack_predict_gwm_rvm = NaN(size(pred_rvm,1),10,10,100);
for i = 1:100
    for j = 1:10
        predictions = NaN(size(pred_rvm,1),10);
        Test_Samples = [pred_rvm(:,:,j,i)];
        for k = 1:10
            predictions(:,k) = predict(stackmodels.gwm_rvm{k,j,i},Test_Samples);           
        end
        stack_predict_gwm_rvm(:,:,j,i) = predictions;
    end
end

% gwm: get stacked brainage estimates
fprintf(' - stacking brain age estimates: all tissues and algorithms.\n')
stack_predict_gwm = NaN(size(pred_rvm,1),10,10,100);
for i = 1:100
    for j = 1:10
        predictions = NaN(size(pred_rvm,1),10);
        Test_Samples = [pred_xgb(:,:,1,j,i) pred_rvm(:,1,j,i) pred_xgb(:,:,2,j,i) pred_rvm(:,2,j,i)];
        for k = 1:10
            predictions(:,k) = predict(stackmodels.gwm{k,j,i},Test_Samples);           
        end
        stack_predict_gwm(:,:,j,i) = predictions;
    end
end

% prepare arrays for brain age calculation
pred_xgb = reshape(pred_xgb, [size(pred_xgb,1),2,2,1000]);
gm_xgtree_array = pred_xgb(:,1,1,:);
gm_xgtree_array = reshape(gm_xgtree_array, [size(gm_xgtree_array,1),1000]);
gm_xglin_array = pred_xgb(:,2,1,:);
gm_xglin_array = reshape(gm_xglin_array, [size(gm_xglin_array,1),1000]);
wm_xgtree_array = pred_xgb(:,1,2,:);
wm_xgtree_array = reshape(wm_xgtree_array, [size(wm_xgtree_array,1),1000]);
wm_xglin_array = pred_xgb(:,2,2,:);
wm_xglin_array = reshape(wm_xglin_array, [size(wm_xglin_array,1),1000]);

pred_rvm = reshape(pred_rvm, [size(pred_rvm,1),2,1000]);
gm_rvm_array = pred_rvm(:,1,:);
gm_rvm_array = reshape(gm_rvm_array, [size(gm_rvm_array,1),1000]);
wm_rvm_array = pred_rvm(:,2,:);
wm_rvm_array = reshape(wm_rvm_array, [size(wm_rvm_array,1),1000]);

% get brain age mean estimates
brainage_gm_xgtree = mean(gm_xgtree_array,2);
brainage_gm_xglin = mean(gm_xglin_array,2);
brainage_gm_rvm = mean(gm_rvm_array,2);
brainage_gm_stack = mean(reshape(stack_predict_gm, [size(stack_predict_gm,1),10000]),2);

brainage_wm_xgtree = mean(wm_xgtree_array,2);
brainage_wm_xglin = mean(wm_xglin_array,2);
brainage_wm_rvm = mean(wm_rvm_array,2);
brainage_wm_stack = mean(reshape(stack_predict_wm, [size(stack_predict_wm,1),10000]),2);

brainage_gwm_xgtree = mean(reshape(stack_predict_gwm_xgtree, [size(stack_predict_gwm_xgtree,1),10000]),2);
brainage_gwm_xglin = mean(reshape(stack_predict_gwm_xglin, [size(stack_predict_gwm_xglin,1),10000]),2);
brainage_gwm_rvm = mean(reshape(stack_predict_gwm_rvm, [size(stack_predict_gwm_rvm,1),10000]),2);
brainage_gwm_stack = mean(reshape(stack_predict_gwm, [size(stack_predict_gwm,1),10000]),2);

% create structure
brainage = struct();
brainage.data = [brainage_gm_xgtree brainage_gm_xglin brainage_gm_rvm brainage_gm_stack brainage_wm_xgtree brainage_wm_xglin brainage_wm_rvm brainage_wm_stack brainage_gwm_xgtree brainage_gwm_xglin brainage_gwm_rvm brainage_gwm_stack];
brainage.varnames = {'brainage_gm_xgtree','brainage_gm_xglin', 'brainage_gm_rvm', 'brainage_gm_stack', 'brainage_wm_xgtree', 'brainage_wm_xglin', 'brainage_wm_rvm', 'brainage_wm_stack', 'brainage_gwm_xgtree', 'brainage_gwm_xglin', 'brainage_gwm_rvm', 'brainage_gwm_stack', ...
    'brainage_gap_gm_xgtree','brainage_gap_gm_xglin', 'brainage_gap_gm_rvm', 'brainage_gap_gm_stack', 'brainage_gap_wm_xgtree', 'brainage_gap_wm_xglin', 'brainage_gap_wm_rvm', 'brainage_gap_wm_stack', 'brainage_gap_gwm_xgtree', 'brainage_gap_gwm_xglin', 'brainage_gap_gwm_rvm', 'brainage_gap_gwm_stack', ...
    'brainage_gap_adj_gm_xgtree','brainage_gap_adj_gm_xglin', 'brainage_gap_adj_gm_rvm', 'brainage_gap_adj_gm_stack', 'brainage_gap_adj_wm_xgtree', 'brainage_gap_adj_wm_xglin', 'brainage_gap_adj_wm_rvm', 'brainage_gap_adj_wm_stack', 'brainage_gap_adj_gwm_xgtree', 'brainage_gap_adj_gwm_xglin', 'brainage_gap_adj_gwm_rvm', 'brainage_gap_adj_gwm_stack'};

% calculate brainage gap
age = covs.table.t2_age;
for i = 1:12
    brainage.data(:,i+12) = brainage.data(:,i) - age;
end
   
% adjust brainage gap by sex, age, age2, TIV (% do not use cases with IQR_poor)
sex = covs.table.sex;
age2 = covs.table.t2_age2;
ac1 = covs.table.t2_ac1;
ac2 = covs.table.t2_ac2;
TIV = covs.table.t2_TIV;

for i = 1:12
    model = fitlm([sex age age2 ac1 ac2 TIV], brainage.data(:,i+12));
    brainage.data(:,i+24) = model.Residuals.Raw;
end

% compare mae and correlations of brainage with chronological age
fprintf(' - calculate correlations between chronological and estimated age.\n')
[corr(age,brainage.data(:,1:4));...
    corr(age,brainage.data(:,5:8));...
    corr(age,brainage.data(:,9:12))]
   
fprintf(' - calculate mean absolute errors between chronological and estimated age.\n')
[mean(abs(brainage.data(:,13:16)));...
   mean(abs(brainage.data(:,17:20)));...
   mean(abs(brainage.data(:,21:24)))]

% ==================================================================================
% === get age predictions from single models (derived from full discovery sample ===
% ==================================================================================
fprintf('Starting with predictions from single models (derived from full discovery sample).\n')

% ----------------------------------------------------------
% ------------------------- xgboost ------------------------
% transform features to pca scores and apply xgboost model

fprintf(' - applying xgb models.\n')
tissue = {'gm','wm'};
for i = 1:length(tissue)

    % define files 
    modelFile = sprintf('results/mri/ml.xgb/singlemodel/xgb_%s_models.RData', tissue{i});
    pcaTrainingFile = sprintf('results/mri/ml.xgb/singlemodel/xgb_%s_pca_training.mat', tissue{i});
    dataFile = sprintf('results/mri/Test_Samples_%s_pca.txt', tissue{i});
    outFile = sprintf('results/mri/Test_Samples_%s_pred.txt', tissue{i});

    % transform to pca scores and export
    eval(sprintf('Test_Samples = %s;', tissue{i}));
    x = load(pcaTrainingFile);
    Test_Samples_centered = Test_Samples-repmat(x.Train_Samples_means,size(Test_Samples,1),1);
    Test_Samples_pca_score = Test_Samples_centered/transpose(x.Train_Samples_pca_coeff);
    dlmwrite(dataFile, Test_Samples_pca_score, 'delimiter', '\t', 'precision', 16)

    % predict brain age
    system(sprintf('Rscript code/mri/ml.xgboost.pred.R %s %s %s', modelFile, dataFile, outFile));
    eval(sprintf('pred_%s = importdata(outFile);', tissue{i}));

    % clean up
    system(sprintf('rm -f %s', dataFile));
    system(sprintf('rm -f %s', outFile));
end

pred_xgb = [pred_gm pred_wm];
if sum(pred_xgb(:)==0) == 0
    fprintf(' - xgb models successfully applied.\n')
else
    fprintf(' - missing predictions identified, exiting.\n')
    exit
end

% -------------------------------------------------------
% ------------------------- rvm -------------------------
% transform features to pca scores and apply rvm models

fprintf(' - applying rvm models.\n')
tissue = {'gm','wm'};
for i = 1:length(tissue)

    % load model
    load(sprintf('results/mri/ml.rvm/singlemodel/rvm_%s.mat', tissue{i}));

    % get estimates
    eval(sprintf('Test_Samples = %s;', tissue{i}));
    Test_Samples_centered = Test_Samples-repmat(RVM.train_means,size(Test_Samples,1),1);
    Test_Samples_pca_score = Test_Samples_centered/RVM.train_pca_coeff';
    [y_mu, ~] = rvm_test(RVM,Test_Samples_pca_score);
    eval(sprintf('pred_%s = y_mu;', tissue{i}));
end

pred_rvm = [pred_gm pred_wm];
if sum(pred_rvm(:)==0) == 0
    fprintf(' - rvm models successfully applied.\n')
else
    fprintf(' - missing predictions identified, exiting.\n')
    exit
end

% --------------------------------
% --- merge into one structure ---

% create structure
brainage_singlemodel = struct();
brainage_singlemodel.data = NaN(size(gm,1),36);
brainage_singlemodel.varnames = {'brainage_gm_xgtree','brainage_gm_xglin', 'brainage_gm_rvm', 'brainage_gm_stack', 'brainage_wm_xgtree', 'brainage_wm_xglin', 'brainage_wm_rvm', 'brainage_wm_stack', 'brainage_gwm_xgtree', 'brainage_gwm_xglin', 'brainage_gwm_rvm', 'brainage_gwm_stack', ...
    'brainage_gap_gm_xgtree','brainage_gap_gm_xglin', 'brainage_gap_gm_rvm', 'brainage_gap_gm_stack', 'brainage_gap_wm_xgtree', 'brainage_gap_wm_xglin', 'brainage_gap_wm_rvm', 'brainage_gap_wm_stack', 'brainage_gap_gwm_xgtree', 'brainage_gap_gwm_xglin', 'brainage_gap_gwm_rvm', 'brainage_gap_gwm_stack', ...
    'brainage_gap_adj_gm_xgtree','brainage_gap_adj_gm_xglin', 'brainage_gap_adj_gm_rvm', 'brainage_gap_adj_gm_stack', 'brainage_gap_adj_wm_xgtree', 'brainage_gap_adj_wm_xglin', 'brainage_gap_adj_wm_rvm', 'brainage_gap_adj_wm_stack', 'brainage_gap_adj_gwm_xgtree', 'brainage_gap_adj_gwm_xglin', 'brainage_gap_adj_gwm_rvm', 'brainage_gap_adj_gwm_stack'};
brainage_singlemodel.data(:,[1 2 5 6]) = pred_xgb;
brainage_singlemodel.data(:,[3 7]) = pred_rvm;

% stack estimates through linear regression
brainage_singlemodel.data(:,4) = stack_coefficients.gm(1) + stack_coefficients.gm(2) * brainage_singlemodel.data(:,1) + stack_coefficients.gm(3) * brainage_singlemodel.data(:,2) + stack_coefficients.gm(4) * brainage_singlemodel.data(:,3);
brainage_singlemodel.data(:,8) = stack_coefficients.wm(1) + stack_coefficients.wm(2) * brainage_singlemodel.data(:,5) + stack_coefficients.wm(3) * brainage_singlemodel.data(:,6) + stack_coefficients.wm(4) * brainage_singlemodel.data(:,7);
brainage_singlemodel.data(:,9) = stack_coefficients.gwm_xgtree(1) + stack_coefficients.gwm_xgtree(2) * brainage_singlemodel.data(:,1) + stack_coefficients.gwm_xgtree(3) * brainage_singlemodel.data(:,5);
brainage_singlemodel.data(:,10) = stack_coefficients.gwm_xglin(1) + stack_coefficients.gwm_xglin(2) * brainage_singlemodel.data(:,2) + stack_coefficients.gwm_xglin(3) * brainage_singlemodel.data(:,6);
brainage_singlemodel.data(:,11) = stack_coefficients.gwm_rvm(1) + stack_coefficients.gwm_rvm(2) * brainage_singlemodel.data(:,3) + stack_coefficients.gwm_rvm(3) * brainage_singlemodel.data(:,7);
brainage_singlemodel.data(:,12) = stack_coefficients.gwm(1) + stack_coefficients.gwm(2) * brainage_singlemodel.data(:,1) + stack_coefficients.gwm(3) * brainage_singlemodel.data(:,2) + stack_coefficients.gwm(4) * brainage_singlemodel.data(:,3) + stack_coefficients.gwm(5) * brainage_singlemodel.data(:,5) + stack_coefficients.gwm(6) * brainage_singlemodel.data(:,6) + stack_coefficients.gwm(7) * brainage_singlemodel.data(:,7);

% calculate gap
age = covs.table.t2_age;
for i = 1:12
    brainage_singlemodel.data(:,i+12) = brainage_singlemodel.data(:,i) - age;
end

% adjust brainage gap by sex, age, age2, assessment center, TIV
sex = covs.table.sex;
age2 = covs.table.t2_age2;
ac1 = covs.table.t2_ac1;
ac2 = covs.table.t2_ac2;
TIV = covs.table.t2_TIV;
for i = 1:12
    model = fitlm([sex age age2 ac1 ac2 TIV], brainage_singlemodel.data(:,i+12));
    brainage_singlemodel.data(:,i+24) = model.Residuals.Raw;
end

% compare mae and correlations of brainage with chronological age
fprintf(' - calculate correlations between chronological and estimated age.\n')
[corr(age,brainage_singlemodel.data(:,1:4));...
    corr(age,brainage_singlemodel.data(:,5:8));...
    corr(age,brainage_singlemodel.data(:,9:12))]

fprintf(' - calculate mean absolute errors between chronological and estimated age.\n')
[mean(abs(brainage_singlemodel.data(:,13:16)));...
   mean(abs(brainage_singlemodel.data(:,17:20)));...
   mean(abs(brainage_singlemodel.data(:,21:24)))]    

fprintf(' - correlating brain age gap from 10000 models with singlemodel estimates.\n')
diag(corr(brainage_singlemodel.data(:,13:24), brainage.data(:,13:24)))'

% save predictions
fprintf('- saving brainage predictions.\n')
save('results/mri/brainage.replication.retest.mat', 'brainage', 'brainage_singlemodel', 'covs', 'gm_logical', 'meta', 'wm_logical')
colNames = [meta.ratings_varnames, brainage.varnames];
brainage_table = [covs.table array2table([meta.ratings, brainage.data],'VariableNames',colNames)];
writetable(brainage_table, 'results/mri/brainage.replication.retest.txt', 'Delimiter', '\t', 'WriteRowNames', 0)
brainage_table = [covs.table array2table([meta.ratings, brainage_singlemodel.data],'VariableNames',colNames)];
writetable(brainage_table, 'results/mri/brainage.replication.retest.singlemodel.txt', 'Delimiter', '\t', 'WriteRowNames', 0)

% quit matlab
fprintf('\n--- Completed: Collect & stack brain age estimates for repeat-imaging visit: replication cohort ---\n')
exit