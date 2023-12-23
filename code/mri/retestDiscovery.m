% =======================================================================================================
% === collect and stack brain age estimates for repeat-imaging visit (cross-validation + singlemodel) ===
% =======================================================================================================
% /opt/matlab/bin/matlab -nodesktop -nodisplay -r "workingDir = pwd; threads = 100; spmPath = '/fast/software/matlab/spm12/'; rvmPath = '/fast/software/matlab/RVM/'; sbPath = '/fast/software/matlab/RVM/SB2_Release_200/'; run code/mri/retestDiscovery.m"
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
discovery = load('results/mri/brainage.discovery.mat');
    gwm = load('results/mri/prepML.discovery.mat','gm','wm');
    discovery.gm = gwm.gm;
    discovery.wm = gwm.wm;
clearvars gwm
load('results/mri/prepML.retest.mat');

% align with discovery dataset - white-British subjects only!
member_logical = ismember(meta.IID, discovery.meta.IID);
meta.IID = meta.IID(member_logical,:);
meta.ratings = meta.ratings(member_logical,:);
meta.files = meta.files(member_logical,:);
covs.table = covs.table(member_logical,:);
gm = gm(member_logical,:);
wm = wm(member_logical,:);

% get indices of subjects in discovery dataset
indices = NaN(size(meta.IID));
for i = 1:length(indices)
    id = meta.IID(i,1);
    for j = 1:length(discovery.meta.IID)
        if discovery.meta.IID(j) == id
            indices(i,1) = j;
            break;
        end
    end
end

% put into matrices identical in size of discovery dataset
meta_ratings = zeros(size(discovery.meta.ratings));
meta_ratings(indices,:) = meta.ratings;
meta_IID = zeros(size(discovery.meta.IID));
meta_IID(indices,:) = meta.IID; 
gm_matrix = zeros(size(discovery.gm));
wm_matrix = zeros(size(discovery.wm));
gm_matrix(indices,:) = gm;
wm_matrix(indices,:) = wm;
covs_table = discovery.covs.table;

% put into structures
covs.table = covs_table;
meta.ratings = meta_ratings;
meta.IID = meta_IID;
gm = gm_matrix;
wm = wm_matrix;

% ----------------------------------------------------------
% ------------------------- xgboost ------------------------
% transform features to pca scores and apply xgboost models
fprintf('Starting with predictions from cross-validation models.\n')

% start parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = threads;
saveProfile(myCluster);
parpool('local',threads);

% apply models to discovery test samples (now containing retest-data)
fprintf(' - applying rvm models.\n')

pred_xgb_cell = {};
tic
parfor j = 1:100
    for i = 1:10

        % echo current status
        % fprintf('Applying model %04d/1000.\n', j*10-10+i);

        % gm: define files 
        modelFile = sprintf('results/mri/ml.xgb/gm/models/xgmodel_%03d_%02d/models.RData',  j, i);
        pcaTrainingFile = sprintf('results/mri/ml.xgb/gm/models/xgmodel_%03d_%02d/pca_training.mat', j, i);
        dataFile = sprintf('results/mri/Test_Samples_gm_%03d_%02d_pca.txt', j, i);
        outFile = sprintf('results/mri/Test_Samples_gm_%03d_%02d_pred.txt', j, i);

        % gm: transform to pca scores
        Test_Samples = gm(discovery.u(:,j)==i,:);
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
        dataFile = sprintf('results/mri/Test_Samples_wm_%03d_%02d_pca.txt', j, i);
        outFile = sprintf('results/mri/Test_Samples_wm_%03d_%02d_pred.txt', j, i);

        % wm: transform to pca scores
        Test_Samples = wm(discovery.u(:,j)==i,:);
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
        pred_xgb_cell{i,j} = tmp;

    end
end
toc

% get estimates
pred_xgb = zeros(size(gm,1),2,2,100);
for j = 1:100
    for i = 1:10
        pred_xgb(discovery.u(:,j)==i,:,:,j) = pred_xgb_cell{i,j};
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

% apply models (takes ~20 mins with 100 cores)
fprintf(' - applying rvm models.\n')

pred_rvm_cell = {};
parfor j = 1:100
    for i = 1:10   
        % fprintf('Applying model %04d/1000.\n', j*10-10+i);

        % gm
        Test_Samples = gm(discovery.u(:,j)==i,:);
        x = parload(sprintf('results/mri/ml.rvm/gm/rvm_%03d_%02d.mat', j, i));
        Test_Samples_centered = Test_Samples-repmat(x.Train_Samples_means,size(Test_Samples,1),1);
        Test_Samples_pca_score = Test_Samples_centered/transpose(x.Train_Samples_pca_coeff);
        [tmp_gm, ~] = rvm_test(x.RVM,Test_Samples_pca_score);
        
        % wm
        Test_Samples = wm(discovery.u(:,j)==i,:);
        x = parload(sprintf('results/mri/ml.rvm/wm/rvm_%03d_%02d.mat', j, i));
        Test_Samples_centered = Test_Samples-repmat(x.Train_Samples_means,size(Test_Samples,1),1);
        Test_Samples_pca_score = Test_Samples_centered/transpose(x.Train_Samples_pca_coeff);
        [tmp_wm, ~] = rvm_test(x.RVM,Test_Samples_pca_score);
        
        % send temporary files to master
        tmp = [tmp_gm tmp_wm];
        pred_rvm_cell{i,j} = tmp;

    end
end

% get estimates
pred_rvm = zeros(size(gm,1),2,100);
for j = 1:100
    for i = 1:10
        pred_rvm(discovery.u(:,j)==i,:,j) = pred_rvm_cell{i,j};
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
stack_predict_gm = NaN(size(pred_rvm,1),100);
for i = 1:100
    for j = 1:10
        predictions = NaN(sum(discovery.u(:,i)==j),1);
        for k = 1:10
            Test_Samples = [pred_xgb(discovery.u(:,i)==j,:,1,i) pred_rvm(discovery.u(:,i)==j,1,i)];
            Test_Samples = Test_Samples(discovery.w{j,i}==k,:);
            predictions(discovery.w{j,i}==k) = predict(stackmodels.gm{k,j,i},Test_Samples);           
        end
        stack_predict_gm(discovery.u(:,i)==j,i) = predictions;
    end
end

% wm: get stacked brainage estimates
fprintf(' - stacking brain age estimates: white matter.\n')
stack_predict_wm = NaN(size(pred_rvm,1),100);
for i = 1:100
    for j = 1:10
        predictions = NaN(sum(discovery.u(:,i)==j),1);
        for k = 1:10
            Test_Samples = [pred_xgb(discovery.u(:,i)==j,:,2,i) pred_rvm(discovery.u(:,i)==j,2,i)];
            Test_Samples = Test_Samples(discovery.w{j,i}==k,:);
            predictions(discovery.w{j,i}==k) = predict(stackmodels.wm{k,j,i},Test_Samples);           
        end
        stack_predict_wm(discovery.u(:,i)==j,i) = predictions;
    end
end

% gwm xgtree: get stacked brainage estimates
fprintf(' - stacking brain age estimates: grey and white matter (xgtree).\n')
stack_predict_gwm_xgtree = NaN(size(pred_rvm,1),100);
for i = 1:100
    for j = 1:10
        predictions = NaN(sum(discovery.u(:,i)==j),1);
        for k = 1:10
            Test_Samples = [pred_xgb(discovery.u(:,i)==j,1,1,i) pred_xgb(discovery.u(:,i)==j,1,2,i)];
            Test_Samples = Test_Samples(discovery.w{j,i}==k,:);
            predictions(discovery.w{j,i}==k) = predict(stackmodels.gwm_xgtree{k,j,i},Test_Samples);           
        end
        stack_predict_gwm_xgtree(discovery.u(:,i)==j,i) = predictions;
    end
end

% gwm xglin: get stacked brainage estimates
fprintf(' - stacking brain age estimates: grey and white matter (xglin).\n')
stack_predict_gwm_xglin = NaN(size(pred_rvm,1),100);
for i = 1:100
    for j = 1:10
        predictions = NaN(sum(discovery.u(:,i)==j),1);
        for k = 1:10
            Test_Samples = [pred_xgb(discovery.u(:,i)==j,2,1,i) pred_xgb(discovery.u(:,i)==j,2,2,i)];
            Test_Samples = Test_Samples(discovery.w{j,i}==k,:);
            predictions(discovery.w{j,i}==k) = predict(stackmodels.gwm_xglin{k,j,i},Test_Samples);           
        end
        stack_predict_gwm_xglin(discovery.u(:,i)==j,i) = predictions;
    end
end

% gwm rvm: get stacked brainage estimates
fprintf(' - stacking brain age estimates: grey and white matter (rvm).\n')
stack_predict_gwm_rvm = NaN(size(pred_rvm,1),100);
for i = 1:100
    for j = 1:10
        predictions = NaN(sum(discovery.u(:,i)==j),1);
        for k = 1:10
            Test_Samples = [pred_rvm(discovery.u(:,i)==j,1,i) pred_rvm(discovery.u(:,i)==j,2,i)];
            Test_Samples = Test_Samples(discovery.w{j,i}==k,:);
            predictions(discovery.w{j,i}==k) = predict(stackmodels.gwm_rvm{k,j,i},Test_Samples);           
        end
        stack_predict_gwm_rvm(discovery.u(:,i)==j,i) = predictions;
    end
end

% gwm: get stacked brainage estimates
fprintf(' - stacking brain age estimates: all tissues and algorithms.\n')
stack_predict_gwm = NaN(size(pred_rvm,1),100);
for i = 1:100
    for j = 1:10
        predictions = NaN(sum(discovery.u(:,i)==j),1);
        for k = 1:10
            Test_Samples = [pred_xgb(discovery.u(:,i)==j,:,1,i) pred_rvm(discovery.u(:,i)==j,1,i) pred_xgb(discovery.u(:,i)==j,:,2,i) pred_rvm(discovery.u(:,i)==j,2,i)];
            Test_Samples = Test_Samples(discovery.w{j,i}==k,:);
            predictions(discovery.w{j,i}==k) = predict(stackmodels.gwm{k,j,i},Test_Samples);           
        end
        stack_predict_gwm(discovery.u(:,i)==j,i) = predictions;
    end
end

% prepare arrays for brain age mean calculation
gm_xgtree_array = pred_xgb(:,1,1,:);
gm_xgtree_array = reshape(gm_xgtree_array, [size(gm_xgtree_array,1),100]);
gm_xglin_array = pred_xgb(:,2,1,:);
gm_xglin_array = reshape(gm_xglin_array, [size(gm_xglin_array,1),100]);
wm_xgtree_array = pred_xgb(:,1,2,:);
wm_xgtree_array = reshape(wm_xgtree_array, [size(wm_xgtree_array,1),100]);
wm_xglin_array = pred_xgb(:,2,2,:);
wm_xglin_array = reshape(wm_xglin_array, [size(wm_xglin_array,1),100]);

gm_rvm_array = pred_rvm(:,1,:);
gm_rvm_array = reshape(gm_rvm_array, [size(gm_rvm_array,1),100]);
wm_rvm_array = pred_rvm(:,2,:);
wm_rvm_array = reshape(wm_rvm_array, [size(wm_rvm_array,1),100]);

% get brain age mean estimates
brainage_gm_xgtree = mean(gm_xgtree_array,2);
brainage_gm_xglin = mean(gm_xglin_array,2);
brainage_gm_rvm = mean(gm_rvm_array,2);
brainage_gm_stack = mean(stack_predict_gm,2);

brainage_wm_xgtree = mean(wm_xgtree_array,2);
brainage_wm_xglin = mean(wm_xglin_array,2);
brainage_wm_rvm = mean(wm_rvm_array,2);
brainage_wm_stack = mean(stack_predict_wm,2);

brainage_gwm_xgtree = mean(stack_predict_gwm_xgtree,2);
brainage_gwm_xglin = mean(stack_predict_gwm_xglin,2);
brainage_gwm_rvm = mean(stack_predict_gwm_rvm,2);
brainage_gwm_stack = mean(stack_predict_gwm,2);

brainage = struct();
brainage.data = [brainage_gm_xgtree brainage_gm_xglin brainage_gm_rvm brainage_gm_stack brainage_wm_xgtree brainage_wm_xglin brainage_wm_rvm brainage_wm_stack brainage_gwm_xgtree brainage_gwm_xglin brainage_gwm_rvm brainage_gwm_stack];
brainage.varnames = {'brainage_gm_xgtree','brainage_gm_xglin', 'brainage_gm_rvm', 'brainage_gm_stack', 'brainage_wm_xgtree', 'brainage_wm_xglin', 'brainage_wm_rvm', 'brainage_wm_stack', 'brainage_gwm_xgtree', 'brainage_gwm_xglin', 'brainage_gwm_rvm', 'brainage_gwm_stack', ...
    'brainage_gap_gm_xgtree','brainage_gap_gm_xglin', 'brainage_gap_gm_rvm', 'brainage_gap_gm_stack', 'brainage_gap_wm_xgtree', 'brainage_gap_wm_xglin', 'brainage_gap_wm_rvm', 'brainage_gap_wm_stack', 'brainage_gap_gwm_xgtree', 'brainage_gap_gwm_xglin', 'brainage_gap_gwm_rvm', 'brainage_gap_gwm_stack', ...
    'brainage_gap_adj_gm_xgtree','brainage_gap_adj_gm_xglin', 'brainage_gap_adj_gm_rvm', 'brainage_gap_adj_gm_stack', 'brainage_gap_adj_wm_xgtree', 'brainage_gap_adj_wm_xglin', 'brainage_gap_adj_wm_rvm', 'brainage_gap_adj_wm_stack', 'brainage_gap_adj_gwm_xgtree', 'brainage_gap_adj_gwm_xglin', 'brainage_gap_adj_gwm_rvm', 'brainage_gap_adj_gwm_stack'};

% now extract retest data
brainage.data = brainage.data(indices,:);
meta.ratings = meta.ratings(indices,:);
meta.IID = meta.IID(indices,:);
covs.table = covs.table(indices,:);

% select individualds with IQR < 3
idx = meta.ratings(:,9) < 3;
brainage.data = brainage.data(idx,:);
meta.ratings = meta.ratings(idx,:);
meta.IID = meta.IID(idx,:);
covs.table = covs.table(idx,:);

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

% save predictions
fprintf('- saving brainage predictions.\n')
save('results/mri/brainage.discovery.retest.mat', 'brainage', 'covs', 'gm_logical', 'meta', 'wm_logical')
colNames = [meta.ratings_varnames, brainage.varnames];
brainage_table = [covs.table array2table([meta.ratings, brainage.data],'VariableNames',colNames)];
writetable(brainage_table, 'results/mri/brainage.discovery.retest.txt', 'Delimiter', '\t', 'WriteRowNames', 0)
fprintf('\n--- Completed: Collect & stack brain age estimates for repeat-imaging visit: discovery cohort ---\n')

