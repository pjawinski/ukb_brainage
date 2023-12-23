% ==================================================================
% === collect and stack brain age estimates for discovery sample ===
% ==================================================================
% /opt/matlab/bin/matlab -nodesktop -nodisplay -r "workingDir = pwd; run code/mri/brainageCollect.m"
fprintf('\n--- Collect & stack brain age estimates: discovery sample ---\n')

% set working directory and load dataset
fprintf(' - loading dataset.\n')
cd(workingDir)
load('results/mri/prepML.discovery.mat');

% get rvm estimates
fprintf(' - importing rvm predictions.\n')
pred_gm_rvm = zeros(1,100);
pred_wm_rvm = zeros(1,100);
for j = 1:100
	for i = 1:10
	load(sprintf('results/mri/ml.rvm/gm/rvm_%03d_%02d.mat',j,i), 'target_rvm')
    	pred_gm_rvm(u(:,j)==i,j) = target_rvm(u(:,j)==i,1);
    	clearvars target_rvm
    load(sprintf('results/mri/ml.rvm/wm/rvm_%03d_%02d.mat',j,i), 'target_rvm')
    	pred_wm_rvm(u(:,j)==i,j) = target_rvm(u(:,j)==i,1);
    	clearvars target_rvm
    end
end
tmp = [pred_wm_rvm pred_gm_rvm];
if sum(tmp(:)==0) == 0
    fprintf(' - no missing rvm predictions identified.\n')
end

% get xgboost estimates
fprintf(' - importing xgb predictions.\n')
pred_gm_xgb = importdata('results/mri/ml.xgb/gm/pred_brainage.txt');
pred_wm_xgb = importdata('results/mri/ml.xgb/wm/pred_brainage.txt');
tmp = [pred_gm_xgb pred_wm_xgb];
if sum(tmp(:)==0) == 0
    fprintf(' - no missing rvm predictions identified.\n')
end

% get brain age estimates
brainage_gm_rvm = mean(pred_gm_rvm,2);
brainage_gm_xgtree = mean(pred_gm_xgb(:,1:100),2);
brainage_gm_xglin = mean(pred_gm_xgb(:,101:200),2);

brainage_wm_rvm = mean(pred_wm_rvm,2);
brainage_wm_xgtree = mean(pred_wm_xgb(:,1:100),2);
brainage_wm_xglin = mean(pred_wm_xgb(:,101:200),2);

%% stack
% (pseudo)randomly generate 10 groups in each fold
fprintf(' - create nested groups for stacking.\n')
w = cell(10,100);
rng(8684 ,'twister')
for i = 1:100
    for j = 1:10
        for k = 1:10:sum(u(:,i)==j)
            w{j,i}(k:k+9,1) = randperm(10)';
        end
        w{j,i}(sum(u(:,i)==j)+1:size(w{j,i}),:)=[];
    end
end

% gm: create linear models for stacking
fprintf(' - stacking brain age estimates: grey matter.\n')
stack_models_gm = cell(10,10,100);
stack_predict_gm = NaN(size(u));
for i = 1:100
    for j = 1:10
        predictions = NaN(size(w{j,i}));
        for k = 1:10
            Train_Samples = [brainage_gm_xgtree(u(:,i)==j) brainage_gm_xglin(u(:,i)==j) brainage_gm_rvm(u(:,i)==j)]; Train_Samples = Train_Samples(w{j,i}~=k,:);
            Train_Targets = covs.table.t1_age; Train_Targets = Train_Targets(u(:,i)==j); Train_Targets = Train_Targets(w{j,i}~=k);
            Test_Samples = [brainage_gm_xgtree(u(:,i)==j) brainage_gm_xglin(u(:,i)==j) brainage_gm_rvm(u(:,i)==j)]; Test_Samples = Test_Samples(w{j,i}==k,:);
            stack_models_gm{k,j,i} = fitlm(Train_Samples, Train_Targets);
            predictions(w{j,i}==k) = predict(stack_models_gm{k,j,i},Test_Samples);           
        end
        stack_predict_gm(u(:,i)==j,i) = predictions;
    end
end

% wm: create linear models for stacking
fprintf(' - stacking brain age estimates: white matter.\n')
stack_models_wm = cell(10,10,100);
stack_predict_wm = NaN(size(u));
for i = 1:100
    for j = 1:10
        predictions = NaN(size(w{j,i}));
        for k = 1:10
            Train_Samples = [brainage_wm_xgtree(u(:,i)==j) brainage_wm_xglin(u(:,i)==j) brainage_wm_rvm(u(:,i)==j)]; Train_Samples = Train_Samples(w{j,i}~=k,:);
            Train_Targets = covs.table.t1_age; Train_Targets = Train_Targets(u(:,i)==j); Train_Targets = Train_Targets(w{j,i}~=k);
            Test_Samples = [brainage_wm_xgtree(u(:,i)==j) brainage_wm_xglin(u(:,i)==j) brainage_wm_rvm(u(:,i)==j)]; Test_Samples = Test_Samples(w{j,i}==k,:);
            stack_models_wm{k,j,i} = fitlm(Train_Samples, Train_Targets);
            predictions(w{j,i}==k) = predict(stack_models_wm{k,j,i},Test_Samples);           
        end
        stack_predict_wm(u(:,i)==j,i) = predictions;
    end
end

% gwm xgtree: create linear models for stacking
fprintf(' - stacking brain age estimates: grey and white matter (xgtree).\n')
stack_models_gwm_xgtree = cell(10,10,100);
stack_predict_gwm_xgtree = NaN(size(u));
for i = 1:100
    for j = 1:10
        predictions = NaN(size(w{j,i}));
        for k = 1:10
            Train_Samples = [brainage_gm_xgtree(u(:,i)==j) brainage_wm_xgtree(u(:,i)==j)]; Train_Samples = Train_Samples(w{j,i}~=k,:);
            Train_Targets = covs.table.t1_age; Train_Targets = Train_Targets(u(:,i)==j); Train_Targets = Train_Targets(w{j,i}~=k);
            Test_Samples = [brainage_gm_xgtree(u(:,i)==j) brainage_wm_xgtree(u(:,i)==j)]; Test_Samples = Test_Samples(w{j,i}==k,:);
            stack_models_gwm_xgtree{k,j,i} = fitlm(Train_Samples, Train_Targets);
            predictions(w{j,i}==k) = predict(stack_models_gwm_xgtree{k,j,i},Test_Samples);           
        end
        stack_predict_gwm_xgtree(u(:,i)==j,i) = predictions;
    end
end

% gwm xglin: create linear models for stacking
fprintf(' - stacking brain age estimates: grey and white matter (xglin).\n')
stack_models_gwm_xglin = cell(10,10,100);
stack_predict_gwm_xglin = NaN(size(u));
for i = 1:100
    for j = 1:10
        predictions = NaN(size(w{j,i}));
        for k = 1:10
            Train_Samples = [brainage_gm_xglin(u(:,i)==j) brainage_wm_xglin(u(:,i)==j)]; Train_Samples = Train_Samples(w{j,i}~=k,:);
            Train_Targets = covs.table.t1_age; Train_Targets = Train_Targets(u(:,i)==j); Train_Targets = Train_Targets(w{j,i}~=k);
            Test_Samples = [brainage_gm_xglin(u(:,i)==j) brainage_wm_xglin(u(:,i)==j)]; Test_Samples = Test_Samples(w{j,i}==k,:);
            stack_models_gwm_xglin{k,j,i} = fitlm(Train_Samples, Train_Targets);
            predictions(w{j,i}==k) = predict(stack_models_gwm_xglin{k,j,i},Test_Samples);           
        end
        stack_predict_gwm_xglin(u(:,i)==j,i) = predictions;
    end
end

% gwm rvm: create linear models for stacking
fprintf(' - stacking brain age estimates: grey and white matter (rvm).\n')
stack_models_gwm_rvm = cell(10,10,100);
stack_predict_gwm_rvm = NaN(size(u));
for i = 1:100
    for j = 1:10
        predictions = NaN(size(w{j,i}));
        for k = 1:10
            Train_Samples = [brainage_gm_rvm(u(:,i)==j) brainage_wm_rvm(u(:,i)==j)]; Train_Samples = Train_Samples(w{j,i}~=k,:);
            Train_Targets = covs.table.t1_age; Train_Targets = Train_Targets(u(:,i)==j); Train_Targets = Train_Targets(w{j,i}~=k);
            Test_Samples = [brainage_gm_rvm(u(:,i)==j) brainage_wm_rvm(u(:,i)==j)]; Test_Samples = Test_Samples(w{j,i}==k,:);
            stack_models_gwm_rvm{k,j,i} = fitlm(Train_Samples, Train_Targets);
            predictions(w{j,i}==k) = predict(stack_models_gwm_rvm{k,j,i},Test_Samples);           
        end
        stack_predict_gwm_rvm(u(:,i)==j,i) = predictions;
    end
end

% gwm: create linear models for stacking
fprintf(' - stacking brain age estimates: all tissues and algorithms.\n')
stack_models_gwm = cell(10,10,100);
stack_predict_gwm = NaN(size(u));
for i = 1:100
    for j = 1:10
        predictions = NaN(size(w{j,i}));
        for k = 1:10
            Train_Samples = [brainage_gm_xgtree(u(:,i)==j) brainage_gm_xglin(u(:,i)==j) brainage_gm_rvm(u(:,i)==j) brainage_wm_xgtree(u(:,i)==j) brainage_wm_xglin(u(:,i)==j) brainage_wm_rvm(u(:,i)==j)]; Train_Samples = Train_Samples(w{j,i}~=k,:);
            Train_Targets = covs.table.t1_age; Train_Targets = Train_Targets(u(:,i)==j); Train_Targets = Train_Targets(w{j,i}~=k);
            Test_Samples = [brainage_gm_xgtree(u(:,i)==j) brainage_gm_xglin(u(:,i)==j) brainage_gm_rvm(u(:,i)==j) brainage_wm_xgtree(u(:,i)==j) brainage_wm_xglin(u(:,i)==j) brainage_wm_rvm(u(:,i)==j)]; Test_Samples = Test_Samples(w{j,i}==k,:);
            stack_models_gwm{k,j,i} = fitlm(Train_Samples, Train_Targets);
            predictions(w{j,i}==k) = predict(stack_models_gwm{k,j,i},Test_Samples);           
        end
        stack_predict_gwm(u(:,i)==j,i) = predictions;
    end
end

% create structure
brainage = struct();
brainage.data = [brainage_gm_xgtree brainage_gm_xglin brainage_gm_rvm mean(stack_predict_gm,2) brainage_wm_xgtree brainage_wm_xglin brainage_wm_rvm mean(stack_predict_wm,2) mean(stack_predict_gwm_xgtree,2) mean(stack_predict_gwm_xglin,2) mean(stack_predict_gwm_rvm,2) mean(stack_predict_gwm,2)];
brainage.varnames = {'brainage_gm_xgtree','brainage_gm_xglin', 'brainage_gm_rvm', 'brainage_gm_stk', 'brainage_wm_xgtree', 'brainage_wm_xglin', 'brainage_wm_rvm', 'brainage_wm_stk', 'brainage_gwm_xgtree', 'brainage_gwm_xglin', 'brainage_gwm_rvm', 'brainage_gwm_stk', ...
    'brainage_gap_gm_xgtree','brainage_gap_gm_xglin', 'brainage_gap_gm_rvm', 'brainage_gap_gm_stk', 'brainage_gap_wm_xgtree', 'brainage_gap_wm_xglin', 'brainage_gap_wm_rvm', 'brainage_gap_wm_stk', 'brainage_gap_gwm_xgtree', 'brainage_gap_gwm_xglin', 'brainage_gap_gwm_rvm', 'brainage_gap_gwm_stk', ...
    'brainage_gap_adj_gm_xgtree','brainage_gap_adj_gm_xglin', 'brainage_gap_adj_gm_rvm', 'brainage_gap_adj_gm_stk', 'brainage_gap_adj_wm_xgtree', 'brainage_gap_adj_wm_xglin', 'brainage_gap_adj_wm_rvm', 'brainage_gap_adj_wm_stk', 'brainage_gap_adj_gwm_xgtree', 'brainage_gap_adj_gwm_xglin', 'brainage_gap_adj_gwm_rvm', 'brainage_gap_adj_gwm_stk'};

%% calculate brainage gap
age = covs.table.t1_age;
for i = 1:12
    brainage.data(:,i+12) = brainage.data(:,i) - age;
end
   
% adjust brainage gap by sex, age, age2, assessment center, TIV
age2 = covs.table.t1_age2;
sex = covs.table.sex;
ac1 = covs.table.t1_ac1;
ac2 = covs.table.t1_ac2;
TIV = meta.ratings(:,1);
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
fprintf('- saving brainage estimates.\n')
save('results/mri/brainage.discovery.mat', 'brainage', 'covs', 'gm_logical', 'meta', 'u', 'v', 'w', 'wm_logical')
colNames = [meta.ratings_varnames, brainage.varnames];
brainage_table = [covs.table array2table([meta.ratings, brainage.data],'VariableNames',colNames)];
writetable(brainage_table, 'results/mri/brainage.discovery.txt', 'Delimiter', '\t', 'WriteRowNames', 0)

% save stack models
fprintf('- storing stackmodels in structure.\n')
stackmodels = struct();
stackmodels.wm = stack_models_wm;
stackmodels.gm = stack_models_gm;
stackmodels.gwm = stack_models_gwm;
stackmodels.gwm_xgtree = stack_models_gwm_xgtree;
stackmodels.gwm_xglin = stack_models_gwm_xglin;
stackmodels.gwm_rvm = stack_models_gwm_rvm;

% average stack coefficients across all models (for singlemodel predictions)
fprintf('- calculating average stack coefficients (for singlemodel predictions).\n')
stack_gm = NaN(4,10,10,100);
stack_wm = NaN(4,10,10,100);
stack_gwm_xgtree = NaN(3,10,10,100);
stack_gwm_xglin = NaN(3,10,10,100);
stack_gwm_rvm = NaN(3,10,10,100);
stack_gwm = NaN(7,10,10,100);

    for i = 1:100
        for j = 1:10
            for k = 1:10
                coeff = table2array(stackmodels.gm{k,j,i}.Coefficients); stack_gm(:,k,j,i) = coeff(:,1);
                coeff = table2array(stackmodels.wm{k,j,i}.Coefficients); stack_wm(:,k,j,i) = coeff(:,1);
                coeff = table2array(stackmodels.gwm_xgtree{k,j,i}.Coefficients); stack_gwm_xgtree(:,k,j,i) = coeff(:,1);
                coeff = table2array(stackmodels.gwm_xglin{k,j,i}.Coefficients); stack_gwm_xglin(:,k,j,i) = coeff(:,1);
                coeff = table2array(stackmodels.gwm_rvm{k,j,i}.Coefficients); stack_gwm_rvm(:,k,j,i) = coeff(:,1);
                coeff = table2array(stackmodels.gwm{k,j,i}.Coefficients); stack_gwm(:,k,j,i) = coeff(:,1);
            end
        end
    end

    % reshape matrices and average coefficients
    stack_gm = reshape(stack_gm, [4,10000])';
    stack_wm = reshape(stack_wm, [4,10000])';
    stack_gwm_xgtree = reshape(stack_gwm_xgtree, [3,10000])';
    stack_gwm_xglin = reshape(stack_gwm_xglin, [3,10000])';
    stack_gwm_rvm = reshape(stack_gwm_rvm, [3,10000])';
    stack_gwm = reshape(stack_gwm, [7,10000])';

    stack_coefficients_gm = mean(stack_gm);
    stack_coefficients_wm = mean(stack_wm);
    stack_coefficients_gwm_xgtree = mean(stack_gwm_xgtree);
    stack_coefficients_gwm_xglin = mean(stack_gwm_xglin);
    stack_coefficients_gwm_rvm = mean(stack_gwm_rvm);
    stack_coefficients_gwm = mean(stack_gwm);

    % store coefficients in structure
    stack_coefficients = struct();
    stack_coefficients.gm = stack_coefficients_gm;
    stack_coefficients.wm = stack_coefficients_wm;
    stack_coefficients.gwm_xgtree = stack_coefficients_gwm_xgtree;
    stack_coefficients.gwm_xglin = stack_coefficients_gwm_xglin;
    stack_coefficients.gwm_rvm = stack_coefficients_gwm_rvm;
    stack_coefficients.gwm = stack_coefficients_gwm;

    stack_coefficients.gm_varnames = {'intercept', 'gm_xgtree', 'gm_xglin', 'gm_rvm'};
    stack_coefficients.wm_varnames = {'intercept', 'wm_xgtree', 'wm_xglin', 'wm_rvm'};
    stack_coefficients.gwm_xgtree_varnames = {'intercept', 'gm_xgtree', 'wm_xgtree'};
    stack_coefficients.gwm_xglin_varnames = {'intercept', 'gm_xglin', 'wm_xglin'};
    stack_coefficients.gwm_rvm_varnames = {'intercept', 'gm_rvm', 'wm_rvm'};
    stack_coefficients.gwm_varnames = {'intercept', 'gm_xgtree', 'gm_xglin', 'gm_rvm', 'wm_xgtree', 'wm_xglin', 'wm_rvm'};

    % save stackmodels
    fprintf('- saving stackmodels.\n')
    save('results/mri/brainage.stackmodels.mat', 'stackmodels', 'stack_coefficients', '-v7.3')

% quit matlab
fprintf('\n--- Completed: Collect & stack brain age estimates: discovery sample ---\n')
exit
