% ===================================================
% === collect brain age estimates for life sample ===
% ===================================================
% brain age models have been applied in the same manner as for the UKB replication sample (singlemodel)
% /opt/matlab/bin/matlab -nodesktop -nodisplay -r "workingDir = pwd; run code/mri/brainageCollect.m"
fprintf('\n--- Collect brain age estimates: LIFE-Adult sample ---\n')

% set working directory
cd(workingDir)

% load predictions made in LIFE dataset
fprintf(' - loading dataset.\n')
T = readtable('data/life/pv532/pv532_brainage.txt', 'TreatAsEmpty', 'NA');
qcX = readtable('data/life/gwas/s304_1_Inds_included.txt', 'TreatAsEmpty', 'NA');
discovery = load('results/mri/brainage.discovery.mat');

% put variables into structures
fprintf(' - putting variables into structures.\n')
brainage = struct();
brainage.data = T{:,20:31};
brainage.varnames = {'brainage_gm_xgtree','brainage_gm_xglin', 'brainage_gm_rvm', 'brainage_gm_stack', 'brainage_wm_xgtree', 'brainage_wm_xglin', 'brainage_wm_rvm', 'brainage_wm_stack', 'brainage_gwm_xgtree', 'brainage_gwm_xglin', 'brainage_gwm_rvm', 'brainage_gwm_stack', ...
    'brainage_gap_gm_xgtree','brainage_gap_gm_xglin', 'brainage_gap_gm_rvm', 'brainage_gap_gm_stack', 'brainage_gap_wm_xgtree', 'brainage_gap_wm_xglin', 'brainage_gap_wm_rvm', 'brainage_gap_wm_stack', 'brainage_gap_gwm_xgtree', 'brainage_gap_gwm_xglin', 'brainage_gap_gwm_rvm', 'brainage_gap_gwm_stack', ...
    'brainage_gap_adj_gm_xgtree','brainage_gap_adj_gm_xglin', 'brainage_gap_adj_gm_rvm', 'brainage_gap_adj_gm_stack', 'brainage_gap_adj_wm_xgtree', 'brainage_gap_adj_wm_xglin', 'brainage_gap_adj_wm_rvm', 'brainage_gap_adj_wm_stack', 'brainage_gap_adj_gwm_xgtree', 'brainage_gap_adj_gwm_xglin', 'brainage_gap_adj_gwm_rvm', 'brainage_gap_adj_gwm_stack'};

meta = struct();
meta.ratings = T{:,10:19};
meta.ratings_varnames = T.Properties.VariableNames(10:19);
meta.IID = T{:,1:3};

covs = struct();
covs.table = T(:,[2:3 4:9]);
covs.table.Properties.VariableNames = {'PSEUDONYM','ADULT_SNP_SAMPLING_ID','ADULT_SNP_ARRAY_QUALITY','MRI_USABLE','sex','age','age2','date'};
covs.colheaders = covs.table.Properties.VariableNames;
covs.IID = T{:,2};

% exclude subjects without genetics, with mri labelled as 'not usable', out of age range (45.1673 - 81.8635), or with poor IQR
fprintf(' - removing individualds who fail qc / who are out of UKB age range.\n')
select_logical = (isnan(covs.table.ADULT_SNP_ARRAY_QUALITY) + (covs.table.MRI_USABLE == 3) + (covs.table.age < min(discovery.covs.table.t1_age)) + (covs.table.age > max(discovery.covs.table.t1_age)) + (meta.ratings(:,9) >= 3)) < 1;
brainage.data = brainage.data(select_logical,:);
meta.ratings = meta.ratings(select_logical,:);
meta.IID = meta.IID(select_logical,:);
covs.table = covs.table(select_logical,:);
covs.IID = covs.IID(select_logical,:);

% remove individuals who failed chromosome X quality checks
arrayID = cell2table(covs.table.ADULT_SNP_SAMPLING_ID, 'VariableNames',{'ADULT_SNP_ARRAY_ID'});
select_logical = ismember(arrayID,qcX);
brainage.data = brainage.data(select_logical,:);
meta.ratings = meta.ratings(select_logical,:);
meta.IID = meta.IID(select_logical,:);
covs.table = covs.table(select_logical,:);
covs.IID = covs.IID(select_logical,:);

%% calculate brainage gap
age = covs.table.age;
for i = 1:12
    brainage.data(:,i+12) = brainage.data(:,i) - age;
end
   
% adjust brainage gap by sex, age, age2, TIV
sex = covs.table.sex;
age2 = covs.table.age2;
TIV = meta.ratings(:,1);
for i = 1:12
    model = fitlm([sex age age2 TIV], brainage.data(:,i+12));
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
save('results/mri/brainage.life.mat', 'brainage', 'covs', 'meta')
colNames = [meta.ratings_varnames, brainage.varnames];
brainage_table = [covs.table array2table([meta.ratings, brainage.data],'VariableNames',colNames)];
writetable(brainage_table, 'results/mri/brainage.life.txt', 'Delimiter', '\t', 'WriteRowNames', 0)

% quit matlab
fprintf('\n--- Completed: Collect brain age estimates: LIFE-Adult sample ---\n')
exit
