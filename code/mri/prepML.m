% =====================================================
% === prepare for machine learning / age prediction ===
% =====================================================
% /opt/matlab/bin/matlab -nodesktop -nodisplay -r "run code/mri/prepML.m"
fprintf('\n--- Preparing datasets for machine learning / age prediction. ---\n')

% set working directory
cd(workingDir)

% ----- discovery sample ----
% load data
fprintf('Processing discovery dataset.\n')
fprintf(' - loading data.\n')
load('results/mri/cat12.r2024.mat')
covs = struct();
covs.table = readtable('results/mri/r2024.vars.txt','TreatAsMissing', 'NA'); % covs = importdata('/volumes/psymol/projects/UK_Biobank/04_data_genetics_linux/00_script/release_feb2020/02_covs.txt', '\t', 1);
covs.colheaders = covs.table.Properties.VariableNames;

% select discovery cohort (r2020 white-British ancestry)
% align covs and cat12 data
fprintf(' - aligning covs and cat12 data.\n')
covs.table = covs.table(covs.table.discovery==1,:);
[~,IA,IB] = intersect(meta.IID, covs.table.IID, 'stable') ;
covs.table = covs.table(IB,:);
meta.IID = meta.IID(IA,:);
meta.files = meta.files(IA,:);
meta.ratings = meta.ratings(IA,:);
gm = gm(IA,:);
wm = wm(IA,:);
if isequal(covs.table.IID,meta.IID) == 1
    fprintf(' - succesfully aligned.\n')
end

% remove no-variance voxels
fprintf(' - removing no-variance voxels.\n')
gm_logical = std(gm)~=0;
gm = gm(:,gm_logical);
wm_logical = std(wm)~=0;
wm = wm(:,wm_logical);

% sort data by age 
fprintf(' - sorting data by age.\n')
[~, I] = sort(covs.table.t1_age);
covs.table = covs.table(I,:);
meta.IID = meta.IID(I,:);
meta.files = meta.files(I,:);
meta.ratings = meta.ratings(I,:);
gm = gm(I,:);
wm = wm(I,:);

% (pseudo)randomly generate 10 groups (wih 100 repeats)
fprintf(' - for machine learning: pseudo-randomly generate 10 groups (with 100 repeats).\n')
u = zeros(size(meta.IID,1),100);
rng(4616,'twister')

for ii = 1:100
    for i = 1:10:size(meta.IID,1)
        u(i:i+9,ii) = randperm(10)';
    end
end
u(size(meta.IID,1)+1:size(u,1),:)=[];

% (pseudo)randomly generate validation group
fprintf(' - for machine learning: pseudo-randomly generate validation groups (with 100 repeats).\n')
v = zeros(10,100);
for j = 1:100
    v(:,j) = 1:10;
    while sum((v(:,j) - [1:10]')==0)>0
        v(:,j) = randperm(10)';
    end
end

% save discovery dataset
fprintf(' - saving discovery dataset in .mat and .txt files (for R import).\n\n')
save('results/mri/prepML.discovery.mat', 'covs', 'meta', 'gm', 'gm_logical', 'wm', 'wm_logical', 'u', 'v')
dlmwrite('results/mri/prepML.discovery.gm.txt', gm, 'delimiter','\t', 'precision', 16)
dlmwrite('results/mri/prepML.discovery.wm.txt', wm, 'delimiter','\t', 'precision', 16)
dlmwrite('results/mri/prepML.discovery.u.txt', u, 'delimiter','\t')
dlmwrite('results/mri/prepML.discovery.v.txt', v, 'delimiter','\t')
writetable(covs.table,'results/mri/prepML.discovery.covs.txt', 'delimiter','\t', 'WriteRowNames', 0)
T = [meta.IID, meta.ratings]; colNames = ['IID', meta.ratings_varnames]; T = array2table(T,'VariableNames',colNames);
writetable(T, 'results/mri/prepML.discovery.meta.txt', 'Delimiter', '\t', 'WriteRowNames', 0)

% ----- replication sample ----
% load
fprintf('Processing replication dataset.\n')
fprintf(' - loading data.\n')
load('results/mri/cat12.r2024.mat')
covs = struct();
covs.table = readtable('results/mri/r2024.vars.txt','TreatAsMissing', 'NA'); % covs = importdata('/volumes/psymol/projects/UK_Biobank/04_data_genetics_linux/00_script/release_feb2020/02_covs.txt', '\t', 1);
covs.colheaders = covs.table.Properties.VariableNames;

% select replication cohort (r2021 multi-ancestry)
% align covs and cat12 data
fprintf(' - aligning covs and cat12 data.\n')
covs.table = covs.table(covs.table.discovery==0,:);
[~,IA,IB] = intersect(meta.IID, covs.table.IID, 'stable') ;
covs.table = covs.table(IB,:);
meta.IID = meta.IID(IA,:);
meta.files = meta.files(IA,:);
meta.ratings = meta.ratings(IA,:);
gm = gm(IA,:);
wm = wm(IA,:);
if isequal(covs.table.IID,meta.IID) == 1
    fprintf(' - succesfully aligned.\n')
end

% remove no-variance voxels (logical from discovery cohort)
fprintf(' - removing no-variance voxels (based on discovery dataset).\n')
gm = gm(:,gm_logical);
wm = wm(:,wm_logical);

% sort data by age 
fprintf(' - sorting data by age.\n')
[~, I] = sort(covs.table.t1_age);
covs.table = covs.table(I,:);
meta.IID = meta.IID(I,:);
meta.files = meta.files(I,:);
meta.ratings = meta.ratings(I,:);
gm = gm(I,:);
wm = wm(I,:);

% save replication dataset
fprintf(' - saving replication dataset in .mat and .txt files (for R import).\n\n')
save('results/mri/prepML.replication.mat', 'covs', 'meta', 'gm', 'gm_logical', 'wm', 'wm_logical')
dlmwrite('results/mri/prepML.replication.gm.txt', gm, 'delimiter','\t', 'precision', 16)
dlmwrite('results/mri/prepML.replication.wm.txt', wm, 'delimiter','\t', 'precision', 16)
writetable(covs.table,'results/mri/prepML.replication.covs.txt', 'delimiter','\t', 'WriteRowNames', 0)
T = [meta.IID, meta.ratings]; colNames = ['IID', meta.ratings_varnames]; T = array2table(T,'VariableNames',colNames);
writetable(T, 'results/mri/prepML.replication.meta.txt', 'Delimiter', '\t', 'WriteRowNames', 0)

% ----- retest data ----
% load
fprintf('Processing retest dataset.\n')
fprintf(' - loading data.\n')
load('results/mri/cat12.r2024.retest.mat')
covs = struct();
covs.table = readtable('results/mri/r2024.vars.txt','TreatAsMissing', 'NA'); % covs = importdata('/volumes/psymol/projects/UK_Biobank/04_data_genetics_linux/00_script/release_feb2020/02_covs.txt', '\t', 1);
covs.colheaders = covs.table.Properties.VariableNames;

% only keep valid retest

% align covs and cat12 data
fprintf(' - aligning covs and cat12 data.\n')
[~,IA,IB] = intersect(meta.IID, covs.table.IID, 'stable') ;
covs.table = covs.table(IB,:);
meta.IID = meta.IID(IA,:);
meta.files = meta.files(IA,:);
meta.ratings = meta.ratings(IA,:);
gm = gm(IA,:);
wm = wm(IA,:);
if isequal(covs.table.IID,meta.IID) == 1
    fprintf(' - succesfully aligned.\n')
end

% remove no-variance voxels (logical from discovery cohort)
fprintf(' - removing no-variance voxels (based on discovery dataset).\n')
gm = gm(:,gm_logical);
wm = wm(:,wm_logical);

% sort data by age 
fprintf(' - sorting data by age.\n')
[~, I] = sort(covs.table.t1_age);
covs.table = covs.table(I,:);
meta.IID = meta.IID(I,:);
meta.files = meta.files(I,:);
meta.ratings = meta.ratings(I,:);
gm = gm(I,:);
wm = wm(I,:);

% save retest dataset
fprintf(' - saving retest dataset in .mat and .txt files (for R import).\n\n')
save('results/mri/prepML.retest.mat', 'covs', 'meta', 'gm', 'gm_logical', 'wm', 'wm_logical')
dlmwrite('results/mri/prepML.retest.gm.txt', gm, 'delimiter','\t', 'precision', 16)
dlmwrite('results/mri/prepML.retest.wm.txt', wm, 'delimiter','\t', 'precision', 16)
writetable(covs.table,'results/mri/prepML.retest.covs.txt', 'delimiter','\t', 'WriteRowNames', 0)
T = [meta.IID, meta.ratings]; colNames = ['IID', meta.ratings_varnames]; T = array2table(T,'VariableNames',colNames);
writetable(T, 'results/mri/prepML.retest.meta.txt', 'Delimiter', '\t', 'WriteRowNames', 0)
fprintf('Completed.\n\n')
exit

