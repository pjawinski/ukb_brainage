% =============================================================
% === PCA in Matlab (singlemodel withouth cross-validation) ===
% =============================================================

% script requires prior definition of...
% matFile, e.g. matFile = '/slow/projects/ukb_brainage/results/mri/prepML.discovery.mat')
% tissue, e.g. tissue = 'gm')
% i (fold in n-fold cross-validation), e.g. i = 7
% j (repeat in n-fold cross-validation), e.g. j = 23
% maxNumCompThreads(5)
tic

% load data
dataset = load(matFile, 'u', 'v', tissue);

% get variables
data = dataset.(tissue);
u = dataset.u;
v = dataset.v;

% Define samples
Train_Samples = data(u(:,j)~= v(i,j),:);
Validation_Samples = data(u(:,j)==v(i,j),:);

% run pca
[Train_Samples_pca_coeff,Train_Samples_pca_score,latent] = pca(Train_Samples);

% get pca scores of the first 500 principal components
Train_Samples_pca_score = Train_Samples_pca_score(:,1:500);
Train_Samples_pca_coeff = Train_Samples_pca_coeff(:,1:500);
Train_Samples_means = mean(Train_Samples);

% calculate pca scores for validation and test sample 
Validation_Samples_centered = Validation_Samples-repmat(Train_Samples_means,size(Validation_Samples,1),1);
Validation_Samples_pca_score = Validation_Samples_centered/Train_Samples_pca_coeff';

% save data
dlmwrite(sprintf('xgb_%s_train_samples_pca.txt', tissue), Train_Samples_pca_score, 'delimiter', '\t', 'precision', 16);
dlmwrite(sprintf('xgb_%s_validation_samples_pca.txt', tissue), Validation_Samples_pca_score, 'delimiter', '\t', 'precision', 16);
save(sprintf('xgb_%s_pca_training.mat', tissue), 'Train_Samples_means', 'Train_Samples_pca_coeff')

% save eigenvalues
eigen = latent/sum(latent);
eigen = eigen(1:500);
eigenFormatted = arrayfun(@(x) sprintf('%.12f', x), eigen, 'UniformOutput', false);
eigenTable = table(eigenFormatted, 'VariableNames', {'Eigenvalues'});
writetable(eigenTable, sprintf('xgb_%s_eigenvalues.txt', tissue), ...
    'Delimiter', '\t', 'WriteVariableNames', true, 'FileType', 'text');

% exit
fprintf(' - matlab pca took %d seconds.\n', round(toc,0))
exit