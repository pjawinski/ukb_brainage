% =====================
% === PCA in Matlab ===
% =====================

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
Train_Samples = data(u(:,j)~=i & u(:,j)~= v(i,j),:);
Validation_Samples = data(u(:,j)==v(i,j),:);
Test_Samples = data(u(:,j)==i,:);

% run pca
[Train_Samples_pca_coeff,Train_Samples_pca_score] = pca(Train_Samples);

% get pca scores of the first 500 principal components
Train_Samples_pca_score = Train_Samples_pca_score(:,1:500);
Train_Samples_pca_coeff = Train_Samples_pca_coeff(:,1:500);
Train_Samples_means = mean(Train_Samples);

% calculate pca scores for validation and test sample 
Validation_Samples_centered = Validation_Samples-repmat(Train_Samples_means,size(Validation_Samples,1),1);
Validation_Samples_pca_score = Validation_Samples_centered/Train_Samples_pca_coeff';
Test_Samples_centered = Test_Samples-repmat(Train_Samples_means,size(Test_Samples,1),1);
Test_Samples_pca_score = Test_Samples_centered/Train_Samples_pca_coeff';
    
% save data
dlmwrite('Train_Samples_pca.txt', Train_Samples_pca_score, 'delimiter', '\t', 'precision', 16);
dlmwrite('Validation_Samples_pca.txt', Validation_Samples_pca_score, 'delimiter', '\t', 'precision', 16);
dlmwrite('Test_Samples_pca.txt', Test_Samples_pca_score, 'delimiter', '\t', 'precision', 16);
save pca_training.mat Train_Samples_means Train_Samples_pca_coeff
	
% exit
fprintf(' - matlab pca for fold %d took %d seconds.\n', i, round(toc,0))
exit