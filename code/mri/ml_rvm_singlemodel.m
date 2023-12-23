% ==================================================
% === run relevance vector machine (singlemodel) ===
% ==================================================
% /opt/matlab/bin/matlab -nodesktop -nodisplay -r "tissue = 'gm'; workingDir = pwd; targetDir = 'results/mri/ml.rvm/singlemodel/'; matFile = 'results/mri/prepML.discovery.mat';  maxNumCompThreads(100); spmPath = '/fast/software/matlab/spm12/'; rvmPath = '/fast/software/matlab/RVM/'; sbPath = '/fast/software/matlab/RVM/SB2_Release_200/'; run code/mri/ml_rvm.m"

% show input variables
fprintf('\n--- RVM settings ---\ntissue: %s\nworkingDir: %s\ntargetDir: %s\nmatFile: %s\nspmPath: %s\nrvmPath: %s\nsbPath: %s\nthreads: %d\n\n',...
    tissue, workingDir, targetDir, matFile, spmPath, rvmPath, sbPath, maxNumCompThreads())

% set working directory
cd(workingDir)

% addpath
addpath('code/functions/')
addpath(spmPath)
addpath(rvmPath)
addpath(sbPath)

% load dataset and get variables
fprintf(' - loading data.\n')
load(matFile)
x = eval(tissue);
y = covs.table.t1_age;
nfolds = max(u(:));
nrepeats = size(u,2);

% create target folder
system(sprintf('mkdir -p %s', targetDir));

% divide into training and test data
Train_Samples = x;
Train_Targets = y;

% run pca and keep first 500 components
fprintf(' - running pca.\n')
[Train_Samples_pca_coeff,Train_Samples_pca_score] = pca(Train_Samples);
Train_Samples_pca_score = Train_Samples_pca_score(:,1:500);
Train_Samples_pca_coeff = Train_Samples_pca_coeff(:,1:500);
Train_Samples_means = mean(Train_Samples);
        
% run rvm
fprintf(' - training rvm model.\n')
[RVM] = rvm_train(Train_Samples_pca_score,Train_Targets);

% make rvm structure sparser
RVM.X = RVM.X(RVM.rv_index(2:end)-1,:);
RVM = rmfield(RVM, 'y_mu');
RVM = rmfield(RVM, 'y_var');
RVM = rmfield(RVM, 'rv_index');
RVM.train_means = mean(Train_Samples);
RVM.train_pca_coeff = Train_Samples_pca_coeff;
RVM.rv_index = (1:size(RVM.rv_mu,1))';

% save model and logfile
fprintf(' - saving model.\n')
save(sprintf('%s/rvm_%s.mat', targetDir, tissue), 'RVM')
LogFileName=sprintf('%s/rvm_%s.log', targetDir, tissue);
system(sprintf('mv logfile.txt %s', LogFileName));

% quit matlab
fprintf('Completed: RVM (singlemodel)\n')
exit

