% ====================================
% === run relevance vector machine ===
% ====================================
% /opt/matlab/bin/matlab -nodesktop -nodisplay -r "tissue = 'gm'; workingDir = pwd; targetDir = 'results/mri/ml.rvm/gm/'; matFile = 'results/mri/prepML.discovery.mat';  maxNumCompThreads(100); spmPath = '/fast/software/matlab/spm12/'; rvmPath = '/fast/software/matlab/RVM/'; sbPath = '/fast/software/matlab/RVM/SB2_Release_200/'; run code/mri/ml_rvm.m"

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
load(matFile)
x = eval(tissue);
y = covs.table.t1_age;
nfolds = max(u(:));
nrepeats = size(u,2);

% loop over repeats and folds
system(sprintf('mkdir -p %s', targetDir));
for j = 1:nrepeats

    % loop over folds
    for i = 1:nfolds

        % divide into training and test data
        tic
        fprintf('\nStarting with rvm_%03d_%02d.\n', j, i)
        Train_Samples = x(u(:,j)~=i,:);
        Train_Targets = y(u(:,j)~=i,:);
        Test_Samples = x(u(:,j)==i,:);

        % create empty outcome vector
        target_rvm = zeros(size(y,1),2);

        % run pca and keep first 500 components
        fprintf(' - running pca.\n')
        [Train_Samples_pca_coeff,Train_Samples_pca_score] = pca(Train_Samples);
        Train_Samples_pca_score = Train_Samples_pca_score(:,1:500);
        Train_Samples_pca_coeff = Train_Samples_pca_coeff(:,1:500);
        Train_Samples_means = mean(Train_Samples);
        
        % run rvm
        fprintf(' - training rvm model.\n')
        [RVM] = rvm_train_5000(Train_Samples_pca_score,Train_Targets);
        
        % make rvm structure sparser
        RVM.X = RVM.X(RVM.rv_index(2:end)-1,:);
        RVM = rmfield(RVM, 'y_mu');
        RVM = rmfield(RVM, 'y_var');
        RVM = rmfield(RVM, 'rv_index');
        RVM.rv_index = (1:size(RVM.rv_mu,1))';

        % make predictions
        Test_Samples_centered = Test_Samples-repmat(mean(Train_Samples),size(Test_Samples,1),1);
        Test_Samples_pca_score = Test_Samples_centered/Train_Samples_pca_coeff';
        [mu, var] = rvm_test(RVM,Test_Samples_pca_score);

        % store results in single variable
        target_rvm(u(:,j)==i,:) = [mu var];

        % save results
        fprintf(' - saving model predictions.\n')
        MatFileName=sprintf('%s/rvm_%03d_%02d.mat', targetDir, j, i);
        LogFileName=sprintf('%s/rvm_%03d_%02d.log', targetDir, j, i);
        parsave(MatFileName, RVM, Train_Samples_means, Train_Samples_pca_coeff, target_rvm);
        system(sprintf('mv logfile.txt %s', LogFileName));

        % display prediction performance
        cor = corr(y(u(:,j)==i,:), mu);
        MAE = mean(abs(y(u(:,j)==i,:)-mu));
        fprintf(' - rvm_%03d_%02d - Pearson: %.4f - MAE: %.4f.\n', j, i, cor, MAE)
        toc
    end
end
fprintf('Completed: RVM')

% quit matlab
exit

