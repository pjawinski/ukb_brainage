% ==========================
% === collect cat12 data ===
% ==========================
% /opt/matlab/bin/matlab -nodesktop -nodisplay -r "scanDir = {'data/t1w/r2020/','data/t1w/r2021/'}; outfile = 'results/mri/cat12.r2020'; spmPath = '/fast/software/matlab/spm12/'; ncores = 50; run code/mri/cat12.collect.m"

% set working directory
cd /slow/projects/ukb_brainage

% add folder with 'read_nifti' function
addpath('code/functions/')
addpath(spmPath)

% find _cat12.zip files
fprintf('[1/5] Finding _cat12.zip files.\n')
for i = 1:numel(scanDir)
    if i == 1
        files = dir(strcat(scanDir{i},'**/*_cat12.zip'));
    else
        tmp = dir(strcat(scanDir{i},'**/*_cat12.zip'));
        files = [files;tmp];
    end
end
fprintf(' - %d files found.\n', numel(files))

% get individual identifier from files.name (numeric in UK Biobank)
IID = zeros(length(files),1);
for i = 1:length(files)
    IID(i,1) = str2num(files(i).name(1:(end-10)));
end

% start parallel pool
fprintf('[2/5] Getting cat12 grey and white matter segmentations.\n')
myCluster = parcluster('local');
myCluster.NumWorkers = ncores;
saveProfile(myCluster);
parpool('local',ncores);

% unzip rsmwp[1/2]T1_orig_defaced.nii, read nifti, and create grey and white matter matrix
files_constant = parallel.pool.Constant(files);
gm = zeros(numel(files), 16128);
wm = zeros(numel(files), 16128);

parfor i = 1:numel(files)
    files_c = files_constant.Value; 
    system(char(strcat({'unzip -joqq '}, files_c(i).folder, '/', files_c(i).name, {' cat12/mri/rsmwp1T1_orig_defaced.nii -d '}, files_c(i).folder, '/')));
    nifti = read_nifti(char(strcat(files_c(i).folder, '/rsmwp1T1_orig_defaced.nii')));
    gm(i,:) = nifti(:)';
    system(char(strcat({'rm -f '}, files_c(i).folder, '/rsmwp1T1_orig_defaced.nii'))); 
    nifti = []
    
    system(char(strcat({'unzip -joqq '}, files_c(i).folder, '/', files_c(i).name, {' cat12/mri/rsmwp2T1_orig_defaced.nii -d '}, files_c(i).folder, '/')));
    nifti = read_nifti(char(strcat(files_c(i).folder, '/rsmwp2T1_orig_defaced.nii')));
    wm(i,:) = nifti(:)';
    system(char(strcat({'rm -f '}, files_c(i).folder, '/rsmwp2T1_orig_defaced.nii'))); 
    nifti = []  
end

% have all files been successfully imported (no empty lines)?
gm_nonEmpty = sum(sum(gm, 2)~=0); % should be 42796 / 3119
wm_nonEmpty = sum(sum(wm, 2)~=0); % should be 42796 / 3119
fprintf(' - %d non-empty grey matter files.\n', gm_nonEmpty)
fprintf(' - %d non-empty white matter files.\n', wm_nonEmpty)

% unzip cat12 files and collect quality ratings
% subject measures
fprintf('[3/5] Getting cat12 meta file info.\n')
meta_ratings = zeros(numel(files),9);
parfor i = 1:numel(files)
    files_c = files_constant.Value;  
    system(char(strcat({'unzip -joqq '}, files(i).folder, '/', files(i).name, {' cat12/report/cat_T1_orig_defaced.mat -d '}, files(i).folder, '/')));
    tmp = parload(strcat(files(i).folder, '/cat_T1_orig_defaced.mat')); 
    tmp = [tmp.S.subjectmeasures.vol_TIV, tmp.S.subjectmeasures.vol_abs_CGW(:,1:4) tmp.S.qualityratings.res_RMS tmp.S.qualityratings.NCR tmp.S.qualityratings.ICR tmp.S.qualityratings.IQR];
    meta_ratings(i,:) = tmp;
    system(char(strcat({'rm -f '}, files(i).folder, '/cat_T1_orig_defaced.mat')));
end

% have all files been successfully opened?
meta_nonEmpty = sum(sum(meta_ratings, 2)~=0); % should be 42796 / 3119
fprintf(' - %d non-empty meta files.\n', meta_nonEmpty)

% merge variable
meta = struct();
meta.files = files;
meta.ratings = meta_ratings;
meta.ratings_varnames = {'TIV', 'CSF', 'GM', 'WM', 'RES', 'RMS','NCR','ICR', 'IQR'};
meta.IID = IID;

% find brain scans with bad quality
meta.ratings(:,10) = meta.ratings(:,9) > 3; 
meta.ratings_varnames = {'TIV', 'CSF', 'GM', 'WM', 'RES', 'RMS','NCR','ICR', 'IQR', 'IQR_poor'};

% save file
fprintf('[4/5] Saving data in file %s.mat\n', outfile)
save(strcat(outfile,'.mat'), 'meta', 'gm', 'wm', '-v7.3')

% write table
% merge variables
fprintf('[5/5] Saving meta data in file %s.txt\n', outfile)
merged = [num2cell(meta.IID), num2cell(meta.ratings)];
colNames = {'IID', 'TIV', 'CSF', 'GM', 'WM', 'VRES', 'RMS','NCR','ICR', 'IQR', 'IQR_poor'};
output = array2table(merged,'VariableNames',colNames);
writetable(output, strcat(outfile,'.txt'), 'Delimiter', '\t', 'WriteRowNames', 0)
exit
