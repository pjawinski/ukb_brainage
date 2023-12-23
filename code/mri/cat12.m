% ===========================
% === CAT12 Preprocessing ===
% ===========================
% /opt/matlab/bin/matlab -nodesktop -nodisplay -r "scanDir = 'data/t1w/r2020/'; filetype = '_20252_2_0.zip'; spmPath = '/fast/software/matlab/spm12/'; ncores = 50; run code/mri/cat12.m"

% set working directory
cd /slow/projects/ukb_brainage

% set Matlab toolbox paths
addpath('code/functions/')
addpath(spmPath)

% set temporary folder
temppath = strcat(scanDir,'tmp/');
system(char(strcat({'mkdir -p '}, temppath)));

% get files that require cat12 preprocessing
filetype_length = length(filetype);
scans = dir(strcat(scanDir, '/*/*', filetype));

clearvars scans2remove
k = 0;
for i = 1:size(scans,1)
    if exist(strcat(scans(i).folder,'/',scans(i).name(1:(end-filetype_length)),'_cat12.zip')) == 2
        k = k + 1;
        scans2remove(k,1) = i;
    end
end
if exist('scans2remove', 'var') == 1
    scans(scans2remove) = []
end

% get subs and folder
subs = extractBefore({scans.name}', filetype);
subs_folder = {scans.folder}';

% remove subs that have been identified as 'unusable'
if exist(strcat(scanDir,'/unusable.txt')) == 2
    unusable = readtable(strcat(scanDir,'/unusable.txt'));
    unusable = table2cell(unusable);
    unusable = cellfun(@num2str,unusable,'un',0);
    idx = ~ismember(subs,unusable);
    subs = subs(idx);
    scans = scans(idx);
    subs_folder = subs_folder(idx);
end

% set up Matlab job (cat12.5 build 1364)
matlabbatch{1}.spm.tools.cat.estwrite.nproc = 0;
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP = 1070;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.darteltpm = {strcat(spmPath,'/toolbox/cat12/templates_1.50mm/Template_1_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shootingtpm = {strcat(spmPath,'/toolbox/cat12/templates_1.50mm/Template_0_IXI555_MNI152_GS.nii')}; 
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regstr = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = 1.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.restypes.fixed = [1 0.1];
matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 2;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 2;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobian.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [0 0];
matlabbatch{2}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{2}.spm.spatial.smooth.dtype = 0;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's';

% start parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = ncores;
saveProfile(myCluster);
parpool('local',ncores);

% submit constant variables to workers
matlabbatch_constant = parallel.pool.Constant(matlabbatch);
subs_constant = parallel.pool.Constant(subs);
subs_folder_constant = parallel.pool.Constant(subs_folder);
temppath_constant = parallel.pool.Constant(temppath);

% initiate spm
spm('Defaults','fMRI');
spm_jobman('initcfg');

% loop
system(char(strcat({'rm -f '}, scanDir, '/errors.txt')));
parfor idx = 1:numel(subs)
    
    % get vars
    temp = temppath_constant.Value; % temp = temppath;
    sub = string(subs_constant.Value(idx)); % sub = string(subs(idx));
    sub_folder = string(subs_folder_constant.Value(idx)); % sub_folder = subs_folder(idx);
    matlabbatch_loop = matlabbatch_constant.Value; % matlabbatch_loop = matlabbatch;

    try
        % copy file to temporary folder and gunzip .nii.gz for spm
        system(strcat({'mkdir -p '}, temp, sub, {'; cp -R '}, sub_folder, {'/. '}, temp, sub));   
        system(strcat({'unzip '}, temp, sub, '/', sub, {'_20252_*_0.zip T1/T1_orig_defaced.nii.gz -d '}, temp, sub,'/')); %system(strcat({'gunzip -fk '}, temp, sub, '/', sub, '.nii.gz')); 

        if exist(strcat(temp, sub, '/T1/T1_orig_defaced.nii.gz'), 'file') == 2
            system(strcat({'gunzip '}, temp, sub, '/T1/T1_orig_defaced.nii.gz'));
            system(strcat({'mkdir -p '}, temp, sub, '/cat12'));
            system(strcat({'mv '}, temp, sub, {'/T1/T1_orig_defaced.nii '}, temp, sub, '/cat12'));

            matlabbatch_loop{1}.spm.tools.cat.estwrite.data = {char(strcat(temp, sub, '/cat12/T1_orig_defaced.nii'))};
            matlabbatch_loop{2}.spm.spatial.smooth.data =  {char(strcat(temp, sub, '/cat12/mri/mwp1T1_orig_defaced.nii'))
                                                        char(strcat(temp, sub, '/cat12/mri/mwp2T1_orig_defaced.nii'))
                                                        char(strcat(temp, sub, '/cat12/mri/wmT1_orig_defaced.nii'))};
            spm_jobman('run',matlabbatch_loop);

            resize_img(char(strcat(temp, sub,'/cat12/mri/smwp1T1_orig_defaced.nii')), [8 8 8], nan(2,3))
            resize_img(char(strcat(temp, sub,'/cat12/mri/smwp2T1_orig_defaced.nii')), [8 8 8], nan(2,3))
            resize_img(char(strcat(temp, sub,'/cat12/mri/swmT1_orig_defaced.nii')), [8 8 8], nan(2,3))

            system(strcat({'rm '}, temp, sub, '/cat12/T1_orig_defaced.nii'));
            system(strcat({'cd '}, temp, sub, {'/; zip -r '}, sub, '_cat12.zip cat12/ -x "*.DS_Store"'));
            system(strcat({'cp '}, temp, sub, '/', sub, {'_cat12.zip '}, sub_folder, '/'));        
        else
            fileID = fopen(strcat(scanDir,'/errors.txt'),'a');
            fprintf(fileID,'%s: /T1/T1_orig_defaced.nii.gz does not exist\n', sub);
            fclose(fileID);
        end
        system(strcat({'rm -R '}, temp, sub));

    catch ME
        fileID = fopen(strcat(scanDir,'/errors.txt'),'a');
        fprintf(fileID,'%s: %s\n', sub, string(ME.message));
        fclose(fileID);
    end 
end

% shut down cluster and remove all jobs created with profile local
p = gcp;
delete(p);
myCluster = parcluster('local');
delete(myCluster.Jobs)
system(char(strcat({'rm -rf '}, temppath)));
exit
