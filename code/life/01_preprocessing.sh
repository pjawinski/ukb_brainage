#!/bin/bash

# Contact: Philippe Jawinski (philippe.jawinski@hu-berlin.de)

# ============================================
# === get list of files for pre-processing ===
# ============================================

# access folder
cd /data/p_02330/

# move imaging files to 01_T1w
mkdir -p T1w
mv LI* T1w

# get list of files
filelist=$(find . -name *.nii.gz)

# how many files? - 3237
echo "$filelist" | wc -l

# how many subjects? - 2651
echo "$filelist" | sed s%./T1w/%%g | sed s%[/].*%%g | sort -u | wc -l

# how many different file types? - 11
filetype=$(echo "$filelist"  | sed s%.*[/]%%g | sort -u)
echo "$filetype"
echo "$filetype" | wc -l

# how many files of each type?
for i in $filetype; do
    count=$(echo "$filelist" | grep $i | wc -l)
    echo "$i - $count"
done

# get list of subjects with all recordings in one row
subs=$(echo "$filelist" | sort | awk '{sic=$1;
    file=$1;
    gsub(/.[/]T1w[/]/, "", sic);
    gsub(/[/].*/, "", sic);
    gsub(/.*[/]/, "", file)}

    sic in subj {
    subj[sic]=subj[sic]"\t"file;
    next}

    {subj[sic]=sic"\t"file;
    next}

    END { for (key in subj) { print subj[key] }}' OFS="\t")

# maximum number of recordings? - 3 (ncols-1)
echo "$subs" | awk 'BEGIN { ncols = 0 } NF > ncols { ncols = NF} END { print ncols}'

# if a subject has more than 1 recording, what types of recordings occur together?
echo "$subs" | awk '{ print $2, $3, $4 }' | sort -u | awk 'NF > 1 { print }'

# for each subject select one recording according to the priority shown below:
prio=$(awk '$2=="S2_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz" || $3=="S2_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz" || $4=="S2_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz" { print $1, "S2_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz"; next }
     $2=="S2_MPRAGE_ADNI_32Ch_PAT2_2.99.nii.gz" || $3=="S2_MPRAGE_ADNI_32Ch_PAT2_2.99.nii.gz" || $4=="S2_MPRAGE_ADNI_32Ch_PAT2_2.99.nii.gz" { print $1, "S2_MPRAGE_ADNI_32Ch_PAT2_2.99.nii.gz"; next }
     $2=="S3_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz" || $3=="S3_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz" || $4=="S3_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz" { print $1, "S3_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz"; next }
     $2=="S4_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz" || $3=="S4_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz" || $4=="S4_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz" { print $1, "S4_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz"; next }
     $2=="S5_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz" || $3=="S5_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz" || $4=="S5_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz" { print $1, "S5_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz"; next }
     $2=="S8_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz" || $3=="S8_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz" || $4=="S8_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz" { print $1, "S8_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz"; next }
     $2=="S9_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz" || $3=="S9_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz" || $4=="S9_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz" { print $1, "S9_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz"; next }
     $2=="S2_MPRAGE_ADNI_32Ch_2.98.nii.gz" || $3=="S2_MPRAGE_ADNI_32Ch_2.98.nii.gz" || $4=="S2_MPRAGE_ADNI_32Ch_2.98.nii.gz" { print $1, "S2_MPRAGE_ADNI_32Ch_2.98.nii.gz"; next }
     $2=="S3_MPRAGE_ADNI_32Ch_2.98.nii.gz" || $3=="S3_MPRAGE_ADNI_32Ch_2.98.nii.gz" || $4=="S3_MPRAGE_ADNI_32Ch_2.98.nii.gz" { print $1, "S3_MPRAGE_ADNI_32Ch_2.98.nii.gz"; next }
     $2=="S3_MPRAGE_ADNI_2.98.nii.gz" || $3=="S3_MPRAGE_ADNI_2.98.nii.gz" || $4=="S3_MPRAGE_ADNI_2.98.nii.gz" { print $1, "S3_MPRAGE_ADNI_2.98.nii.gz"; next }
     $2=="S4_MPRAGE_ADNI_32Ch_2.98.nii.gz" || $3=="S4_MPRAGE_ADNI_32Ch_2.98.nii.gz" || $4=="S4_MPRAGE_ADNI_32Ch_2.98.nii.gz" { print $1, "S4_MPRAGE_ADNI_32Ch_2.98.nii.gz"}' OFS='\t' <(echo "$subs"))

# all subjects preserved? - yep (2651)!
echo "$prio" | wc -l

# all "S2_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz" recordings preserved? - yep (2067)!
grep S2_MPRAGE_ADNI_32Ch_PAT2_2.98.nii.gz <<< "$prio" | wc -l

# create subs_to_process.txt
prio=($prio)
rm -f subs_to_process.txt
touch subs_to_process.txt
for i in $(seq 1 2 ${#prio[@]}); do
    sic=${prio[$i-1]}
    file=${prio[$i]}
    awk -v sic=$sic -v file=$file 'index($1, sic) != 0 && index($1, file) != 0 { print; exit }' <(echo "$filelist") >> "subs_to_process.txt"
done

# check number of rows - correct (2651)
cat subs_to_process.txt | wc -l

# create backup_file
cp subs_to_process_full.txt

# ===============================================
# === Access HCP server and run preprocessing ===
# ===============================================

# access hcp
ssh comps06h12
tmux
cd /data/p_02330

# start Matlab
matlab -nodesktop -nodisplay 

%% set Matlab toolbox paths
addpath /data/p_02330/functions/
addpath /data/p_02330/toolbox/spm12/

% set paths and read subs_to_process.txt
path = strcat(pwd);
temppath = strcat(path,'/temp');
system(char(strcat({'mkdir -p '}, temppath)));
fid = fopen(strcat(path,'/subs_to_process.txt')); % subjects = textread(strcat(path,'subs_to_process.txt'));                                     
subjects = textscan(fid, '%s');                         
fclose(fid);
subjects = subjects{1};

% set up Matlab job (cat12.5 build 1364)
matlabbatch{1}.spm.tools.cat.estwrite.nproc = 0;
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP = 1070;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.darteltpm = {'/data/p_02330/toolbox/spm12/toolbox/cat12/templates_1.50mm/Template_1_IXI555_MNI152.nii'}; % matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.darteltpm = {'/Applications/MATLAB_R2018a.app/toolbox/spm12/toolbox/cat12/templates_1.50mm/Template_1_IXI555_MNI152.nii'};
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shootingtpm = {'/data/p_02330/toolbox/spm12/toolbox/cat12/templates_1.50mm/Template_0_IXI555_MNI152_GS.nii'}; % matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shootingtpm = {'/Applications/MATLAB_R2018a.app/toolbox/spm12/toolbox/cat12/templates_1.50mm/Template_0_IXI555_MNI152_GS.nii'};
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

% set maximum number of compute threads
maxNumCompThreads(1);

% start parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = 20;
saveProfile(myCluster);
parpool('local',20);

% submit constant variables to workers
matlabbatch_constant = parallel.pool.Constant(matlabbatch);
subjects_constant = parallel.pool.Constant(subjects);
dir_constant = parallel.pool.Constant(path);
tempdir_constant = parallel.pool.Constant(temppath);

% initiate spm
spm('Defaults','fMRI');
spm_jobman('initcfg');

% test without loop
% exist(strcat('/home/jawinskp/toolbox/spm12/toolbox/cat12/templates_1.50mm/Template_1_IXI555_MNI152.nii'),'file')
% exist(strcat('/home/jawinskp/toolbox/spm12/toolbox/cat12/templates_1.50mm/Template_0_IXI555_MNI152_GS.nii'),'file')
idx=1;
dir = path;
tempdir = temppath;
parts = strsplit(string(subjects(idx)), {'./', '.gz'});
element = strcat('/', parts(2));
parts = strsplit(element, {'/T1w', '/NIFTI/'});
element_folder = strcat('/T1w', parts(2), '/NIFTI');
parts = strsplit(element, {'/NIFTI/', '.gz'});
element_nifti = parts(2);
matlabbatch_loop = matlabbatch;

parfor idx = 1:numel(subjects)
    
    dir = dir_constant.Value;
    tempdir = tempdir_constant.Value;
    parts = strsplit(string(subjects_constant.Value(idx)), {'./', '.gz'});
    element = strcat('/', parts(2));
    parts = strsplit(element, {'/T1w', '/NIFTI/'});
    element_folder = strcat('/T1w', parts(2), '/NIFTI');
    parts = strsplit(element, {'/NIFTI/', '.gz'});
    element_nifti = parts(2);
    matlabbatch_loop = matlabbatch_constant.Value;

    try
        system(strcat({'mkdir -p '}, tempdir, element_folder, {'/; cp -R '}, path, element_folder, {'/. '}, tempdir, element_folder));
        system(strcat({'gunzip -fk '}, tempdir, element));    

        if exist(strcat(tempdir, element_folder, '/', element_nifti), 'file') == 2
            system(strcat({'mkdir -p '}, tempdir, element_folder, '/cat12'));
            system(strcat({'mv '}, tempdir, element_folder, '/', element_nifti, {' '}, tempdir, element_folder, '/cat12'));

            matlabbatch_loop{1}.spm.tools.cat.estwrite.data = {char(strcat(tempdir, element_folder, '/cat12/', element_nifti))};
            matlabbatch_loop{2}.spm.spatial.smooth.data =  {char(strcat(tempdir, element_folder, '/cat12/mri/mwp1', element_nifti))
                                                        char(strcat(tempdir, element_folder, '/cat12/mri/mwp2', element_nifti))
                                                        char(strcat(tempdir, element_folder, '/cat12/mri/wm', element_nifti))};
            spm_jobman('run',matlabbatch_loop);

            resize_img(char(strcat(tempdir, element_folder,'/cat12/mri/smwp1', element_nifti)), [8 8 8], nan(2,3))
            resize_img(char(strcat(tempdir, element_folder,'/cat12/mri/smwp2', element_nifti)), [8 8 8], nan(2,3))
            resize_img(char(strcat(tempdir, element_folder,'/cat12/mri/swm', element_nifti)), [8 8 8], nan(2,3))

            system(strcat({'rm '}, tempdir, element_folder, '/cat12/', element_nifti));
            char_element_nifti = char(element_nifti);
            system(strcat({'cd '}, tempdir, element_folder, {'/; zip -r '}, char_element_nifti(1:end-4), '_cat12.zip cat12/ -x "*.DS_Store"'));
            system(strcat({'mkdir -p '}, dir, '/cat12', element_folder, {'/; cp '}, tempdir, element_folder, '/', char_element_nifti(1:end-4), {'_cat12.zip '}, dir, '/cat12', element_folder, '/'));
        end
    
        charelement = char(element);
        system(strcat({'rm -R '}, tempdir, element_folder));    
        system(strcat('sed -i -e "s%', string(subjects_constant.Value(idx)), '%%g"', {' '}, dir, '/subs_to_process.txt')); % remove subj from to do list
        system(string(strcat('sed -i -e "/^[[:space:]]*$/d"', {' '}, dir, '/subs_to_process.txt'))); % remove nasty white space

    catch ME
    system(strcat({'echo "'}, string(strcat(element,{': '}, ME.message)), {'" >> '}, dir, '/error_messages.txt')); 
    end
end 

% shut down cluster and remove all jobs created with profile local
p = gcp;
delete(p);
myCluster = parcluster('local');
delete(myCluster.Jobs)

% exit matlab
exit

# ========================================================
# === get cat12 data required for brain age estimation ===
# ========================================================

# how many cat12.zip files have been created? - 2651 expected, 2650 observed
find cat12/ -name *cat12.zip | wc -l

# get list of subjects whose job execution failed - 1 subject with S2_MPRAGE_ADNI_32Ch_PAT2_2.98.nii (restarting job execution failed again)
case=$(head -1 error_messages.txt | sed s%/T1w/%%g | sed s%[/].*%%g)

# are there other recordings available for this subject? - Nope.
find T1w/ -name *.nii.gz | grep ${case}

# start Matlab and get cat12 meta infos
matlab -nodesktop -nodisplay 

% add folder with 'read_nifti' function
addpath /data/p_02330/functions/
addpath /data/p_02330/toolbox/spm12/

% set source and destination path
path = '/data/p_02330/cat12/'


% get list of _cat12.zip files
files = dir(strcat(path,'**/*_cat12.zip'));

% get individual identifier from files.name (numeric in UK Biobank)
IID = {};
for i = 1:length(files)
    parts = strsplit(string(files(i).folder), {'/'});
    IID{i,1} = parts(6);
end

% correct typos in meta.IID
for i = 1:size(meta.IID,1)
    if meta.IID{i,1} == "LI01927958"
        meta.IID{i,1} = "LI01927058";
    end
    if meta.IID{i,1} == "LI02033591"
        meta.IID{i,1} = "LI02022591";
    end
end


% unzip rsmwp[1/2]*.nii, read nifti, and create grey and white matter matrix
gm = zeros(length(files), 16128);
wm = zeros(length(files), 16128);

tic
for i = 1:length(files)
    system(char(strcat({'unzip -joqq '}, files(i).folder, '/', files(i).name, {' cat12/mri/rsmwp1*.nii -d '}, files(i).folder, '/')));
    filename = strcat('/rsmwp1', files(i).name(1:end-10), '.nii');
    nifti = read_nifti(strcat(files(i).folder, filename));
    gm(i,:) = nifti(:)';
    system(char(strcat({'rm -f '}, files(i).folder, '/rsmwp1*.nii'))); 
    filename = []; nifti = [];
    
    system(char(strcat({'unzip -joqq '}, files(i).folder, '/', files(i).name, {' cat12/mri/rsmwp2*.nii -d '}, files(i).folder, '/')));
    filename = strcat('/rsmwp2', files(i).name(1:end-10), '.nii');
    nifti = read_nifti(strcat(files(i).folder, filename));
    wm(i,:) = nifti(:)';
    system(char(strcat({'rm -f '}, files(i).folder, '/rsmwp2*.nii'))); 
    filename = []; nifti = [];

    fprintf(strcat(sprintf('%0.2f',i/length(files)*100), ' percent done.\n'));
end
toc

% have all files been successfully imported (no empty lines)?
sum(sum(gm, 2)~=0) % should be 2650
sum(sum(wm, 2)~=0) % should be 2650

% unzip cat12/report/*.mat and collect qualityratings and subjectmeasures
meta_ratings = zeros(length(files),9);

tic
for i = 1:length(files)
    system(char(strcat({'unzip -joqq '}, files(i).folder, '/', files(i).name, {' cat12/report/cat_*.mat -d '}, files(i).folder, '/')));
    logfilename = strcat('/cat_', files(i).name(1:end-10), '.mat');
    load(strcat(files(i).folder, logfilename)); 
    meta_ratings(i,:) = [S.subjectmeasures.vol_TIV, S.subjectmeasures.vol_abs_CGW(:,1:4) S.qualityratings.res_RMS S.qualityratings.NCR S.qualityratings.ICR S.qualityratings.IQR];
    clearvars S
    system(char(strcat({'rm -f '}, strcat(files(i).folder, logfilename))));
    fprintf(strcat(sprintf('%0.2f',i/length(files)*100), ' percent done.\n'));
end
toc

% have all files been successfully imported (no empty lines)? - Yep.
sum(sum(meta_ratings, 2)~=0) % should be 2650

% put into structure
meta = struct();
meta.files = files;
meta.ratings = meta_ratings;
meta.ratings_varnames = {'TIV', 'CSF', 'GM', 'WM', 'RES', 'RMS','NCR','ICR', 'IQR'};
meta.IID = IID;

% identify recordings with poor image quality rating (IQR > 3)
[TF, LTHRESH, UTHRESH, CENTER] = isoutlier(meta.ratings(:,9), 'quartiles');
sum(meta.ratings(:,9) > UTHRESH)
sum(meta.ratings(:,9) > 3)
meta.ratings(:,10) = meta.ratings(:,9) > 3;
meta.ratings_varnames = {'TIV', 'CSF', 'GM', 'WM', 'RES', 'RMS','NCR','ICR', 'IQR', 'IQR_poor'};

% save cat12 data required for brain age estimation
save(strcat('/data/p_02330/01_cat12_data.mat'), 'meta', 'gm', 'wm', '-v7.3')

