#!/bin/bash

# ========================================================================
# ======= 2020 release: download t1 Files listed in basket 2005558 =======
# ========================================================================

# set working directory
cd /slow/projects/ukb_brainage/

# descend into basket directory and create file with .bulk data links
cd data/basket/20200217_2005558/data/
wget -O ukbconv -nd biobank.ndph.ox.ac.uk/ukb/util/ukbconv
./ukbconv ukb40487.enc_ukb bulk -o"t1w" -i<(echo '20252_2_0')

# get some stats
cat t1w.bulk | wc -l # how many files listed? - 42184
awk '{ print $1}' t1w.bulk | sort | uniq | wc -l # how many files of different people? - 40682
awk '$2=="20252_2_0" { count++; next} END { print count}' t1w.bulk # # how many files of field 20252_2_0? - 40681 
awk 'NR==FNR && $2=="20252_2_0" { ID[$1]; next} NR==FNR { next } $2=="20252_3_0" && !($1 in ID) { print }' t1w.bulk t1w.bulk # one subject only has repeat-imaging visit (20252_3_0)

# create destination folder, extract file links for initial imaging visit, and copy download key
cd - 
mkdir -p data/t1w/r2020
awk '$2=="20252_2_0" { print }' data/basket/20200217_2005558/data/t1w.bulk > data/t1w/r2021/t1w.bulk
cp data/basket/20200217_2005558/data/*.key data/t1w/r2021/

# descend into download folder and prepare download
cd data/t1w/r2020/
wget -O ukbfetch -nd biobank.ndph.ox.ac.uk/showcase/util/ukbfetch
chmod 700 ukbfetch
key=$(ls *.key)
ntotal=$(cat t1w.bulk | wc -l)

# download data
# load files in line x to x+500 (8 instances in parallel)
# -s Flag for starting line in t1w.bulk
# -o Flag for out log file name
N=8
(for i in $(eval echo "{1..$ntotal..500}"); do 
   ((j=j%N)); ((j++==0)) && wait
   ./ukbfetch -a${key} -bt1w.bulk -s${i} -m500 -ot1w.${i} &
done)

# have all files been downloaded?
# yep! - otherwise store missing files in new .bulk file and repeat download
awk 'NR==FNR { file[$1]; ndownloaded++; next } { ntotal++ } !($1 in file) { print }  END { print ntotal, ndownloaded }
    ' <(ls *.zip | sed "s%_20252_2_0.zip%%g") t1w.bulk 

# move t1w files into subject folders
shopt -s nullglob
for filename in *.zip; do
    curfile=$(basename $filename)
    curfile_first=${curfile:0:7}
    mkdir -p $curfile_first
    mv $filename $curfile_first
done

# run cat12 preprocessing
cd -
/opt/matlab/bin/matlab -nodesktop -nodisplay -r "scanDir = 'data/t1w/r2020/'; filetype = '_20252_2_0.zip'; spmPath = '/fast/software/matlab/spm12/'; ncores = 50; run code/mri/cat12.m"

# have all files been successfully processed? - Nope
find data/t1w/r2020 -name "*20252_2_0.zip" | wc -l
find data/t1w/r2020 -name "*cat12.zip" | wc -l

# which files have not been processed successfully? 
subs=$(find data/t1w/r2020 -mindepth 1 -maxdepth 1 -type d -not -exec sh -c 'ls -1 "{}"|egrep -i -q ".*cat12.zip"' ';' -print | sed "s%.*/%%g")
subs=($subs)
echo ${subs[@]}

# check whether files are "unusable"
> data/t1w/r2020/unusable.txt
for i in ${subs[@]}; do
    if [ $(unzip -l data/t1w/r2020/${i}/${i}_20252_2_0.zip | grep unusable | wc -l) == 1 ]; then
        echo "$i is unusable"; echo $i >> data/t1w/r2020/unusable.txt
    else
        echo "$i is usable"
    fi
done

# run Matlab script once again
/opt/matlab/bin/matlab -nodesktop -nodisplay -r "scanDir = 'data/t1w/r2020/'; filetype = '_20252_2_0.zip'; spmPath = '/fast/software/matlab/spm12/'; ncores = 50; run code/mri/cat12.m"

# change access rights
chmod 770 */*cat12.zip

# ========================================================================
# ======= 2021 release: download T1 files listed in basket 2007685 =======
# ========================================================================

# set working directory
cd /slow/projects/ukb_brainage/

# descend into basket directory and create file with .bulk data links
cd data/basket/20210205_2007685/data/
./ukbconv ukb45233.enc_ukb bulk -o"t1w" -i<(echo '20252_2_0')

# get some stats
cat t1w.bulk | wc -l # how many files listed? - 47248
awk '{ print $1}' t1w.bulk | sort | uniq | wc -l # how many files of different people? - 44178
awk '$2=="20252_2_0" { count++; next} END { print count}' t1w.bulk # how many files of field 20252_2_0? - 44175
awk 'NR==FNR && $2=="20252_2_0" { ID[$1]; next} NR==FNR { next } $2=="20252_3_0" && !($1 in ID) { print }' t1w.bulk t1w.bulk # three subjects only have repeat-imaging visit (20252_3_0)

# get list of incremental files
cd -
r2020=$(find data/t1w/r2020/* -name "*_20252_2_0.zip" | sed "s%.*/%%g" | sed "s%_20252_2_0.zip%%g") # r2020=$(find data/t1w/00_T1w_wave1/* data/t1w/00_T1w_wave2/* data/t1w/00_T1w_wave2_suppl/* -name "*_20252_2_0.zip" | sed "s%.*/%%g" | sed "s%_20252_2_0.zip%%g")
cd data/basket/20210205_2007685/data/
awk 'NR==FNR { file[$1]; next } !($1 in file) && $2=="20252_2_0" { print }' <(echo "${r2020}") t1w.bulk > t1w.incremental.bulk

# create destination folder, extract file links for initial imaging visit, and copy download key
cd - 
mkdir -p data/t1w/r2021
\cp data/basket/20210205_2007685/data/t1w.incremental.bulk data/t1w/r2021/
\cp data/basket/20210205_2007685/data/*.key data/t1w/r2021/

# descend into download folder and prepare download
cd data/t1w/r2021/
wget -O ukbfetch -nd biobank.ndph.ox.ac.uk/showcase/util/ukbfetch
chmod 700 ukbfetch
key=$(ls *.key)
ntotal=$(cat t1w.incremental.bulk | wc -l)

# download data
# load files in line x to x+500 (8 instances in parallel)
# -s Flag for starting line in t1w.bulk
# -o Flag for out log file name
N=8
(for i in $(eval echo "{1..$ntotal..500}"); do 
   ((j=j%N)); ((j++==0)) && wait
   ./ukbfetch -a${key} -bt1w.incremental.bulk -s${i} -m500 -ot1w.${i} &
done)

# have all files been downloaded?
# yep! - otherwise store missing files in new .bulk file and repeat download
awk 'NR==FNR { file[$1]; ndownloaded++; next } { ntotal++ } ($1 in file) { print }  END { print ntotal, ndownloaded }
    ' <(ls *.zip | sed "s%_20252_2_0.zip%%g") t1w.incremental.bulk 

# move t1w files into subject folders
shopt -s nullglob
for filename in *.zip; do
    curfile=$(basename $filename)
    curfile_first=${curfile:0:7}
    mkdir -p $curfile_first
    mv $filename $curfile_first
done

# run cat12 preprocessing
cd -
/opt/matlab/bin/matlab -nodesktop -nodisplay -r "scanDir = 'data/t1w/r2021/'; filetype = '_20252_2_0.zip'; spmPath = '/fast/software/matlab/spm12/'; ncores = 50; run code/mri/cat12.m"

# have all files been successfully processed? - Nope
find data/t1w/r2021 -name "*20252_2_0.zip" | wc -l
find data/t1w/r2021 -name "*cat12.zip" | wc -l

# which files have not been processed successfully? 
subs=$(find data/t1w/r2021 -mindepth 1 -maxdepth 1 -type d -not -exec sh -c 'ls -1 "{}"|egrep -i -q ".*cat12.zip"' ';' -print | sed "s%.*/%%g")
subs=($subs)
echo ${subs[@]}

# check whether files are "unusable"
> data/t1w/r2021/unusable.txt
for i in ${subs[@]}; do
    if [ $(unzip -l data/t1w/r2021/${i}/${i}_20252_2_0.zip | grep unusable | wc -l) == 1 ]; then
        echo "$i is unusable"; echo $i >> data/t1w/r2021/unusable.txt
    else
        echo "$i is usable"
    fi
done

# run Matlab script once again
/opt/matlab/bin/matlab -nodesktop -nodisplay -r "scanDir = 'data/t1w/r2021/'; filetype = '_20252_2_0.zip'; spmPath = '/fast/software/matlab/spm12/'; ncores = 50; run code/mri/cat12.m"

# change access rights
chmod 770 */*cat12.zip

# =======================================================================
# ======= 2022 release: get additional T1 files in basket 2016290 =======
# =======================================================================
# !! no additional files for the initial imaging visit

# set working directory
cd /slow/projects/ukb_brainage/

# descend into basket directory and create file with .bulk data links
cd data/basket/20220914_2016290/data/
./ukbconv ukb669463.enc_ukb bulk -o"t1w" -i<(echo '20252_2_0')

# get some stats
cat t1w.bulk | wc -l # how many files listed? - 49121
awk '{ print $1}' t1w.bulk | sort | uniq | wc -l # how many files of different people? - 44181
awk '$2=="20252_2_0" { count++; next} END { print count}' t1w.bulk # how many files of field 20252_2_0? - 44165
awk 'NR==FNR && $2=="20252_2_0" { ID[$1]; next} NR==FNR { next } $2=="20252_3_0" && !($1 in ID) { print }' t1w.bulk t1w.bulk # sixteen subjects only have repeat-imaging visit (20252_3_0)

# get list of incremental files
cd -
r2021=$(find data/t1w/r2020/* data/t1w/r2021/* -name "*_20252_2_0.zip" | sed "s%.*/%%g" | sed "s%_20252_2_0.zip%%g") # r2021=$(find data/t1w/00_T1w_wave1/* data/t1w/00_T1w_wave2/* data/t1w/00_T1w_wave2_suppl/* data/t1w/00_T1w_wave3/* -name "*_20252_2_0.zip" | sed "s%.*/%%g" | sed "s%_20252_2_0.zip%%g")
cd data/basket/20220914_2016290/data/
awk 'NR==FNR { file[$1]; next } !($1 in file) && $2=="20252_2_0" { print }' <(echo "${r2021}") t1w.bulk # no additional files

# ====================================================================================
# ======= 2022 release: get T1 files of repeat-imaging-visit in basket 2016290 =======
# ====================================================================================

# set working directory
cd /slow/projects/ukb_brainage/

# descend into basket directory and count repeat imaging files
cd data/basket/20220914_2016290/data/
awk '$2=="20252_3_0" { count++; next} END { print count}' t1w.bulk # - 4956

# create destination folder, extract file links for initial imaging visit, and copy download key
cd - 
mkdir -p data/t1w/r2022_retest
awk '$2=="20252_3_0" { print }' data/basket/20220914_2016290/data/t1w.bulk > data/t1w/r2022_retest/t1w.bulk
\cp data/basket/20220914_2016290/data/*.key data/t1w/r2022_retest/

# descend into download folder and prepare download
cd data/t1w/r2022_retest/
wget -O ukbfetch -nd biobank.ndph.ox.ac.uk/showcase/util/ukbfetch
chmod 700 ukbfetch
key=$(ls *.key)
ntotal=$(cat t1w.bulk | wc -l)

# download data
# load files in line x to x+500 (8 instances in parallel)
# -s Flag for starting line in t1w.bulk
# -o Flag for out log file name
N=8
(for i in $(eval echo "{1..$ntotal..500}"); do 
   ((j=j%N)); ((j++==0)) && wait
   ./ukbfetch -a${key} -bt1w.bulk -s${i} -m500 -ot1w.${i} &
done)

# have all files been downloaded?
# yep! - otherwise store missing files in new .bulk file and repeat download
awk 'NR==FNR { file[$1]; ndownloaded++; next } { ntotal++ } ($1 in file) { print }  END { print ntotal, ndownloaded }
    ' <(ls *.zip | sed "s%_20252_3_0.zip%%g") t1w.bulk 

# move t1w files into subject folders
shopt -s nullglob
for filename in *.zip; do
    curfile=$(basename $filename)
    curfile_first=${curfile:0:7}
    mkdir -p $curfile_first
    mv $filename $curfile_first
done

# run cat12 preprocessing
cd -
/opt/matlab/bin/matlab -nodesktop -nodisplay -r "scanDir = 'data/t1w/r2022_retest/'; filetype = '_20252_3_0.zip'; spmPath = '/fast/software/matlab/spm12/'; ncores = 50; run code/mri/cat12.m"

# have all files been successfully processed? - Nope
find data/t1w/r2022_retest -name "*20252_3_0.zip" | wc -l
find data/t1w/r2022_retest -name "*cat12.zip" | wc -l

# which files have not been processed successfully? 
subs=$(find data/t1w/r2022_retest -mindepth 1 -maxdepth 1 -type d -not -exec sh -c 'ls -1 "{}"|egrep -i -q ".*cat12.zip"' ';' -print | sed "s%.*/%%g")
subs=($subs)
echo ${subs[@]}

# check whether files are "unusable"
> data/t1w/r2022_retest/unusable.txt
for i in ${subs[@]}; do
    if [ $(unzip -l data/t1w/r2022_retest/${i}/${i}_20252_2_0.zip | grep unusable | wc -l) == 1 ]; then
        echo "$i is unusable"; echo $i >> data/t1w/r2022_retest/unusable.txt
    else
        echo "$i is usable"
    fi
done

# run Matlab script once again
/opt/matlab/bin/matlab -nodesktop -nodisplay -r "scanDir = 'data/t1w/r2020/'; filetype = '_20252_2_0.zip'; spmPath = '/fast/software/matlab/spm12/'; ncores = 50; run code/mri/cat12.m"

# change access rights
chmod 770 data/t1w/r2022_retest/*/*cat12.zip

# ==================================================================
# ======= Download FreeSurfer Files listed in basket 2013814 =======
# ==================================================================

# set working directory
cd /slow/projects/ukb_brainage/

# descend into basket directory and create file with .bulk data links
cd data/basket/20210703_2013814/data/
./ukbconv ukb40487.enc_ukb bulk -o"surface" -i<(echo '20263')
cd - 

# create destination folder, extract file links for initial imaging visit, and copy download key
mkdir -p data/surface
awk '$2=="20263_2_0" { print }' data/basket/20210703_2013814/data/surface.bulk > data/surface/surface_45k.bulk

# how many items listed? - 46,484 (46,484 last updated 28 Apr 2021: https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20263)
cat data/surface/surface_45k.bulk | wc -l

# how many files of different people? 43,140 (43,140 last updated 28 Apr 2021: https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20263)
awk '{ print $1}' data/surface/surface_45k.bulk | sort | uniq | wc -l

# how many files of field 2025632_0? 43,075
awk '$2=="20263_2_0" { count++; next} END { print count}' data/surface/surface_45k.bulk

# who misses 20263_2_0? - 65 subjects
awk 'NR==FNR && $2=="20263_2_0" { ID[$1]; next} NR==FNR { next } $2=="20263_3_0" && !($1 in ID) { print }' data/surface/surface_45k.bulk data/surface/surface_45k.bulk | wc -l

# create list of subjects for download (2_0 imaging visit)
awk '$2=="20263_2_0" { print; next} { next }' data/surface/surface_45k.bulk > data/surface/surface_45k_batch1.bulk

# create list of retest files
awk '$2=="20263_3_0" { print } ' data/surface/surface_45k.bulk > data/surface/surface_45k_batch1_retest.bulk
cat data/surface/surface_45k_batch1_retest.bulk | wc -l

# create download folder and change directory
mkdir -p data/surface/00_batch1
cd data/surface/00_batch1
xwget -O ukbfetch -nd biobank.ndph.ox.ac.uk/showcase/util/ukbfetch
chmod 700 ukbfetch
cp data/basket/20210703_2013814/data/*.key data/surface/00_batch1
cp ../surface_45k_batch1.bulk .

# start download (6 simultaneously) + one for retest
cat surface_45k_batch1.bulk | wc -l

tmux new -s surface_downl1
key=$(ls *.key)
for i in {1..8000..500}; do 
./ukbfetch -a${key} -bsurface_45k_batch1.bulk -s${i} -m500 -osurface_${i}
done

tmux new -s surface_downl2
key=$(ls *.key)
for i in {8001..16000..500}; do 
./ukbfetch -a${key} -bsurface_45k_batch1.bulk -s${i} -m500 -osurface_${i}
done

tmux new -s surface_downl3
key=$(ls *.key)
for i in {16001..24000..500}; do 
./ukbfetch -a${key} -bsurface_45k_batch1.bulk -s${i} -m500 -osurface_${i}
done

tmux new -s surface_downl4
key=$(ls *.key)
for i in {24001..32000..500}; do 
./ukbfetch -a${key} -bsurface_45k_batch1.bulk -s${i} -m500 -osurface_${i}
done

tmux new -s surface_downl5
key=$(ls *.key)
for i in {32001..40000..500}; do 
./ukbfetch -a${key} -bsurface_45k_batch1.bulk -s${i} -m500 -osurface_${i}
done

tmux new -s surface_downl6
key=$(ls *.key)
for i in {40001..43075..500}; do 
./ukbfetch -a${key} -bsurface_45k_batch1.bulk -s${i} -m500 -osurface_${i}
done

# all files downloaded? nope.
cd /slow/projects/ukb_brainage/data/surface/00_batch1
awk 'NR==FNR { file[$1]; count1++; next } { count2++ } !($1 in file) { print } END { print count1, count2 }' <(find . -name "*.zip" | sed "s%_20263_2_0.zip%%g" | sed "s%.*/%%g") surface_45k_batch1.bulk
awk 'NR==FNR { file[$1]; count1++; next } { count2++ } !($1 in file) { print } END { print count1, count2 }' <(find . -name "*.zip" | sed "s%_20263_2_0.zip%%g" | sed "s%.*/%%g") surface_45k_batch1.bulk > surface_45k_batch2.bulk

# download missing files
cat surface_45k_batch2.bulk | wc -l

tmux new -s surface_downl1
key=$(ls *.key)
for i in {1..98..500}; do 
./ukbfetch -a${key} -bsurface_45k_batch2.bulk -s${i} -m500 -osurface_${i}
done

# all files downloaded? yep.
cd /slow/projects/ukb_brainage/data/surface/00_batch1
awk 'NR==FNR { file[$1]; count1++; next } { count2++ } !($1 in file) { print } END { print count1, count2 }' <(find . -name "*.zip" | sed "s%_20263_2_0.zip%%g" | sed "s%.*/%%g") surface_45k_batch1.bulk

# move into folder 
shopt -s nullglob
for filename in *.zip; do
    curfile=$(basename $filename)
    curfile_first=${curfile:0:7}
    mkdir -p $curfile_first
    mv $filename $curfile_first
done
chmod -R 770 *
