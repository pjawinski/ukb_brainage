#!/bin/bash

# =================================
# === prepare phenome-wide scan ===
# =================================

# set working directory
cd /slow/projects/01_UKB/11b_brainage/

# create conda environment for phesant
conda create -p envs/phesant -c conda-forge python=3.9 r-base=4.1.1 r-optparse r-MASS r-lmtest r-nnet r-forestplot r-data.table r-dplyr r-htmlwidgets r-plotly r-ggrepel pandoc
conda activate envs/phesant

# save python environment in yml file
conda env export -p envs/phesant > envs/phesant.yml

# add name and remove prefix
awk 'NR==1 { print "name: phesant"; next } $1=="prefix:" { $2="phesant"; print; next} { print }' envs/phesant.yml > envs/phesant.yml.tmp; \mv envs/phesant.yml.tmp envs/phesant.yml

# clone phesant from GitHub
cd /fast/software
git clone https://github.com/MRCIEU/PHESANT.git
cd PHESANT

# run test testWAS (https://github.com/MRCIEU/PHESANT/tree/master/testWAS)
cd /fast/software/PHESANT/WAS
testDir="../testWAS/"

Rscript phenomeScan.r \
--phenofile="${testDir}data/phenotypes.csv" \
--traitofinterestfile="${testDir}data/exposure.csv" \
--variablelistfile="${testDir}variable-lists/outcome-info.tsv" \
--datacodingfile="${testDir}variable-lists/data-coding-ordinal-info.txt" \
--traitofinterest="exposure" \
--resDir="${testDir}results/" \
--userId="userId"

# shortcut
Rscript phenomeScan.r --test

# run part 1 of 3
Rscript phenomeScan.r \
--phenofile="${testDir}data/phenotypes.csv" \
--traitofinterestfile="${testDir}data/exposure.csv" \
--variablelistfile="${testDir}variable-lists/outcome-info.tsv" \
--datacodingfile="${testDir}variable-lists/data-coding-ordinal-info.txt" \
--traitofinterest="exposure" \
--resDir="${testDir}results/" \
--userId="userId" \
--partIdx=1 \
--numParts=3

# run results processing
cd ../resultsProcessing/

Rscript mainCombineResults.r \
--resDir="../testWAS/results/" \
--variablelistfile="../testWAS/variable-lists/outcome-info.tsv"

# prepare phenotype file in ukb basket (white-British only)
cd /slow/projects/01_UKB/11b_brainage/data/baskets/20200409_2007685/data
idfile='/fast/UK_Biobank/04_data_genetics_linux/00_script/release_feb2020/02_MRI_sample_wba.txt'
awk -F',' 'NR==FNR { id[$1]; next } FNR==1 { print } $1 in id { print }' OFS='\t' <(awk 'NR > 1 { print "\""$1"\"" }' $idfile) ukb41573.csv > ukb41573_imaging_wba.csv
awk -F',' 'NR==1 { gsub(/-/,"."); gsub(/[.]/,"_"); gsub(/,"/,",\"x"); print; next} NR > 1 { print }' ukb41573_imaging_wba.csv > ukb41573_imaging_wba_phesant.csv

# prepare phenotype file in ukb basket 2021 (multi-ancestry)
cd /fast/UK_Biobank/02_data_standard/20210205_2007685/data
./ukbconv ukb45233.enc_ukb csv
idfile='MRI_sample.txt'
awk -F',' 'NR==FNR { id[$1]; next } FNR==1 { print } $1 in id { print }' OFS='\t' <(awk 'NR > 1 { print "\""$1"\"" }' ${idfile}) ukb45233.csv > ukb45233_imaging.csv
awk -F',' 'NR==1 { gsub(/-/,"."); gsub(/[.]/,"_"); gsub(/,"/,",\"x"); print; next} NR > 1 { print }' ukb45233_imaging.csv > ukb45233_imaging_phesant.csv


