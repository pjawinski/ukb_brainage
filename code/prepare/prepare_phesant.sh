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


# ================================
# === deprecated code snippets ===
# ================================

# R packages
R -e ' packageurls = c("https://cran.r-project.org/src/contrib/Archive/getopt/getopt_1.20.0.tar.gz",
"https://cran.r-project.org/src/contrib/Archive/optparse/optparse_1.3.2.tar.gz",
"https://cran.r-project.org/src/contrib/Archive/MASS/MASS_7.3-45.tar.gz",
"https://cran.r-project.org/src/contrib/Archive/lattice/lattice_0.20-31.tar.gz",
"https://cran.r-project.org/src/contrib/Archive/zoo/zoo_1.7-13.tar.gz",
"https://cran.r-project.org/src/contrib/Archive/lmtest/lmtest_0.9-34.tar.gz",
"https://cran.r-project.org/src/contrib/Archive/nnet/nnet_7.3-12.tar.gz",
"https://cran.r-project.org/src/contrib/Archive/magrittr/magrittr_1.5.tar.gz",
"https://cran.r-project.org/src/contrib/Archive/backports/backports_1.1.0.tar.gz",
"https://cran.r-project.org/src/contrib/Archive/checkmate/checkmate_1.8.2.tar.gz",
"https://cran.r-project.org/src/contrib/Archive/forestplot/forestplot_1.7.tar.gz"
"https://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.10.4.tar.gz")
for (url in packageurls) { 
install.packages(url, repos=NULL, type="source")
}
'

