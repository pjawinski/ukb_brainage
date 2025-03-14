#!/bin/bash

# ===================================
# === prepare ld score regression ===
# ===================================

# for some reason ldsc freezes at "Reading sumstats from [...]", so installed pandas 0.24.2 and numpy 1.16.5 (works on cluster1,so tried that on cluster3) 
# pip install 'pandas==0.24.2' 'numpy==1.16.5'

# set working directory
cd /slow/projects/ukb_brainage

# activate ldsc on cluster3
conda activate ldsc 

# save python environment in yml file
conda env export > envs/ldsc.yml

# add name and remove prefix
awk 'NR==1 { print "name: ldsc"; next } $1=="prefix:" { $2="ldsc"; print; next} { print }' envs/ldsc.yml > envs/ldsc.yml.tmp; \mv envs/ldsc.yml.tmp envs/ldsc.yml

# ===============
# === install ===
# ===============

# set working directory
cd /slow/projects/ukb_brainage

# create conda environment
if [ ! -d "envs/ldsc" ]; then
  conda env create --file envs/ldsc.yml -p envs/ldsc
fi
conda activate envs/ldsc

# download LDSC
cd /fast/software/
git clone https://github.com/bulik/ldsc.git
cd ldsc

# get LD scores and hapmap3 snplist
mkdir -p /fast/software/ldsc/resources
cd /fast/software/ldsc/resources
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
wget http://ldsc.broadinstitute.org/static/media/w_hm3.noMHC.snplist.zip # cp "$HOME/ldsc/resources/w_hm3.noMHC.snplist.zip" .
unzip w_hm3.noMHC.snplist.zip
tar -jxvf eur_w_ld_chr.tar.bz2
bunzip2 w_hm3.snplist.bz2
chmod -R 770 *

# download files for partitioned heritability
# Recommendations for cell-type specific analyses (readme_baseline_versions)
# 1. We recommend that for identifying critical tissues/cell-types via P-value of tau, it is best to use the baseline model, specifically baseline v1.2.
# 2. We recommend that for estimating heritability enrichment (i.e., %h2/%SNPs) of any annotation, including tissue-specific annotations, it is best to use baselineLD v2.2.    
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/readme_baseline_versions 
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_baseline_v1.2_ldscores.tgz 
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_baselineLD_v2.2_ldscores.tgz
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_weights_hm3_no_MHC.tgz
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_frq.tgz
mkdir 1000G_Phase3_baseline_v1.2; tar zxvf 1000G_Phase3_baseline_v1.2_ldscores.tar -C 1000G_Phase3_baseline_v1.2/
mkdir 1000G_Phase3_baseline_v1.2/baseline_v1.2/tmp; tar zxvf 1000G_Phase3_baseline_v1.2/baseline_v1.2/1000G_Phase3_baseline_v1.2_ldscores.tgz -C 1000G_Phase3_baseline_v1.2/baseline_v1.2/tmp/ # Caveat: .tgz within .tar file; extract!
cd 1000G_Phase3_baseline_v1.2/baseline_v1.2; files=$(ls); for i in $files; do cmp $i tmp/$i; done # all files identical
rm -f baseline* ; cd -; rm -rf 1000G_Phase3_baseline_v1.2/baseline_v1.2
mkdir 1000G_Phase3_baseline_v2.2; tar zxvf 1000G_Phase3_baseline_v2.2_ldscores.tar -C 1000G_Phase3_baseline_v2.2/
tar zxvf 1000G_Phase3_weights_hm3_no_MHC.tgz
tar zxvf 1000G_Phase3_frq.tgz
chmod -R 770 *

# download data for cell-type specific analysis
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/LDSC_SEG_ldscores/Multi_tissue_gene_expr_1000Gv3_ldscores.tgz
tar -xvzf Multi_tissue_gene_expr_1000Gv3_ldscores.tgz
sed 's%Multi_tissue_gene_expr_1000Gv3_ldscores/%/fast/software/ldsc/resources/Multi_tissue_gene_expr_1000Gv3_ldscores/%g' Multi_tissue_gene_expr.ldcts > Multi_tissue_gene_expr_fullpaths.ldcts
chmod -R 770 *

# download data for cell-type group analysis - 1000G Phase3
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/LDSC_SEG_ldscores/1000G_Phase3_cell_type_groups.tgz
tar -xvzf 1000G_Phase3_cell_type_groups.tgz
mv cell_type_groups 1000G_Phase3_cell_type_groups
chmod -R 770 1000G_Phase3_cell_type_groups
sed 's%Multi_tissue_gene_expr_1000Gv3_ldscores/%/fast/software/ldsc/resources/Multi_tissue_gene_expr_1000Gv3_ldscores/%g' 1000G_Phase3_cell_type_groups > Multi_tissue_gene_expr_fullpaths.ldcts
chmod -R 770 *

# create symbolic links in bin folder
cd /fast/software/bin # cd /fast/software/bin
ln -s /fast/software/ldsc/munge_sumstats.py munge_sumstats.py
ln -s /fast/software/ldsc/ldsc.py ldsc.py

# create symbolic link in project folder
cd /slow/projects/ukb_brainage/data
ln -s /fast/software/ldsc ldsc


