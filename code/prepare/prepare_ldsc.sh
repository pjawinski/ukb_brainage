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
cd /home/groups/markett/ukb_brainage

# create conda environment
if [ ! -d "envs/ldsc" ]; then
  conda env create --file envs/ldsc.yml -p envs/ldsc
fi
conda activate envs/ldsc

# download LDSC
cd /home/groups/markett/software
git clone https://github.com/bulik/ldsc.git
cd ldsc

# get LD scores and hapmap3 snplist
mkdir -p /home/groups/markett/software/ldsc/resources
cd /home/groups/markett/software/ldsc/resources
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
wget http://ldsc.broadinstitute.org/static/media/w_hm3.noMHC.snplist.zip; cp $HOME/ldsc/resources/w_hm3.noMHC.snplist.zip .
unzip w_hm3.noMHC.snplist.zip

tar -jxvf eur_w_ld_chr.tar.bz2
bunzip2 w_hm3.snplist.bz2
chmod -R 770 *

old="w_hm3.snplist"
new="w_hm3.noMHC.snplist"
cmp --silent $old $new || echo "files are different"

awk 'NR==1 { next } NR==FNR { snp[NR]=$1"\t"$2"\t"$3; next} ($1"\t"$2"\t"$3 != snp[FNR]) { print }' $old $new 
awk 'NR==FNR { count1++; next} { count2++ } END { print count1, count2 }' $old $new 

# download files for partitioned heritability
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_baselineLD_v2.2_ldscores.tgz
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_weights_hm3_no_MHC.tgz
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_frq.tgz
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_cell_type_groups.tgz

tar zxvf 1000G_Phase3_baselineLD_v2.2_ldscores.tgz; mkdir baselineLD; mv baselineLD* baselineLD/
tar zxvf 1000G_Phase3_weights_hm3_no_MHC.tgz
tar zxvf 1000G_Phase3_frq.tgz
tar zxvf 1000G_Phase3_cell_type_groups.tgz
chmod -R 770 *

wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase1_cell_type_groups.tgz
tar zxvf 1000G_Phase1_cell_type_groups.tgz
chmod -R 770 *

# create symbolic links in bin folder
cd /home/groups/markett/software/bin # cd /fast/software/bin
ln -s /home/groups/markett/software/ldsc/munge_sumstats.py munge_sumstats.py
ln -s /home/groups/markett/software/ldsc/ldsc.py ldsc.py

# create symbolic link in project folder
cd /home/groups/markett/ukb_brainage/data
ln -s /home/groups/markett/software/ldsc ldsc



