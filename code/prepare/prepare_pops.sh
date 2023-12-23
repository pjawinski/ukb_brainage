#!/bin/bash

# ========================================
# === prepare POPS gene prioritization ===
# ========================================

# prepare environment
cd /slow/projects/ukb_brainage

# create conda 'default' project environment
mamba create -p envs/pops -c r -c conda-forge mamba python=3.7 scipy==1.5.2 pandas==1.0.5 numpy==1.19.5 scikit-learn==0.24.2 # mamba create -p envs/pops -c r -c conda-forge mamba python=3.7 scipy==1.7.1 pandas==1.3.4 numpy==1.21.2 scikit-learn==0.23.2

# save python environment in yml file
conda env export --no-builds -p envs/pops > envs/pops.yml

# add name and remove prefix
awk 'NR==1 { print "name: pops"; next } $1=="prefix:" { $2="pops"; print; next} { print }' envs/pops.yml > tmp && \mv tmp envs/pops.yml

# install software
cd /fast/software
git clone https://github.com/FinucaneLab/pops
cd pops

# run example
# munge features
conda activate /slow/projects/ukb_brainage/envs/pops
python munge_feature_directory.py \
 --gene_annot_path example/data/utils/gene_annot_jun10.txt \
 --feature_dir example/data/features_raw/ \
 --save_prefix example/data/features_munged/pops_features \
 --max_cols 500

# Run PoPS
python pops.py \
 --gene_annot_path example/data/utils/gene_annot_jun10.txt \
 --feature_mat_prefix example/data/features_munged/pops_features \
 --num_feature_chunks 2 \
 --magma_prefix example/data/magma_scores/PASS_Schizophrenia \
 --control_features_path example/data/utils/features_jul17_control.txt \
 --out_prefix example/out/PASS_Schizophrenia

# create symbolic links
ln -s /fast/software/pops/pops.py /fast/software/bin/pops.py
ln -s /fast/software/pops/munge_feature_directory.py /fast/software/bin/munge_feature_directory.py
chmod 750 *
ls -lh

# get article feature matrix
mkdir -p data/features_raw data/features_munged data/utils
wget -O data/features_raw/PoPS.features.txt.gz https://www.dropbox.com/sh/o6t5jprvxb8b500/AABK2h7UIw3K85oHE85Ey2ZNa/data/PoPS.features.txt.gz
wget -O data/utils/magma_0kb.genes.annot https://www.dropbox.com/sh/o6t5jprvxb8b500/AAAKoDqP6xdOmONoFLmQrNG8a/data/magma_0kb.genes.annot
wget -O data/utils/gene_loc.txt https://www.dropbox.com/sh/o6t5jprvxb8b500/AABtdAOOakFs6n0yW-9RUNR7a/data/gene_loc.txt
wget -O data/utils/control.features https://www.dropbox.com/sh/o6t5jprvxb8b500/AAAlq-Pux0eGdlIGB3IxPPj4a/data/control.features
chmod 750 data/features_raw/PoPS.features.txt.gz
chmod 750 data/utils/magma_0kb.genes.annot
chmod 750 data/utils/gene_loc.txt
cp example/data/utils/gene_annot_jun10.txt data/utils/
gunzip data/features_raw/PoPS.features.txt.gz

# munge files
python munge_feature_directory.py \
 --gene_annot_path data/utils/gene_annot_jun10.txt \
 --feature_dir data/features_raw/ \
 --save_prefix data/features_munged/pops_features \
 --max_cols 500

# =====================
# === install MAGMA ===
# =====================

# prepare environment
cd /fast/software

# download magma
mkdir -p magma
cd magma
wget https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.10.zip
unzip magma_v1.10.zip

# create symbolic link
ln -s /fast/software/magma/magma /fast/software/bin/magma

