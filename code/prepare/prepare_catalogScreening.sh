#!/bin/bash

# ============================
# === prepare GWAS catalog ===
# ============================

# set working directory
mkdir -p /fast/software/gwas_catalog
cd /fast/software/gwas_catalog

# download most recent GWAS catalog
# wget -P data/ --content-disposition https://www.ebi.ac.uk/gwas/api/search/downloads/full
wget --content-disposition https://www.ebi.ac.uk/gwas/api/search/downloads/full

# create symbolic link in project folder
cd /slow/projects/ukb_brainage/data
ln -s /home/groups/markett/software/gwas_catalog gwas_catalog