#!/bin/bash
# ===========================
# === prepare finemapping ===
# ===========================

# prepare FINEMAP
cd /fast/software
mkdir finemap; cd finemap
wget http://christianbenner.com/finemap_v1.4.2_x86_64.tgz
tar -xvzf finemap_v1.4.1_x86_64.tgz
ln -s /fast/software/finemap/finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64 /fast/software/bin/finemap

# prepare LDstore2
cd /fast/software
mkdir ldstore; cd ldstore
wget http://www.christianbenner.com/ldstore_v2.0_x86_64.tgz
tar -xvzf ldstore_v2.0_x86_64.tgz
ln -s /fast/software/ldstore/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 /fast/software/bin/ldstore

# prepare environment
cd /slow/projects/ukb_brainage

	# create conda 'default' project environment
	mamba create -p envs/finemap -c conda-forge pigz python=3.9 r-base r-bigsnpr r-data.table r-dplyr r-Rfast r-stringr r-susieR pandoc
	conda activate envs/finemap

	# save python environment in yml file
	conda env export -p envs/finemap > envs/finemap.yml

	# add name and remove prefix
	awk 'NR==1 { print "name: finemap"; next } $1=="prefix:" { $2="finemap"; print; next} { print }' envs/finemap.yml > tmp && \mv tmp envs/finemap.yml