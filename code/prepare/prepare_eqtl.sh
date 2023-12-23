#!/bin/bash

# =================================
# === prepare GTEx eQTL mapping ===
# =================================

# enter working directory
cd /home/groups/markett/ukb_brainage # cd /slow/projects/ukb_brainage

# create conda environment
if [ ! -d "envs/eqtl" ]; then
  conda env create --file envs/eqtl.yml -p envs/eqtl
fi
conda activate envs/eqtl

# add gene biotype
mkdir /home/groups/markett/software/gtex # mkdir /fast/software/gtex
cd /home/groups/markett/software/gtex # cd /fast/software/gtex

# download refseq identifiers
wget https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/README_eQTL_v8.txt
wget https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar
wget https://storage.googleapis.com/gtex_analysis_v8/reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz
wget https://storage.googleapis.com/gtex_analysis_v8/multi_tissue_qtl_data/GTEx_Analysis_v8.metasoft.txt.gz
tar -xvf GTEx_Analysis_v8_eQTL.tar

# select multi-tissue eqtls if data is available for at least 10 tissues (92% of all eqtls) and if mval is >= 0.9 in more than 50% of available tissues
# a very small fraction of remaining entries (0.1%) shows PVALUE_RE2 > 5E-8
# set inclusion criterion to PVALUE_RE2 < 5e-8 (add 0 to p-values (+0); noticed that 'e' scientific notation - but not 'E' - appears to make awk identify numbers smaller than 1e-308 as strings. so that they are not adequately dealt with)
awk -F'\t' 'NR==1 { print; next }
	{ count=0; for (i=66;i<144;i++) {if ($i>=0.9 && $i!="NA") {count++}} }
	$2>=10 && count>=0.5*$2 && $9+0<5E-8' OFS='\t' <(gzip -dc "GTEx_Analysis_v8.metasoft.txt.gz") > "GTEx_Analysis_v8.metasoft.selected.txt"
gzip -f "GTEx_Analysis_v8.metasoft.selected.txt"

# change access rights
chmod -R 770 *

# create symbolic link in project folder
cd /home/groups/markett/ukb_brainage # cd /slow/projects/ukb_brainage
ln -s /home/groups/markett/software/gtex data/gtex # ln -s /fast/software/gtex data/gtex

