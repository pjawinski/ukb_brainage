#!/bin/bash

# =====================
# === install gwama ===
# =====================

# set working directory
mkdir /fast/software/gwama
cd /fast/software/gwama

# download mrmega
wget https://www.geenivaramu.ee/tools/GWAMA_v2.2.2.zip
wget https://www.geenivaramu.ee/tools/samples.zip
wget https://www.geenivaramu.ee/tools/dbsnp37.txt.gz
wget https://www.geenivaramu.ee/tools/MANH.R
wget https://www.geenivaramu.ee/tools/QQ.R
wget https://www.geenivaramu.ee/tools/SNPTEST2GWAMA.pl
wget https://www.geenivaramu.ee/tools/SNPTEST2_2_GWAMA.pl
wget https://www.geenivaramu.ee/tools/SNPTEST2.5_2_GWAMA.pl
wget https://www.geenivaramu.ee/tools/PLINK2GWAMA.pl

# install
unzip GWAMA_v2.2.2.zip 
make
ln -s /fast/software/gwama/GWAMA /fast/software/bin/GWAMA