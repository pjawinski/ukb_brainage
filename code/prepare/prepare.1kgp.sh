#!/bin/bash

# ===================================
# === get 1k genomes project data ===
# ===================================

# download phase3 v5a
mkdir -p /fast/software/1kgp/v5a && cd /fast/software/1kgp/v5a
for i in {1..22}; do (
	echo starting with chromosome "${i}"
	wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr"${i}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
	wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr"${i}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi ) & 
done
wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1a.20130502.genotypes.vcf.gz &
wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1a.20130502.genotypes.vcf.gz.tbi &
wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chrY.phase3_integrated_v1a.20130502.genotypes.vcf.gz &
wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chrY.phase3_integrated_v1a.20130502.genotypes.vcf.gz.tbi
wait
chmod 550 ALL*
rm -f wget-log*

# keep meta data only (to speed up 1kgp allele + frequency matching)
rm -f *meta*
FILES=$(ls *.vcf.gz)
N=24
time (
for FILE in ${FILES}; do
	((i=i%N)); ((i++==0)) && wait
	(echo "Extracting meta data of $FILE"
	outfile=$(echo "${FILE}" | sed 's/vcf.gz/meta.vcf/g')
	zcat "${FILE}" | awk '{ print $1, $2, $3, $4, $5, $6, $7, $8, $9 }' > "${outfile}"
	gzip -f "${outfile}") &
done
wait
)