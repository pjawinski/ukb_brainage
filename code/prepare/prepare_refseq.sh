#!/bin/bash

# =======================================================
# === prepare refseq file for gene biotype annotation ===
# =======================================================

# set working directory
cd /Users/philippe/Desktop/projects/ukb_brainage

# create refseq folder
mkdir data/refseq
cd data/refseq

# download refseq identifiers
# take care: latest README file differs from original one (105.20201022) | wget -O README.txt https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/README.txt
wget -O GRCh37_latest_genomic.gff.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/105.20201022/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz # original code: wget -O GRCh37_latest_genomic.gff.gz https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz

# get chromomsome translation
wget -O GRCh37_NCBI2UCSC.txt https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh37_NCBI2UCSC.txt

# remove ^M
cat -v GRCh37_NCBI2UCSC.txt | sed "s/\^M//g" > GRCh37_NCBI2UCSC.tmp.txt; \mv GRCh37_NCBI2UCSC.tmp.txt GRCh37_NCBI2UCSC.txt

# rename chr1..22 X Y MT
awk -F'\t' 'NR==FNR { gsub(/chr/,""); print $0 } $2=="M" { $2=="MT"}' OFS='\t' GRCh37_NCBI2UCSC.txt > GRCh37_NCBI2UCSC_chr122XYMT.txt

# get chromsome (translate), gene, biotype, description, synonyms from resources/GRCh37_latest_genomic.gff.gz
# 48,531 rows left
awk -F'\t' 'NR==FNR { chr[$1]=$2; next }
!($1 in chr) { next } $9!~/Name/ || $9!~/biotype/ { next }
{ gene=$9; sub(/.*Name=/,"",gene); sub(/[;].*/,"",gene)}
{ biotype=$9; sub(/.*gene_biotype=/,"",biotype); sub(/[;].*/,"",biotype) } $9!~/biotype/
{ description=$9; sub(/.*description=/,"",description); sub(/[;].*/,"",description) } $9!~/description/ { description="" }
{ synonym=$9; sub(/.*gene_synonym=/,"",synonym); sub(/[;].*/,"",synonym) } $9!~/synonym/ { synonym="" }
{ chrom=chr[$1] } !($1 in chr) { chrom="" } { print $1, $2, $3, $4, $5, $6, $7, $8, chrom, gene, biotype, description, synonym, $9 }' OFS='\t' GRCh37_NCBI2UCSC_chr122XYMT.txt <(gzip -dc GRCh37_latest_genomic.gff.gz) > GRCh37_latest_genomic.edit.gff

# number of genes and pseudogenes in resources/GRCh37_latest_genomic.edit.gff 
# 31427 17104
awk -F'\t' '$3=="gene" { count1++ } $3=="pseudogene" { count2++ } END { print count1, count2 }' GRCh37_latest_genomic.edit.gff

# sort by name, start, and stop coordinate and remove duplicate genes
# 43,682 rows left
cat GRCh37_latest_genomic.edit.gff | wc -l
cat GRCh37_latest_genomic.edit.gff | awk -F'\t' '{print $10}' | sort -u | wc -l
cat GRCh37_latest_genomic.edit.gff | sort -k10,10 -k1,1 -k4,4n -k5,5n > GRCh37_latest_genomic.edit.tmp.gff
\mv GRCh37_latest_genomic.edit.tmp.gff GRCh37_latest_genomic.edit.gff
awk -F'\t' '!($10 in gene) { print; gene[$10]; next}' OFS='\t' GRCh37_latest_genomic.edit.gff > GRCh37_latest_genomic.edit.tmp.gff
\mv GRCh37_latest_genomic.edit.tmp.gff GRCh37_latest_genomic.edit.gff

# synonyms
awk -F'\t' '{ print $13 }' GRCh37_latest_genomic.edit.gff > GRCh37_latest_genomic.edit.synonyms.gff

# gzip
chmod 770 *
gzip -f GRCh37_latest_genomic.edit.gff


