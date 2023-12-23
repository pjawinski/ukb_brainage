#!/bin/bash

# ==========================
# === prepare glist_hg19 ===
# ==========================

# add gene biotype
mkdir -p /fast/software/glist_hg19
cd /fast/software/glist_hg19

# get gene list
wget -O glist-hg19-gcta.txt https://cnsgenomics.com/software/gcta/res/glist-hg19.txt 

# make non-numeric chromomosomes numeric (required for gcta software)
awk '{print $1}' glist-hg19-gcta.txt | sort -g | uniq -c 
awk '$1=="X" { print 23, $2, $3, $4; next} $1=="Y" { print 24, $2, $3, $4; next} $1=="XY" { print 25, $2, $3, $4; next} { print }' glist-hg19-gcta.txt > glist-hg19-gcta.tmp.txt
\mv glist-hg19-gcta.tmp.txt glist-hg19-gcta.txt

# create symbolic link
cd /slow/projects/ukb_brainage/data
ln -s /home/groups/markett/software/glist_hg19/ glist_hg19 

# make glist-hg19.txt tab-delimited
cd /home/groups/markett/ukb_brainage
awk '{ print $1, $2, $3, $4 }' OFS='\t' data/glist_hg19/glist-hg19-gcta.txt > data/glist_hg19/glist-hg19-gcta.tsv

# how many matches between glist-hg19-gcta.tsv and GRCh37_latest_genomic.edit.gff?
# - 22968 out of 26292
awk -F'\t' 'NR==FNR { gene[$4]=$0; next } $10 in gene { print gene[$10], $3, $9, $4, $5, $10, $11, $12}' OFS='\t' data/glist_hg19/glist-hg19-gcta.tsv <(gzip -dc data/refseq/GRCh37_latest_genomic.edit.gff.gz) | wc -l
cat data/glist_hg19/glist-hg19-gcta.tsv | wc -l

# reason for not finding gene? - gene synonym used in some cases
awk -F'\t' 'NR==FNR { gene[$10]; next } !($4 in gene) { print }' OFS='\t' <(gzip -dc data/refseq/GRCh37_latest_genomic.edit.gff.gz) data/glist_hg19/glist-hg19-gcta.tsv | head
grep -w ACPT <(gzip -dc data/refseq/GRCh37_latest_genomic.edit.gff.gz)
grep -w ADCK3 <(gzip -dc data/refseq/GRCh37_latest_genomic.edit.gff.gz)

# how many genes and pseudpgenes in glist-hg19.tsv?
# - 22171 genes
# - 797 peusdogenes
awk -F'\t' 'NR==FNR { gene[$4]=$0; next } $10 in gene { print $3 }' OFS='\t' data/glist_hg19/glist-hg19-gcta.tsv <(gzip -dc data/refseq/GRCh37_latest_genomic.edit.gff.gz) | sort | uniq -c

# what kind of biotypes in glist-hg19.tsv? - various biotypes
# (MAGMA only includes genes with biotype protein-coding - see /Users/philippe/Documents/magma_v1.07b_mac/NCBI37.3/README)
# 1 RNase_MRP_RNA
# 1 RNase_P_RNA
# 2 Y_RNA
# 12 antisense_RNA
# 24 guide_RNA
# 2109 lncRNA
# 1607 miRNA
# 1 misc_RNA
# 1 ncRNA_pseudogene
# 2 other
# 18009 protein_coding
# 37 pseudogene
# 36 snRNA
# 362 snoRNA
# 1 telomerase_RNA
# 759 transcribed_pseudogene
# 4 vault_RNA
awk -F'\t' 'NR==FNR { gene[$4]=$0; next } $10 in gene { print $11 }' OFS='\t' data/glist_hg19/glist-hg19-gcta.tsv <(gzip -dc data/refseq/GRCh37_latest_genomic.edit.gff.gz) | sort | uniq -c

# total number of genes and pseudogenes in data/refseq/GRCh37_latest_genomic.edit.gff.gz
# 28025 gene
# 15657 pseudogene
awk -F'\t' '{ print $3 }' <(gzip -dc data/refseq/GRCh37_latest_genomic.edit.gff.gz) | sort | uniq -c

# total number of protein-coding genes in data/refseq/GRCh37_latest_genomic.edit.gff.gz
# 19346
awk -F'\t' '$11=="protein_coding" { print $3 }' <(gzip -dc data/refseq/GRCh37_latest_genomic.edit.gff.gz) | wc -l

# show chromosomes in data/refseq/GRCh37_latest_genomic.edit.gff.gz
awk -F'\t' '{ print $9 }' <(gzip -dc data/refseq/GRCh37_latest_genomic.edit.gff.gz) | sort | uniq -c
awk -F'\t' '{ gsub(/[_].*/,"",$9) } { print $9 }' <(gzip -dc data/refseq/GRCh37_latest_genomic.edit.gff.gz) | sort -g | uniq -c

# extract with biotype 'protein_coding' on chr{{1..22,X,Y,M}} from data/refseq/GRCh37_latest_genomic.edit.gff.gz
chr_select=$(for i in {1..22} X Y XY M; do echo "$i"; done)
awk -F'\t' 'NR==FNR { chr[$1]=$1; next } { gsub(/[_].*/,"",$9) } ($9 in chr) && $3=="gene" && $11=="protein_coding" { print chr[$9], $4, $5, $10 }' <(echo "${chr_select}") <(gzip -dc data/refseq/GRCh37_latest_genomic.edit.gff.gz) > data/glist_hg19/glist-hg19-refseq.txt

# duplicate genes? - nope.
cat data/glist_hg19/glist-hg19-refseq.txt | wc -l
cat data/glist_hg19/glist-hg19-refseq.txt | awk '{print $4}' | sort -u | wc -l

# recode X, Y, and M
awk '{print $1}' data/glist_hg19/glist-hg19-refseq.txt | sort -g | uniq -c 
awk '$1=="X" { print 23, $2, $3, $4; next} $1=="Y" { print 24, $2, $3, $4; next} $1=="M" { print 26, $2, $3, $4; next} { print }' data/glist_hg19/glist-hg19-refseq.txt > data/glist_hg19/glist-hg19-refseq.tmp.txt
\mv data/glist_hg19/glist-hg19-refseq.tmp.txt data/glist_hg19/glist-hg19-refseq.txt


