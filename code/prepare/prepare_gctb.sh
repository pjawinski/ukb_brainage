#!/bin/bash

# =================================
# === download and install gctb ===
# =================================

# set directory
mkdir -p /fast/software/gctb
cd /fast/software/gctb

# download gctb
wget https://cnsgenomics.com/software/gctb/download/gctb_2.5.2_Linux.zip
unzip gctb_2.5.2_Linux.zip
ln -s /fast/software/gctb/gctb_2.5.2_Linux/gctb /fast/software/bin/gctb

# download resources
mkdir -p /fast/software/gctb/resources
cd /fast/software/gctb/resources
wget https://sbayes.pctgplots.cloud.edu.au/data/SBayesRC/resources/v2.0/LD/Imputed/ukbEUR_Imputed.zip
wget https://sbayes.pctgplots.cloud.edu.au/data/SBayesRC/resources/v2.0/LD/HapMap3/ukbEUR_HM3.zip
wget https://sbayes.pctgplots.cloud.edu.au/data/SBayesRC/resources/v2.0/Annotation/annot_baseline2.2.zip
unzip ukbEUR_Imputed.zip
unzip https://sbayes.pctgplots.cloud.edu.au/data/SBayesRC/resources/v2.0/LD/HapMap3/ukbEUR_HM3.zip
unzip annot_baseline2.2.zip 
wget https://zenodo.org/records/3350914/files/ukbEURu_hm3_sparse.zip
unzip ukbEURu_hm3_sparse.zip
ukbEURu_hm3_sparse_mldm_list.txt
> ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_shrunk_sparse.list
for i in {1..22}; do
  echo "/fast/software/gctb/resources/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_chr${i}_v3_50k.ldm.sparse" >> "/fast/software/gctb/resources/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_shrunk_sparse.list"
done

# get additional baseline2.2 annotations /fast/software/gctb/resources/
cd /fast/software/gctb/resources/
wget https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/baselineLF_v2.2.UKB.tar.gz
tar -xzvf baselineLF_v2.2.UKB.tar.gz --wildcards '*.annot.gz'
cd baselineLF_v2.2.UKB
chmod 770 *.annot.gz
pigz -dc baselineLF2.2.UKB.{1..22}.annot.gz |
  awk 'NR==1 { output=$3"\tIntercept"; for (i=5; i<=NF; i++) { output=output"\t"$i }; print output; next }
      $1!="CHR" { output=$3"\t"1; for (i=5; i<=NF; i++) { output=output"\t"$i }; print output; next }' OFS='\t' > xf
pigz -f annot_baselineLF_v2.2.UKB
chmod 770 annot_baselineLF_v2.2.UKB.gz
rm -f *annot.gz

# create ld matrices in discovery sample
# ! caveat: next time, check first if all variants are in annotation dataaset (fortunately, they all were present in annot_baselineLF_v2.2.UKB.gz)
mkdir -p ukbEUR_32k
cd /fast/software/gctb/resources/ukbEUR_32k
wget https://cnsgenomics.com/software/gctb/download/ref4cM_v37.pos

(for i in {1..22}; do (
  mkdir "/slow/projects/ukb_brainage/data/genetics/chr${i}/imp_mri_qc/bed_01/"
  plink2 \
    --pfile "/slow/projects/ukb_brainage/data/genetics/chr${i}/imp_mri_qc/chr${i}_mri_qc" \
    --make-bed \
    --maf 0.01 \
    --out "/slow/projects/ukb_brainage/data/genetics/chr${i}/imp_mri_qc/bed_01/chr${i}_mri_qc"
  gctb \
    --bfile "/slow/projects/ukb_brainage/data/genetics/chr${i}/imp_mri_qc/bed_01/chr${i}_mri_qc" \
    --make-block-ldm \
    --block-info ref4cM_v37.pos \
    --out ldm \
    --thread 4
    ) &
done
wait)

# merge ldm.info
gctb --ldm ldm --merge-block-ldm-info --out ldm

# get annotations for snp set
awk 'NR==1 { next } NR==FNR {snp[$2]; next } FNR==1 { print; next} $1 in snp { print }' ldm/snp.info <(pigz -dc /fast/software/gctb/resources/baselineLF_v2.2.UKB/annot_baselineLF_v2.2.UKB.gz) > ldm/annot_baselineLF_v2.2.UKB_32k
pigz ldm/annot_baselineLF_v2.2.UKB_32k

# are all autosomal GWAS variants in ld dataset? - 3 variants missing outside of block boundaries
pigz -dc /slow/projects/ukb_brainage/results/gap_gm/discovery/gwas/sumstats.txt.gz | awk '$1!="X" && $1!="Y" && $1!="XY" && $1!="MT" && $1!="CHR" { print $1 }' | sort | uniq -c
pigz -dc /slow/projects/ukb_brainage/results/gap_gm/discovery/gwas/sumstats.txt.gz | awk '$1!="X" && $1!="Y" && $1!="XY" && $1!="MT" && $1!="CHR" { print $1 }' | wc -l # 9,373,077 variants
awk 'NR>1 { print }' /fast/software/gctb/resources/ukbEUR_32k/ldm/snp.info | wc -l # 9,373,075 variants
awk 'NR==1 { next } NR==FNR { id[$2]; next } FNR==1 { print; next } !($3 in id) && $1!="X" && $1!="Y" && $1!="XY" && $1!="MT" { print }
  ' /fast/software/gctb/resources/ukbEUR_32k/ldm/snp.info <(pigz -dc /slow/projects/ukb_brainage/results/gap_gm/discovery/gwas/sumstats.txt.gz) # 3 autosomal variants missing
awk 'NR==1 { print; next } $2==5 && $3<10000000 { print }' /fast/software/gctb/resources/ukbEUR_32k/ref4cM_v37.pos # variants outside of block boundaries (first block on chr5 starts at 11,940; variants have position 11,882 11,883 and 11,889)

# run eigen-decomposition in each block (MHC block takes 96 hours; next time check if option '--multiThreadEigen true' [spotted in c++ source files] allows multiple cores per block)
pid=$(echo $$) 
taskset -cp 0-55 $pid
> jobs.txt
for i in {1..591}; do
  echo "gctb --ldm ldm --make-ldm-eigen --block ${i} --out ldm" >> jobs.txt
done
parallel --progress --verbose -j 56 < jobs.txt

# get ld file for finemapping
cd /fast/software/gctb/resources/ukbEUR_32k
gctb \
--get-ld \
--ldm-eigen ldm \
--rsq 0.5 \
--out ldm/pairwise0.5 \
--thread 30
