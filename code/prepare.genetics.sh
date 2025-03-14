#!/bin/bash

# ==============================
# ====== Prepare genetics ======
# ==============================

# set working directory and create download folder
projectDir="/slow/projects/ukb_brainage"
mkdir -p ${projectDir}/data/genetics && cd ${projectDir}/data/genetics

# download utlities
mkdir -p util
wget -P util/ biobank.ndph.ox.ac.uk/showcase/util/gfetch
wget -P util/ biobank.ndph.ox.ac.uk/showcase/util/ukbgene
chmod 770 util/*

# manually put keyfile from basket into current folder and save filename in variable 'keyfile'
cp ${projectDir}/data/basket/20200217_2005558/data/*.key . # cp ${projectDir}/data/basket/20220914_2016290/data/*.key
keyfile=$(ls *.key)

# download meta files
mkdir -p meta ukb_snp_bim ukb_imp_mfi ukb_imp_bgi
./util/ukbgene rel -a.${keyfile}; mv *rel* meta/
wget -P meta/ biobank.ndph.ox.ac.uk/ukb/ukb/docs/ukb_genetic_data_description.txt
wget -P meta/ biobank.ctsu.ox.ac.uk/ukb/ukb/auxdata/ukb_snp_qc.txt
wget biobank.ctsu.ox.ac.uk/ukb/ukb/auxdata/ukb_snp_bim.tar; tar -xvf ukb_snp_bim/ukb_snp_bim.tar -C ukb_snp_bim/; rm -f ukb_snp_bim.tar
wget biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/ukb_imp_bgi.tgz; tar -zxvf ukb_imp_bgi.tgz -C ukb_imp_bgi/; rm -f ukb_imp_bgi.tgz 
wget biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/ukb_imp_mfi.tgz; tar -zxvf ukb_imp_mfi.tgz -C ukb_imp_mfi/; rm -f ukb_imp_mfi.tgz

# get genotype calls and imputations (5 downloads in parallel)
N=5; (for i in {1..22} X XY; do 
   ((j=j%N)); ((j++==0)) && wait (
      mkdir -p chr${i}/cal chr${i}/imp
      cd chr${i}/cal
         ../../util/ukbgene cal -c${i} -a../../${keyfile} # .bed file
         ../../util/ukbgene cal -c${i} -a../../${keyfile} -m; \mv *.fam ukb_cal_chr${i}_v2.fam # .fam file
         \cp ../../ukb_snp_bim/ukb_snp_chr${i}_v2.bim .; \mv *.bim ukb_cal_chr${i}_v2.bim
      cd ../imp
         ../../util/ukbgene imp -c${i} -a../../${keyfile} # .bgen file
         ../../util/ukbgene imp -c${i} -a../../${keyfile} -m # .sample file
      cd ../..
      ) & # get individual ids (.sample file)
done
wait)

(for i in Y MT; do (
   mkdir -p chr${i}/cal
   cd chr${i}/cal
      ../../util/ukbgene cal -c${i} -a../../${keyfile} # .bed file
      ../../util/ukbgene cal -c${i} -a../../${keyfile} -m; \mv *.fam ukb_cal_chr${i}_v2.fam # .fam file
      \cp ../../ukb_snp_bim/ukb_snp_chr${i}_v2.bim .; \mv *.bim ukb_cal_chr${i}_v2.bim
      cd -
   ) &
done
wait)

# create imp_mri .pgen files with selected MRI cases
subjects="${projectDir}/results/mri/iid.r2024.unrelated.txt"
N=2; (for i in {1..20} X XY; do 
   ((j=j%N)); ((j++==0)) && wait
      (sampleFile=$(ls chr${i}/imp/*.sample)
      mkdir -p chr${i}/imp_mri
      plink2 --bgen chr${i}/imp/ukb_imp_chr${i}_v3.bgen \
      --sample ${sampleFile} \
      --keep ${subjects} \
      --make-pgen \
      --out chr${i}/imp_mri/chr${i}_mri
      ) &
done
wait)

for i in Y MT; do
   mkdir -p chr${i}/imp_mri
   plink2 --bed chr${i}/cal/ukb_cal_chr${i}_v2.bed \
   --bim chr${i}/cal/ukb_cal_chr${i}_v2.bim \
   --fam $(ls chr${i}/cal/ukb*s4*.fam) \
   --keep ${subjects} \
   --make-pgen \
   --out chr${i}/imp_mri/chr${i}_mri
done

# get list of biallelic snps (by removing variants with identical positions or identical identifiers) and with info > 0.8
targetDir="qc_snplist"
mkdir -p $targetDir

(for i in {1..22} X XY; do (
awk 'NR==1 { next }
     NR==FNR && (seenPos[$1"_"$2]++) { dups[$1"_"$2]; next }
     NR==FNR && (seenId[$3]++) { dups[$3]; next }
     FNR==1 { print; next}
     NR!=FNR && !($1"_"$2 in dups) && !($3 in dups) { print }' chr${i}/imp_mri/chr${i}_mri.pvar chr${i}/imp_mri/chr${i}_mri.pvar \
     > ${targetDir}/chr${i}_nodups.txt
awk 'NR==FNR && $8 >= 0.8 { hq[$2]; next}
     NR==FNR { next } 
     FNR==1 { print  }
     $3 in hq { print }' ukb_imp_mfi/ukb_mfi_chr${i}_v3.txt ${targetDir}/chr${i}_nodups.txt  \
      > ${targetDir}/chr${i}_info_08_nodups.txt
awk 'NR==FNR && $8 >= 0.9 { hq[$2]; next}
     NR==FNR { next } 
     FNR==1 { print  }
     $3 in hq { print }' ukb_imp_mfi/ukb_mfi_chr${i}_v3.txt ${targetDir}/chr${i}_nodups.txt  \
      > ${targetDir}/chr${i}_info_09_nodups.txt
echo chr${i} completed.
) &
done
wait)

awk 'NR==1 { print } FNR==1 { next } { print }' $(for i in {1..22} X XY; do echo ${targetDir}/chr${i}_info_08_nodups.txt; done) > ${targetDir}/ukb_info_08_nodups.txt
awk 'NR==1 { next } { print $3 }' ${targetDir}/ukb_info_08_nodups.txt > ${targetDir}/ukb_info_08_nodups.snplist
awk 'NR==1 { print } FNR==1 { next } { print }' $(for i in {1..22} X XY; do echo ${targetDir}/chr${i}_info_09_nodups.txt; done) > ${targetDir}/ukb_info_09_nodups.txt
awk 'NR==1 { next } { print $3 }' ${targetDir}/ukb_info_09_nodups.txt > ${targetDir}/ukb_info_09_nodups.snplist
rm -f ${targetDir}/chr*

# ==========================================
# === create files for discovery dataset ===
# ==========================================

# create imp_mri_qc .pgen files
subjects="${projectDir}/results/mri/iid.discovery.txt"
snps="qc_snplist/ukb_info_08_nodups.snplist" 

N=2; (for i in {1..22} X XY; do 
   ((j=j%N)); ((j++==0)) && wait (
      mkdir -p chr${i}/imp_mri_qc
      plink2 \
      --pfile chr${i}/imp_mri/chr${i}_mri \
      --keep ${subjects} \
      --extract ${snps} \
      --maf 0.001 \
      --make-pgen \
      --out chr${i}/imp_mri_qc/chr${i}_mri_qc
   ) &
done
wait)

for i in Y MT; do
   mkdir -p chr${i}/imp_mri_qc
   plink2 \
   --pfile chr${i}/imp_mri/chr${i}_mri \
   --keep ${subjects} \
   --rm-dup exclude-all \
   --maf 0.001 \
   --geno 0.05 \
   --make-pgen \
   --out chr${i}/imp_mri_qc/chr${i}_mri_qc
done

# convert to .bed files
N=2; (for i in {1..22} X XY Y MT; do 
   ((j=j%N)); ((j++==0)) && wait (
      mkdir -p chr${i}/imp_mri_qc/bed
      plink2 \
      --pfile chr${i}/imp_mri_qc/chr${i}_mri_qc \
      --make-bed \
      --out chr${i}/imp_mri_qc/bed/chr${i}_mri_qc
   ) &
done
wait)

# for gcta software: create alternative .bim files with X, Y, XY, and MT converted to numeric values
for i in X Y XY MT; do
   \cp chr${i}/imp_mri_qc/bed/chr${i}_mri_qc.bim chr${i}/imp_mri_qc_${anc}/bed/chr${i}_mri_qc_original.bim
   awk '$1=="X" {$1=23; print}
        $1=="Y" {$1=24; print}
        $1=="XY" {$1=25; print}
        $1=="MT" {$1=26; print}' OFS='\t' chr${i}/imp_mri_qc/bed/chr${i}_mri_qc.bim > chr${i}/imp_mri_qc/bed/chr${i}_mri_qc_numeric.bim
done

# XY must be 23 for fastBAT gene-based analysis
awk '$1=="XY" {$1=23; print}' OFS='\t' chrXY/imp_mri_qc/bed/chrXY_mri_qc_original.bim > chrXY/imp_mri_qc/bed/chrXY_mri_qc_numeric_23.bim

# create bgen-1.2 for PRSice
snps="qc_snplist/ukb_info_09_nodups.snplist"
N=2; (for i in {1..22} X XY; do 
   ((j=j%N)); ((j++==0)) && wait
      (mkdir -p chr${i}/imp_mri_qc/bgen
      plink2 \
      --pfile chr${i}/imp_mri_qc/chr${i}_mri_qc \
      --export bgen-1.2 \
      --out chr${i}/imp_mri_qc/bgen/chr${i}_mri_qc
   ) &
done
wait)

# ================================================
# === create files for ukb replication dataset ===
# ================================================

# AFR AMR CSA EAS EUR MID NA
for anc in AFR AMR CSA EAS EUR MID NA; do
   count=$(cat $projectDir/results/mri/iid.replication.${anc}.txt | wc -l)
   count=$(expr $count - 1)
   echo "$anc $count"
done

# create imp_mri_qc .pgen files
snps="qc_snplist/ukb_info_08_nodups.snplist" 
for anc in AFR AMR CSA EAS EUR MID; do
   subjects="$projectDir/results/mri/iid.replication.${anc}.txt"
   N=2; (for i in {1..22} X XY; do 
      ((j=j%N)); ((j++==0)) && wait
         (mkdir -p chr${i}/imp_mri_qc_${anc}
         plink2 --pfile chr${i}/imp_mri/chr${i}_mri \
         --keep ${subjects} \
         --extract ${snps} \
         --maf 0.001 \
         --make-pgen \
         --out chr${i}/imp_mri_qc_${anc}/chr${i}_mri_qc
         ) &
   done
   wait)

   for i in Y MT; do
      mkdir -p chr${i}/imp_mri_qc_${anc}
      plink2 \
      --pfile chr${i}/imp_mri/chr${i}_mri \
      --keep ${subjects} \
      --rm-dup exclude-all \
      --maf 0.001 \
      --geno 0.05 \
      --make-pgen \
      --out chr${i}/imp_mri_qc_${anc}/chr${i}_mri_qc
   done

   # convert to .bed files
   N=2; (for i in {1..22} X XY Y MT; do
      ((j=j%N)); ((j++==0)) && wait
         (mkdir -p chr${i}/imp_mri_qc_${anc}/bed
         plink2 \
         --pfile chr${i}/imp_mri_qc_${anc}/chr${i}_mri_qc \
         --make-bed \
         --out chr${i}/imp_mri_qc_${anc}/bed/chr${i}_mri_qc
         ) &
   done
   wait)

   # create bgen-1.2 with maf 0.01 for PRSice2
   N=2; (for i in {1..22} X XY Y MT; do 
      ((j=j%N)); ((j++==0)) && wait
         (mkdir -p chr${i}/imp_mri_qc_${anc}/bgen
         plink2 \
         --pfile chr${i}/imp_mri_qc_${anc}/chr${i}_mri_qc \
         --export bgen-1.2 \
         --out chr${i}/imp_mri_qc_${anc}/bgen/chr${i}_mri_qc
         ) &
   done
   wait)

   # for gcta software: create alternative .bim files with X, Y, XY, and MT converted to numeric values
   for i in X Y XY MT; do
      \cp chr${i}/imp_mri_qc_${anc}/bed/chr${i}_mri_qc.bim chr${i}/imp_mri_qc_${anc}/bed/chr${i}_mri_qc_original.bim
      awk '$1=="X" {$1=23; print}
           $1=="Y" {$1=24; print}
           $1=="XY" {$1=25; print}
           $1=="MT" {$1=26; print}' OFS='\t' chr${i}/imp_mri_qc_${anc}/bed/chr${i}_mri_qc.bim > chr${i}/imp_mri_qc_${anc}/bed/chr${i}_mri_qc_numeric.bim
   done
done

# =====================================================================
# === create files for joined EUR discovery and replication dataset ===
# =====================================================================

# create imp_mri_qc .pgen files
snps="qc_snplist/ukb_info_08_nodups.snplist" 
subjects="$projectDir/results/mri/iid.EURjoined.txt"
N=2; (for i in {1..22} X XY; do 
   ((j=j%N)); ((j++==0)) && wait
      (mkdir -p chr${i}/imp_mri_qc_EURjoined
      plink2 --pfile chr${i}/imp_mri/chr${i}_mri \
      --keep ${subjects} \
      --extract ${snps} \
      --maf 0.001 \
      --make-pgen \
      --out chr${i}/imp_mri_qc_EURjoined/chr${i}_mri_qc
      ) &
done
wait)

for i in Y MT; do
   mkdir -p chr${i}/imp_mri_qc_EURjoined
   plink2 \
   --pfile chr${i}/imp_mri/chr${i}_mri \
   --keep ${subjects} \
   --rm-dup exclude-all \
   --maf 0.001 \
   --geno 0.05 \
   --make-pgen \
   --out chr${i}/imp_mri_qc_EURjoined/chr${i}_mri_qc
done

# convert to .bed files
N=2; (for i in {1..22} X XY Y MT; do
   ((j=j%N)); ((j++==0)) && wait
      (mkdir -p chr${i}/imp_mri_qc_EURjoined/bed
      plink2 \
      --pfile chr${i}/imp_mri_qc_EURjoined/chr${i}_mri_qc \
      --make-bed \
      --out chr${i}/imp_mri_qc_EURjoined/bed/chr${i}_mri_qc
      ) &
done
wait)

# create bgen-1.2 with maf 0.01 for PRSice2
N=2; (for i in {1..22} X XY Y MT; do 
   ((j=j%N)); ((j++==0)) && wait
      (mkdir -p chr${i}/imp_mri_qc_EURjoined/bgen
      plink2 \
      --pfile chr${i}/imp_mri_qc_EURjoined/chr${i}_mri_qc \
      --export bgen-1.2 \
      --out chr${i}/imp_mri_qc_EURjoined/bgen/chr${i}_mri_qc
      ) &
done
wait)

# get .bgi files
N=2; (for i in {1..22} X XY Y MT; do 
   ((j=j%N)); ((j++==0)) && wait
   bgenix -index -g chr${i}/imp_mri_qc_EURjoined/bgen/chr${i}_mri_qc.bgen &
done
wait)

# for gcta software: create alternative .bim files with X, Y, XY, and MT converted to numeric values
for i in X Y XY MT; do
   \cp chr${i}/imp_mri_qc_EURjoined/bed/chr${i}_mri_qc.bim chr${i}/imp_mri_qc_EURjoined/bed/chr${i}_mri_qc_original.bim
   awk '$1=="X" {$1=23; print}
        $1=="Y" {$1=24; print}
        $1=="XY" {$1=25; print}
        $1=="MT" {$1=26; print}' OFS='\t' chr${i}/imp_mri_qc_EURjoined/bed/chr${i}_mri_qc.bim > chr${i}/imp_mri_qc_EURjoined/bed/chr${i}_mri_qc_numeric.bim
done

# XY must be 23 for fastBAT gene-based analysis
awk '$1=="XY" {$1=23; print}' OFS='\t' chrXY/imp_mri_qc_EURjoined/bed/chrXY_mri_qc_original.bim > chrXY/imp_mri_qc_EURjoined/bed/chrXY_mri_qc_numeric_23.bim


# ===================================================================
# === calculate genetic principal components for discovery sample ===
# ===================================================================

# select same snps used for PCA in UK Biobank (Bycroft et al.)
# [make sure you're in data/genetics/ folder]
mkdir -p pca
awk '$118==1 { print $1 }' meta/ukb_snp_qc.txt > pca/pca.snplist

# extract PCA snps
pcaSNPs="pca/pca.snplist"
subjects="${projectDir}/results/mri/iid.discovery.txt"

for i in {1..22}; do
   mkdir -p pca/chr
   plink2 \
   --pfile chr${i}/imp_mri/chr${i}_mri \
   --keep $subjects \
   --extract $pcaSNPs \
   --export bgen-1.2 \
   --out pca/chr/chr${i}
done

# Create single file with pca snps
mkdir -p pca/bgen
cat-bgen -g $(for i in {1..22}; do echo pca/chr/chr${i}.bgen; done) -og pca/bgen/merged.bgen
cp pca/chr/chr1.sample pca/bgen/merged.sample
bgenix -g pca/bgen/merged.bgen -index

# convert bgen to pgen
mkdir -p pca/pgen/
plink2 \
--bgen pca/bgen/merged.bgen \
--make-pgen \
--out pca/pgen/merged

# perform pca
mkdir -p pca/results
plink2 \
--pfile pca/pgen/merged \
--seed 1585839423 \
--pca 20 approx var-wts \
--out pca/results/pca

# add principal components to vars file
mkdir -p ${projectDir}/results/genetics
header=$(head -1 ${projectDir}/results/mri/r2024.vars.txt)
header=$(echo "$header"$'\t'"$(head -1 pca/results/pca.eigenvec | awk '{ output = $2; for(i=3;i<=NF;i++) { output = output"\t"$i }; print output }')")
awk -v header="$header" '
   BEGIN { print header } 
   NR==1 { numPC=NF-1 }
   FNR==1 { next }
   NR==FNR { sub(/[_].*/,"",$1); id[$1]=$2; for(i=3;i<=NF;i++) { id[$1]=id[$1]"\t"$i }; next }
   $1 in id { print $0, id[$1]; next }
   { output=$0; for(i=1;i<=numPC;i++) { output=output"\tNA" }; print output }
   ' OFS='\t' pca/results/pca.eigenvec ${projectDir}/results/mri/r2024.vars.txt \
   > ${projectDir}/results/genetics/r2024.vars.pca.txt

# =========================================
# === create genetic relatedness matrix ===
# =========================================

# create folder [make sure you're in data/genetics/ folder]
mkdir -p grm/chr

# make chromosome-wise grms
for i in {1..22}; do
   gcta64 \
   --pfile chr${i}/imp_mri_qc/chr${i}_mri_qc \
   --chr ${i} \
   --maf 0.01 \
   --make-grm \
   --out grm/chr/chr${i} \
   --thread-num 112
done

# merge
ls grm/chr/*.grm.bin | sed -e 's/\..*$//' > grm/files4merge.txt
gcta64 \
--mgrm grm/files4merge.txt \
--make-grm \
--out grm/ukb_merged \
--thread-num 112
rm -f grm/files4merge.txt

