#!/bin/bash

# ===================================================================
# === ApoE-e4 association analysis (rs429358$ x rs7412 haplotype) ===
# ===================================================================

# get arguments
trait="$1" # trait="gap_gm"
targetDir="$2" # targetDir="results/${trait}/apoe"
geneticsDir="$3" # geneticsDir="data/genetics"
maf=$4 # maf=0.01
threads=$5 # threads=50

# echo settings
echo $'\n'"--- ApoE-e4 association settings ---"
echo "trait: "${trait}
echo "targetDir: "${targetDir}
echo "geneticsDir: "${geneticsDir}
echo "maf: "${maf}
echo "threads: "${threads}$'\n'

# set targetdir and make folder
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# write apoe snps to vcf (requires phased data)
plink2 --pfile ${geneticsDir}/chr19/imp_mri_qc/chr19_mri_qc \
--extract <(echo rs429358$'\n'rs7412) \
--export vcf vcf-dosage=HDS \
--freq \
--out ${targetDir}/apoe

# ApoE-e4 scoring
echo "Score ApoE-e4 alleles."
scriptDir=$(dirname "$0")
Rscript ${scriptDir}/apoe.R "${targetDir}/apoe.vcf" "${targetDir}/apoe.e4.vcf"

# run association
echo "Run association analysis."
plink2 --vcf ${targetDir}/apoe.e4.vcf \
--double-id \
--glm no-x-sex hide-covar \
--maf 0.01 \
--pheno data/${trait}/${trait}.txt \
--covar data/${trait}/covs.txt \
--covar-variance-standardize \
--vif 100000 \
--threads ${threads} \
--freq \
--out ${targetDir}/apoe.e4

# how many covariates
ncovs=$(head -1 data/${trait}/covs.txt | awk '{print NF-2}')

# merge freq and gwas results 
header="CHR"$'\t'BP$'\t'"ID"$'\t'"A1"$'\t'"A2"$'\t'"A1_FREQ"$'\t'"BETA"$'\t'"SE"$'\t'"P"$'\t'"N"
awk -v header="$header" 'BEGIN { print header } \
NR==FNR { if ($4==$6) { chr[$3]=$1; pos[$3]=$2; snp[$3]=$3; a1[$3]=$6; a2[$3]=$5; b[$3]=$9; se[$3]=$10; p[$3]=$12; n[$3]=$8; next } else if ($5==$6) { chr[$3]=$1; pos[$3]=$2; snp[$3]=$3; a1[$3]=$6; a2[$3]=$4; b[$3]=$9; se[$3]=$10; p[$3]=$12; n[$3]=$8; next }} \
$2 in chr { if ($4==a1[$2]) { print chr[$2], pos[$2], snp[$2], a1[$2], a2[$2], $5, b[$2], se[$2], p[$2], n[$2] } else { print chr[$2], pos[$2], snp[$2], a1[$2], a2[$2], 1-$5, b[$2], se[$2], p[$2], n[$2] }}' OFS='\t' ${targetDir}/apoe.e4.*glm.linear ${targetDir}/apoe.e4.afreq > ${targetDir}/apoe.sumstats.txt
echo "Adding allele frequencies for apoe-e4 completed."
	
# add Z column
echo "Adding Z score column."
header="CHR"$'\t'BP$'\t'"ID"$'\t'"A1"$'\t'"A2"$'\t'"A1_FREQ"$'\t'"BETA"$'\t'"SE"$'\t'"Z"$'\t'"P"$'\t'"N"
awk -v header="$header" 'BEGIN { print header } FNR>1 { print $1, $2, $3, $4, $5, $6, $7, $8, $7/$8, $9, $10}' OFS='\t' ${targetDir}/apoe.sumstats.txt > ${targetDir}/apoe.sumstats.tmp.txt
mv -f ${targetDir}/apoe.sumstats.tmp.txt ${targetDir}/apoe.sumstats.txt

# add ETA2 column 
echo "Adding ETA2 column."
header="CHR"$'\t'BP$'\t'"ID"$'\t'"A1"$'\t'"A2"$'\t'"A1_FREQ"$'\t'"BETA"$'\t'"SE"$'\t'"Z"$'\t'"P"$'\t'"N"$'\t'"ETA2"
awk -v header="$header" -v ncovs="${ncovs}" 'BEGIN { print header } FNR>1 { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $9^2/($9^2+($11-2-ncovs)) }' OFS='\t' ${targetDir}/apoe.sumstats.txt > ${targetDir}/apoe.sumstats.tmp.txt
mv -f ${targetDir}/apoe.sumstats.tmp.txt ${targetDir}/apoe.sumstats.txt

# concatenate plink results
echo "Concatenating plink results."
header="CHROM"$'\t'BP$'\t'"ID"$'\t'"REF"$'\t'"ALT"$'\t'"A1"$'\t'"TEST"$'\t'"OBC_CT"$'\t'"BETA"$'\t'"SE"$'\t'"T_STAT"$'\t'"P"
awk -v header="$header" 'BEGIN { print header } FNR>1' ${targetDir}/*.glm.linear > ${targetDir}/apoe.sumstats_plink.txt

# clean up folder
echo "Cleaning up folder."
rm -rf ${targetDir}/*.afreq
rm -rf ${targetDir}/*.linear
rm -rf ${targetDir}/*.vcf
chmod 770 ${targetDir}/*
echo "--- ApoE-e4 association analysis completed. --- "

