#!/bin/bash

# ====================================
# === Get 95% credible set of SNPs ===
# ====================================

# get arguments
trait="$1" # trait=gap_wm
targetDir="$2" # targetDir="results/${trait}/finemap"
geneticsDir="$3" # geneticsDir="data/genetics"
LDsample="$4" # LDsample="data/${trait}/${trait}.txt"
conditionalFile="$5" # conditionalFile="results/${trait}/conditional/conditional.cleaned.tophits.annovar.txt"
sumstats="$6"  # sumstats="results/${trait}/gwas/sumstats.txt.gz"
sumstatsPLINK="$7" # sumstatsPLINK="results/${trait}/gwas/sumstats_plink.txt.gz"
pthresh="$8" # pthresh=5E-8

# echo settings
echo $'\n'"--- Fine-mapping Settings ---"
echo "trait: "${trait}
echo "targetDir: "${targetDir}
echo "geneticsDir: "${geneticsDir}
echo "LDsample: "${LDsample}
echo "conditionalFile: "${conditionalFile}
echo "sumstats: "${sumstats}
echo "sumstatsPLINK: "${sumstatsPLINK}
echo "pthresh: "${pthresh}$'\n'

# set targetdir and make folder
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# get list of index snps
echo "Getting index snps and respective chromosomes."
snps=$(awk -F'\t' -v pthresh=${pthresh} '$2 == 1 && $17 < pthresh && $31 < pthresh { print $5 }' "${conditionalFile}"); snps=($snps)
chr=$(awk -F'\t' -v pthresh=${pthresh} '$2 == 1 && $17 < pthresh && $31 < pthresh { print $3 }' "${conditionalFile}"); chr=($chr)
pos=$(awk -F'\t' -v pthresh=${pthresh} '$2 == 1 && $17 < pthresh && $31 < pthresh { print $4 }' "${conditionalFile}"); pos=($pos)

# get participant IDs (serving as reference sample for clumping)
echo "Getting participant IDs of LD sample."
awk '{ print $1, $2 }' "${LDsample}" > "${targetDir}/finemap.subs.txt"

# loop across loci
echo "Starting to finemap all loci."
(
for (( i=0; i<${#snps[@]}; i++ )); do (

	# create folder, get all snps in 10 Mbp flanking region, and determine LD with index snp
	mkdir -p ${targetDir}/${chr[$i]}_${pos[$i]}
	awk -F'\t' -v chr="${chr[$i]}" -v pos="${pos[$i]}" 'NR==1 { next } $1==chr && ($2 > pos - 10000000) && ($2 < pos + 10000000) { print $3 }' <(gzip -dc ${sumstats}) > $targetDir/${chr[$i]}_${pos[$i]}/flankingSNPs.txt
	plink \
	--bfile "${geneticsDir}/chr${chr[$i]}/imp_mri_qc/bed/chr${chr[$i]}_mri_qc" \
	--keep "${targetDir}/finemap.subs.txt" \
	--extract "$targetDir/${chr[$i]}_${pos[$i]}/flankingSNPs.txt" \
	--r2 \
	--ld-snp "${snps[$i]}" \
	--ld-window-kb 10000 \
	--ld-window 999999 \
	--ld-window-r2 0 \
	--out ${targetDir}/${chr[$i]}_${pos[$i]}/flankingSNPs
	awk '{print $1,$2,$3,$4,$5,$6,$7}' $targetDir/${chr[$i]}_${pos[$i]}/flankingSNPs.ld > $targetDir/${chr[$i]}_${pos[$i]}/flankingSNPs.ld.tmp; \mv $targetDir/${chr[$i]}_${pos[$i]}/flankingSNPs.ld.tmp $targetDir/${chr[$i]}_${pos[$i]}/flankingSNPs.ld

	# get furthest observation with LD > 0.1
	minBP=$(awk -v pos="${pos[$i]}" 'BEGIN {BP=pos} NR==1 {next} $5<BP && $7 >= 0.1 {BP=$5;next} END {print BP}' $targetDir/${chr[$i]}_${pos[$i]}/flankingSNPs.ld)
	maxBP=$(awk -v pos="${pos[$i]}" 'BEGIN {BP=pos} NR==1 {next} $5>BP && $7 >= 0.1 {BP=$5;next} END {print BP}' $targetDir/${chr[$i]}_${pos[$i]}/flankingSNPs.ld)

	# create .z file including statistics of SNPs in range of minBP and maxBP
	# info: allele1 and allele2 need to be in line with .bgen file
	# info: after PLINK conversion from pgen to bgen, ALT appears to be allele1 and REF appears to be allele2
	# info: make sure that allele2 is effect allele (in line with output of SNPTEST named by FINEMAP)
	header="rsid chromosome position allele1 allele2 beta se"
	awk -F'\t' -v header="${header}" -v chr="${chr[$i]}" -v minBP="${minBP}" -v maxBP="${maxBP}" '
		BEGIN { print header } NR==1 { next }
		$1==chr && ($2 >= minBP) && ($2 < maxBP) && $4==$6 { print $3, $1, $2, $5, $4, $9, $10}
		$1==chr && ($2 >= minBP) && ($2 < maxBP) && $5==$6 { print $3, $1, $2, $5, $4, -$9, $10 }' OFS='\t' \
		<(gzip -dc ${sumstatsPLINK}) > $targetDir/${chr[$i]}_${pos[$i]}/ldstore.z

	header="rsid chromosome position allele1 allele2 maf beta se"
	awk -F'\t' -v header="${header}" 'BEGIN { print header }
		NR==FNR && $6<=0.5 {maf[$3]=$6; next}
		NR==FNR && $6 >0.5 {maf[$3]=1-$6; next}
		$1 in maf {print $1, $2, $3, $4, $5, maf[$1], $6, $7}' \
		<(gzip -dc ${sumstats}) $targetDir/${chr[$i]}_${pos[$i]}/ldstore.z > $targetDir/${chr[$i]}_${pos[$i]}/ldstore.z.tmp
	\mv $targetDir/${chr[$i]}_${pos[$i]}/ldstore.z.tmp $targetDir/${chr[$i]}_${pos[$i]}/ldstore.z

	# create .ld file with ldstore
	n=$(tail -n +3 ${geneticsDir}/chr${chr[$i]}/imp_mri_qc/bgen/chr${chr[$i]}_mri_qc.sample | wc -l)
	echo "z;bgen;bgi;bcor;ld;n_samples" > ${targetDir}/${chr[$i]}_${pos[$i]}/ldstore.master
	echo "${targetDir}/${chr[$i]}_${pos[$i]}/ldstore.z;${geneticsDir}/chr${chr[$i]}/imp_mri_qc/bgen/chr${chr[$i]}_mri_qc.bgen;${geneticsDir}/chr${chr[$i]}/imp_mri_qc/bgen/chr${chr[$i]}_mri_qc.bgen.bgi;${targetDir}/${chr[$i]}_${pos[$i]}/ldstore.bcor;${targetDir}/${chr[$i]}_${pos[$i]}/ldstore.ld;$n" >> $targetDir/${chr[$i]}_${pos[$i]}/ldstore.master
	ldstore \
	--in-files ${targetDir}/${chr[$i]}_${pos[$i]}/ldstore.master \
	--read-only-bgen \
	--write-bcor

	# run finemapping
	n=$(zcat ${sumstats} | awk -F'\t' -v chr="${chr[$i]}" '$1==chr { print $11; exit}')
	echo "z;bcor;snp;config;cred;log;k;n_samples" > ${targetDir}/${chr[$i]}_${pos[$i]}/finemap.master
	echo "${targetDir}/${chr[$i]}_${pos[$i]}/ldstore.z;${targetDir}/${chr[$i]}_${pos[$i]}/ldstore.bcor;${targetDir}/${chr[$i]}_${pos[$i]}/finemap.snp;${targetDir}/${chr[$i]}_${pos[$i]}/finemap.config;${targetDir}/${chr[$i]}_${pos[$i]}/finemap.cred;${targetDir}/${chr[$i]}_${pos[$i]}/finemap.log;${targetDir}/${chr[$i]}_${pos[$i]}/finemap.k;$n" >> $targetDir/${chr[$i]}_${pos[$i]}/finemap.master
	finemap \
	--sss \
	--n-causal-snps 10 \
	--prob-cred-set 0.95 \
	--in-files ${targetDir}/${chr[$i]}_${pos[$i]}/finemap.master \
	--dataset 1 \
	--log ${targetDir}/${chr[$i]}_${pos[$i]}/finemap.log

	# get list of snps with a cumulative posterior inclusion probability of 95%
	snpTotal=$(grep "Number of SNPs" ${targetDir}/${chr[$i]}_${pos[$i]}/finemap.log_sss | awk -F" : " '{ print $2}')
	snpCausal=$(grep "Post-expected # of causal SNPs" ${targetDir}/${chr[$i]}_${pos[$i]}/finemap.log_sss | awk -F" : " '{ print $2}')
	regionalh2=$(grep "Regional SNP heritability" ${targetDir}/${chr[$i]}_${pos[$i]}/finemap.log_sss | awk '{ print $6,$12}' | sed 's/)//g')
	bestK=$(grep -A 11 "Post-Pr(# of causal SNPs is k)" ${targetDir}/${chr[$i]}_${pos[$i]}/finemap.log_sss | tail -10 | awk 'BEGIN {maxProb=$3;k=$1} $3>maxProb {maxProb=$3;k=$1} END {print k}')

	# output credible set
	scriptDir=$(dirname "$0") # scriptDir="code/genetics/finemap" 
	Rscript "${scriptDir}/finemap.R" \
	"${targetDir}/${chr[$i]}_${pos[$i]}/finemap.cred${bestK}" \
	"${targetDir}/${chr[$i]}_${pos[$i]}/finemap.snp" \
	"${targetDir}/${chr[$i]}_${pos[$i]}/flankingSNPs.ld" \
	"${snpTotal}" \
	"${snpCausal}" \
	"${bestK}" \
	"${regionalh2}" \
	"${minBP}" \
	"${maxBP}" \
	"${targetDir}/${chr[$i]}_${pos[$i]}/finemapDF.txt" \
	"${targetDir}/${chr[$i]}_${pos[$i]}/finemapLS.txt"

	# move to targetDir
	\cp "${targetDir}/${chr[$i]}_${pos[$i]}/finemapDF.txt" "${targetDir}/finemap.${chr[$i]}_${pos[$i]}_finemapDF.txt"
	\cp "${targetDir}/${chr[$i]}_${pos[$i]}/finemapLS.txt" "${targetDir}/finemap.${chr[$i]}_${pos[$i]}_finemapLS.txt"
) &
done
wait
)

# summarise
awk 'NR==1 {print;next} FNR==1 {next} {print}' $targetDir/finemap.*DF.txt > $targetDir/finemap.df.txt
cat $targetDir/finemap.df.txt | awk -F'\t' 'NR==1 {print; next} $2=="X" {$2=23} $2=="Y" {$2=24} $2=="XY" {$2=25} $2=="MT" {$2=26} {print}' OFS='\t' > $targetDir/finemap.df.tmp
(head -n 1 $targetDir/finemap.df.tmp && tail -n +2 $targetDir/finemap.df.tmp | sort -t$'\t' -k2,2g -k3,3g -k10,10g -k20,20g) \
	| awk -F'\t' 'NR==1 {print; next} $2==23 {$2="X"} $2==24 {$2="Y"} $2==25 {$2="XY"} $2==26 {$2="MT"} {print}' OFS='\t' \
	> $targetDir/finemap.df.txt; \rm -f $targetDir/finemap.df.tmp 

awk 'NR==1 {print;next} FNR==1 {next} {print}' $targetDir/finemap.*LS.txt > $targetDir/finemap.ls.txt
cat $targetDir/finemap.ls.txt | awk -F'\t' 'NR==1 {print; next} $2=="X" {$2=23} $2=="Y" {$2=24} $2=="XY" {$2=25} $2=="MT" {$2=26} {print}' > $targetDir/finemap.ls.tmp
(head -n 1 $targetDir/finemap.ls.tmp && tail -n +2 $targetDir/finemap.ls.tmp | sort -t$'\t' -k2,2g -k3,3g) \
	| awk -F'\t' 'NR==1 {print; next} $2==23 {$2="X"} $2==24 {$2="Y"} $2==25 {$2="XY"} $2==26 {$2="MT"} {print}' OFS='\t' \
	> $targetDir/finemap.ls.txt; \rm -f $targetDir/finemap.ls.tmp 

# clean up
echo "Cleaning up."
rm -f ${targetDir}/*finemapDF*
rm -f ${targetDir}/*finemapLS*
rm -f ${targetDir}/finemap.subs.txt
mkdir -p "${targetDir}/finemap.loci"
for (( i=0; i<${#snps[@]}; i++ )); do
	mv "${targetDir}/${chr[$i]}_${pos[$i]}" "${targetDir}/finemap.loci/"
done
tar -cvzf "${targetDir}/finemap.loci.tar.gz" --directory=${targetDir} "finemap.loci" --remove-files
chmod -R 770 ${targetDir}/*
echo "--- Fine-Mapping finished. --- "
