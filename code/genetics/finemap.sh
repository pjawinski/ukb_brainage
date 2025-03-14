#!/bin/bash

# ====================================
# === Get 95% credible set of SNPs ===
# ====================================

# get arguments
trait="${1}" # trait=gap_wm
targetDir="${2}" # targetDir="results/${trait}/finemap"
bedFiles="${3}" # bedFiles='data/genetics/chr${i}/imp_mri_qc_EURjoined/bed/chr${i}_mri_qc'
bgenFiles="${4}" # bgenFiles='data/genetics/chr${i}/imp_mri_qc_EURjoined/bgen/chr${i}_mri_qc'
LDsample="${5}" # LDsample="data/${trait}/${trait}.txt"
conditionalFile="${6}" # conditionalFile="results/${trait}/conditional/conditional.cleaned.tophits.annovar.txt"
sumstats="${7}" # sumstats="results/${trait}/gwas/sumstats.txt.gz"
sumstatsCols="${8}" # sumstatsCols="CHR,BP,ID,A1,A2,A1_FREQ,BETA,SE,P,N"
pthresh="${9}" # pthresh=5E-8
nthreads="${10}" # nthreads=56

# echo settings
echo $'\n'"--- Fine-mapping Settings ---"
echo "trait: ${trait}"
echo "targetDir: ${targetDir}"
echo "bedFiles: ${bedFiles}"
echo "bgenFiles: ${bgenFiles}"
echo "LDsample: ${LDsample}"
echo "conditionalFile: ${conditionalFile}"
echo "sumstats: ${sumstats}"
echo "sumstatsCols: ${sumstatsCols}"
echo "pthresh: ${pthresh}"
echo "nthreads: ${nthreads}"$'\n'
sleep 10

# set targetdir and make folder
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# Get relevant columns from sumstats file
echo "Getting CHR, BP, ID, A1, A2, A1_FREQ, BETA, SE, P, and N from sumstats File."
awk -v cols="${sumstatsCols}" '
  BEGIN { ncols=split(cols,colnames,","); print "CHR\tBP\tID\tA1\tA2\tA1_FREQ\tBETA\tSE\tP\tN" } 
  NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } } }
  NR > 1 { print $colidx[1], $colidx[2], $colidx[3], $colidx[4], $colidx[5], $colidx[6], $colidx[7], $colidx[8], $colidx[9], $colidx[10]}' OFS="\t" <(gzip -dc "${sumstats}") > "${targetDir}"/finemap.sumstats.txt

# get list of index snps
echo "Getting index snps and respective chromosomes."
snps=$(awk -F'\t' -v pthresh="${pthresh}" -v cols="ID,LEAD_SNP,CHR,BP,P,LEAD_SNP_pJ" 'BEGIN { ncols=split(cols,colnames,",") } NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } }; next } $colidx[1] == $colidx[2] && $colidx[5] < pthresh && $colidx[6] < pthresh { print $colidx[1] }' "${conditionalFile}"); snps=($snps)
chr=$(awk -F'\t' -v pthresh="${pthresh}" -v cols="ID,LEAD_SNP,CHR,BP,P,LEAD_SNP_pJ" 'BEGIN { ncols=split(cols,colnames,",") } NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } }; next } $colidx[1] == $colidx[2] && $colidx[5] < pthresh && $colidx[6] < pthresh { print $colidx[3] }' "${conditionalFile}"); chr=($chr)
pos=$(awk -F'\t' -v pthresh="${pthresh}" -v cols="ID,LEAD_SNP,CHR,BP,P,LEAD_SNP_pJ" 'BEGIN { ncols=split(cols,colnames,",") } NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } }; next } $colidx[1] == $colidx[2] && $colidx[5] < pthresh && $colidx[6] < pthresh { print $colidx[4] }' "${conditionalFile}"); pos=($pos)

# get participant IDs (serving as reference sample for clumping)
echo "Getting participant IDs of LD sample."
awk '{ print $1, $2 }' "${LDsample}" > "${targetDir}"/finemap.subs.txt

# loop across loci
echo "Starting to finemap all loci."
(for (( i=0; i<${#snps[@]}; i++ )); do (

	# create folder, get all snps in 10 Mbp flanking region, and determine LD with index snp
	# note: this can now also be done by applying plink2 to pgen dosage data (use --r2-unphased for agreement with ldstore)
	mkdir -p "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"
	awk -F'\t' -v chr="${chr[$i]}" -v pos="${pos[$i]}" 'NR==1 { next } $1==chr && ($2 > pos - 10000000) && ($2 < pos + 10000000) { print $3 }' "${targetDir}"/finemap.sumstats.txt > "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/flankingSNPs.txt
	plink \
	--bfile $(eval echo "$(echo "${bedFiles}" | sed 's/${i}/${chr[$i]}/g')")  \
	--keep "${targetDir}"/finemap.subs.txt \
	--extract "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/flankingSNPs.txt \
	--r2 \
	--ld-snp "${snps[$i]}" \
	--ld-window-kb 10000 \
	--ld-window 999999 \
	--ld-window-r2 0.0 \
	--out "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/flankingSNPs
	awk '{print $1,$2,$3,$4,$5,$6,$7}' "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/flankingSNPs.ld > "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/flankingSNPs.ld.tmp
	\mv "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/flankingSNPs.ld.tmp "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/flankingSNPs.ld

	# get furthest observation with LD > 0.1
	minBP=$(awk -v pos="${pos[$i]}" 'BEGIN {BP=pos} NR==1 {next} $5<BP && $7 >= 0.1 {BP=$5;next} END {print BP}' "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/flankingSNPs.ld)
	maxBP=$(awk -v pos="${pos[$i]}" 'BEGIN {BP=pos} NR==1 {next} $5>BP && $7 >= 0.1 {BP=$5;next} END {print BP}' "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/flankingSNPs.ld)

	# make sure that window is at least 500 kbp
	if [ $((${pos[$i]} - minBP)) -lt 250000 ]; then minBP=$((${pos[$i]} - 250000)); fi
	if [ $((maxBP - ${pos[$i]})) -lt 250000 ]; then maxBP=$((${pos[$i]} + 250000)); fi

	# create .z file including statistics of SNPs in range of minBP and maxBP
	# info: allele1 and allele2 need to be in line with .bgen file
	# info: after PLINK conversion from pgen to bgen, ALT appears to be allele1 and REF appears to be allele2
	# info: make sure that allele2 is effect allele (in line with output of SNPTEST named by FINEMAP)
	header="rsid chromosome position allele1 allele2 maf beta se"
	awk -F'\t' -v header="${header}" -v chr="${chr[$i]}" -v minBP="${minBP}" -v maxBP="${maxBP}" '
		BEGIN { print header }
		NR==FNR { snp[$2]=$5"_"$6; next}
		FNR==1 { next }
		$6<=0.5 { maf=$6 }
		$6>0.5 { maf=1-$6 }
		$1==chr && ($2 >= minBP) && ($2 <= maxBP) && snp[$3]==$4"_"$5 { print $3, $1, $2, $4, $5, maf, -$7, $8}
		$1==chr && ($2 >= minBP) && ($2 <= maxBP) && snp[$3]==$5"_"$4 { print $3, $1, $2, $5, $4, maf, $7, $8 }' \
		"$(eval echo "$(echo "${bedFiles}" | sed 's/${i}/${chr[$i]}/g')").bim" "${targetDir}"/finemap.sumstats.txt > "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/ldstore.z

	# create .ld file with ldstore
	n=$(tail -n +3 $(eval echo "$(echo "${bgenFiles}" | sed 's/${i}/${chr[$i]}/g')".sample) | wc -l)
	echo "z;bgen;bgi;bcor;ld;n_samples" > "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/ldstore.master
	echo "${targetDir}/${chr[$i]}_${pos[$i]}/ldstore.z;$(eval echo "$(echo "${bgenFiles}" | sed 's/${i}/${chr[$i]}/g')".bgen);$(eval echo "$(echo "${bgenFiles}" | sed 's/${i}/${chr[$i]}/g')".bgen.bgi);${targetDir}/${chr[$i]}_${pos[$i]}/ldstore.bcor;${targetDir}/${chr[$i]}_${pos[$i]}/ldstore.ld;$n" >> "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/ldstore.master
	ldstore \
	--in-files "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/ldstore.master \
	--read-only-bgen \
	--n-threads ${nthreads} \
	--write-bcor

	# run finemapping
	n=$(awk -F'\t' -v chr="${chr[$i]}" 'BEGIN { max=0 } $1==chr && $10>max { max=$10 } END { print max }' "${targetDir}"/finemap.sumstats.txt)
	echo "z;bcor;snp;config;cred;log;k;n_samples" > "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/finemap.master
	echo "${targetDir}/${chr[$i]}_${pos[$i]}/ldstore.z;${targetDir}/${chr[$i]}_${pos[$i]}/ldstore.bcor;${targetDir}/${chr[$i]}_${pos[$i]}/finemap.snp;${targetDir}/${chr[$i]}_${pos[$i]}/finemap.config;${targetDir}/${chr[$i]}_${pos[$i]}/finemap.cred;${targetDir}/${chr[$i]}_${pos[$i]}/finemap.log;${targetDir}/${chr[$i]}_${pos[$i]}/finemap.k;$n" >> "${targetDir}/${chr[$i]}_${pos[$i]}/finemap.master"
	finemap \
	--sss \
	--n-causal-snps 10 \
	--prob-cred-set 0.95 \
	--in-files "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/finemap.master \
	--force-n-samples 1 \
	--n-threads ${nthreads} \
	--dataset 1 \
	--log "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/finemap.log

	# get list of snps with a cumulative posterior inclusion probability of 95%
	snpTotal=$(grep "Number of SNPs" "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/finemap.log_sss | awk -F" : " '{ print $2}')
	snpCausal=$(grep "Post-expected # of causal SNPs" "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/finemap.log_sss | awk -F" : " '{ print $2}')
	regionalh2=$(grep "Regional SNP heritability" "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/finemap.log_sss | awk '{ print $6,$12}' | sed 's/)//g')
	bestK=$(grep -A 11 "Post-Pr(# of causal SNPs is k)" "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/finemap.log_sss | tail -10 | awk 'BEGIN {maxProb=$3;k=$1} $3>maxProb {maxProb=$3;k=$1} END {print k}')

	# output credible set
	scriptDir=$(dirname "$0") # scriptDir="code/genetics" 
	Rscript "${scriptDir}"/finemap.R \
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
awk 'NR==1 {print;next} FNR==1 {next} {print}' "${targetDir}"/finemap.*DF.txt > "${targetDir}"/finemap.df.txt
awk -F'\t' 'NR==1 {print; next} $2=="X" {$2=23} $2=="Y" {$2=24} $2=="XY" {$2=25} $2=="MT" {$2=26} {print}' OFS='\t' "${targetDir}"/finemap.df.txt > "${targetDir}"/finemap.df.tmp
(head -n 1 "${targetDir}"/finemap.df.tmp && tail -n +2 "${targetDir}"/finemap.df.tmp | sort -t$'\t' -k2,2g -k3,3g -k12,12g -k22,22g) \
	| awk -F'\t' 'NR==1 {print; next} $2==23 {$2="X"} $2==24 {$2="Y"} $2==25 {$2="XY"} $2==26 {$2="MT"} {print}' OFS='\t' \
	> "${targetDir}"/finemap.df.txt; \rm -f "${targetDir}"/finemap.df.tmp 

awk 'NR==1 {print;next} FNR==1 {next} {print}' "${targetDir}"/finemap.*LS.txt > "${targetDir}"/finemap.ls.txt
awk -F'\t' 'NR==1 {print; next} $2=="X" {$2=23} $2=="Y" {$2=24} $2=="XY" {$2=25} $2=="MT" {$2=26} {print}' "${targetDir}"/finemap.ls.txt > "${targetDir}"/finemap.ls.tmp
(head -n 1 "${targetDir}"/finemap.ls.tmp && tail -n +2 "${targetDir}"/finemap.ls.tmp | sort -t$'\t' -k2,2g -k3,3g) \
	| awk -F'\t' 'NR==1 {print; next} $2==23 {$2="X"} $2==24 {$2="Y"} $2==25 {$2="XY"} $2==26 {$2="MT"} {print}' OFS='\t' \
	> "${targetDir}"/finemap.ls.txt; \rm -f "${targetDir}"/finemap.ls.tmp 

# clean up
echo "Cleaning up."
rm -f "${targetDir}"/*finemapDF*
rm -f "${targetDir}"/*finemapLS*
rm -f "${targetDir}"/finemap.sumstats.txt
rm -f "${targetDir}"/finemap.subs.txt
mkdir -p "${targetDir}"/finemap.loci
for (( i=0; i<${#snps[@]}; i++ )); do
	mv "${targetDir}"/"${chr[$i]}"_"${pos[$i]}" "${targetDir}"/finemap.loci/
done
tar --use-compress-program="pigz --best --recursive" -cf "${targetDir}"/finemap.loci.tar.gz --directory="${targetDir}" finemap.loci --remove-files
chmod -R 770 "${targetDir}"/*
echo "--- Competed: Fine-Mapping --- "
