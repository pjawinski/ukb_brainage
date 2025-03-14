#!/bin/bash

# ========================================
# === Get 95% credible set of variants ===
# ========================================

# get arguments
trait="${1}" # trait=gap_wm
targetDir="${2}" # targetDir="results/${trait}/susieR"
bedFiles="${3}" # bedFiles='data/genetics/chr${i}/imp_mri_qc_EURjoined/bed/chr${i}_mri_qc'
bgenFiles="${4}" # bgenFiles='data/genetics/chr${i}/imp_mri_qc_EURjoined/bgen/chr${i}_mri_qc'
LDsample="${5}" # LDsample="data/genetics/chr1/imp_mri_qc_EURjoined/chr1_mri_qc.psam"
conditionalFile="${6}" # conditionalFile="results/${trait}/gwama/eur/conditional/conditional.cleaned.tophits.annovar.txt"
sumstats="${7}" # sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz"
sumstatsCols="${8}" # sumstatsCols="CHR,BP,ID,A1,A2,A1_FREQ,BETA,SE,P,N" # requires the following columns: "CHR,POS,SNPID,A1,A2,A1_FREQ,BETA,SE,P,N" 
pthresh="${9}" # pthresh=5E-8
nthreads="${10}" # nthreads=56

# echo settings
echo $'\n'"--- Get 95% credible set of variants | Settings ---"
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
  NR > 1 { print $colidx[1], $colidx[2], $colidx[3], $colidx[4], $colidx[5], $colidx[6], $colidx[7], $colidx[8], $colidx[9], $colidx[10]}' OFS="\t" <(gzip -dc "${sumstats}") > "${targetDir}"/susieR.sumstats.txt

# get list of index snps
echo "Getting index snps and respective chromosomes."
snps=$(awk -F'\t' -v pthresh="${pthresh}" -v cols="ID,LEAD_SNP,CHR,BP,P,LEAD_SNP_pJ" 'BEGIN { ncols=split(cols,colnames,",") } NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } }; next } $colidx[1] == $colidx[2] && $colidx[5] < pthresh && $colidx[6] < pthresh { print $colidx[1] }' "${conditionalFile}"); snps=($snps)
chr=$(awk -F'\t' -v pthresh="${pthresh}" -v cols="ID,LEAD_SNP,CHR,BP,P,LEAD_SNP_pJ" 'BEGIN { ncols=split(cols,colnames,",") } NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } }; next } $colidx[1] == $colidx[2] && $colidx[5] < pthresh && $colidx[6] < pthresh { print $colidx[3] }' "${conditionalFile}"); chr=($chr)
pos=$(awk -F'\t' -v pthresh="${pthresh}" -v cols="ID,LEAD_SNP,CHR,BP,P,LEAD_SNP_pJ" 'BEGIN { ncols=split(cols,colnames,",") } NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } }; next } $colidx[1] == $colidx[2] && $colidx[5] < pthresh && $colidx[6] < pthresh { print $colidx[4] }' "${conditionalFile}"); pos=($pos)

# get participant IDs (serving as reference sample for clumping)
echo "Getting participant IDs of LD sample."
awk '{ print $1, $2 }' "${LDsample}" > "${targetDir}"/susieR.subs.txt

# loop across loci
echo "Starting to finemap all loci."
(for (( i=0; i<${#snps[@]}; i++ )); do (

	# create folder, get all snps in 10 Mbp flanking region, and determine LD with index snp
	# note: this can now also be done by applying plink2 to pgen dosage data (use --r2-unphased for agreement with ldstore)
	mkdir -p "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"
	awk -F'\t' -v chr="${chr[$i]}" -v pos="${pos[$i]}" 'NR==1 { next } $1==chr && ($2 > pos - 10000000) && ($2 < pos + 10000000) { print $3 }' "${targetDir}"/susieR.sumstats.txt > "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/flankingSNPs.txt
	plink \
	--bfile $(eval echo "$(echo "${bedFiles}" | sed 's/${i}/${chr[$i]}/g')")  \
	--keep "${targetDir}"/susieR.subs.txt \
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
		"$(eval echo "$(echo "${bedFiles}" | sed 's/${i}/${chr[$i]}/g')").bim" "${targetDir}"/susieR.sumstats.txt > "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/ldstore.z

	# create .ld file with ldstore
	n=$(tail -n +3 $(eval echo "$(echo "${bgenFiles}" | sed 's/${i}/${chr[$i]}/g')".sample) | wc -l)
	echo "z;bgen;bgi;bcor;ld;n_samples" > "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/ldstore.master
	echo "${targetDir}/${chr[$i]}_${pos[$i]}/ldstore.z;$(eval echo "$(echo "${bgenFiles}" | sed 's/${i}/${chr[$i]}/g')".bgen);$(eval echo "$(echo "${bgenFiles}" | sed 's/${i}/${chr[$i]}/g')".bgen.bgi);${targetDir}/${chr[$i]}_${pos[$i]}/ldstore.bcor;${targetDir}/${chr[$i]}_${pos[$i]}/ldstore.ld;$n" >> "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/ldstore.master
	ldstore \
	--in-files "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/ldstore.master \
	--read-only-bgen \
	--n-threads ${nthreads} \
	--write-text

	# output credible set
	n=$(awk -v chr="${chr[$i]}" 'BEGIN {n=-1} NR==1 {next} $1==chr && $10>n {n=$10;next} END {print n}' "${targetDir}"/susieR.sumstats.txt)
	scriptDir=$(dirname "$0") # scriptDir="code/genetics"
	Rscript "${scriptDir}"/susieR.R \
		"${targetDir}/${chr[$i]}_${pos[$i]}/ldstore.ld" \
		"${targetDir}/${chr[$i]}_${pos[$i]}/ldstore.z" \
		"${snps[i]}" \
		"${n}" \
		"${targetDir}/susieR.${chr[$i]}_${pos[$i]}.DF.txt" \
		"${targetDir}/susieR.${chr[$i]}_${pos[$i]}.LS.txt"

) &
done
wait)

# pause for a moment
sleep 10

# summarise
awk -F'\t' 'NR==1 {print;next} FNR==1 {next} {print}' "${targetDir}"/susieR.*DF.txt > "${targetDir}"/susieR.df.txt
awk -F'\t' 'NR==1 {print;next} FNR==1 {next} {print}' "${targetDir}"/susieR.*LS.txt > "${targetDir}"/susieR.ls.txt
awk -F'\t' 'NR==1 {print;next} FNR==1 {next} $11>=0.5 {print}' "${targetDir}"/susieR.*DF.txt > "${targetDir}"/susieR.df.purity50.txt
awk -F'\t' 'NR==1 {print;next} FNR==1 {next}  $8>=0.5 {print}' "${targetDir}"/susieR.*LS.txt > "${targetDir}"/susieR.ls.purity50.txt

# clean up
echo "Cleaning up."
rm -f "${targetDir}"/*DF*
rm -f "${targetDir}"/*LS*
rm -f "${targetDir}"/susieR.sumstats.txt
rm -f "${targetDir}"/susieR.subs.txt
rm -f "${targetDir}"/*/susieR.subs.txt
mkdir -p "${targetDir}"/susieR.loci
for (( i=0; i<${#snps[@]}; i++ )); do
	pigz -f "${targetDir}"/"${chr[$i]}"_"${pos[$i]}"/ldstore.ld
	mv "${targetDir}"/"${chr[$i]}"_"${pos[$i]}" "${targetDir}"/susieR.loci/
done
tar --use-compress-program="pigz --best --recursive" -cf "${targetDir}"/susieR.loci.tar.gz --directory="${targetDir}" susieR.loci --remove-files
chmod -R 770 "${targetDir}"/*
echo "--- Competed: Get 95% credible set of variants --- "
