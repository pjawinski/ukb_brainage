#!/bin/bash

# ==============================
# === SMR eqtl/sqtl analysis ===
# ==============================

# get arguments
trait="$1" # trait="gap_gm"
targetDir="$2" # targetDir="results/${trait}/smr"
geneticsDir="$3" # geneticsDir="data/genetics"
sumstats="$4" # sumstats="results/${trait}/gwas/sumstats.txt.gz"
conditionalFile="$5" # conditionalFile="results/${trait}/conditional/conditional.cleaned.tophits.annovar.txt"
pthresh=$6 # pthresh=1E-6
LDsample="$7" # LDsample="data/${trait}/${trait}.txt"
threads=$8 # threads=100
clumpR2="$9" # clumpR2=0.8
smrCorrect="${10}" # smrCorrect="fdr"

# echo settings
echo $'\n'"--- Summary-data based mendelian randomization ---"
echo "trait: "${trait}
echo "targetDir: "${targetDir}
echo "geneticsDir: "${geneticsDir}
echo "sumstats: "${sumstats}
echo "conditionalFile: "${conditionalFile}
echo "pthresh: "${pthresh}
echo "LDsample: "${LDsample}
echo "threads: "${threads}
echo "clumpR2: "${clumpR2}
echo "smrCorrect: "${smrCorrect}$'\n'

# set targetdir and make folder
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# create input file from sumstats
echo "Creating input file."
header="SNP"$'\t'"A1"$'\t'"A2"$'\t'"freq"$'\t'"b"$'\t'"se"$'\t'"p"$'\t'"N"
awk -v header="$header" 'BEGIN { print header } NR > 1 { print $3, $4, $5, $6, $7, $8, $10, $11 }' OFS='\t' <(gzip -dc $sumstats) > ${targetDir}/smr.input.txt

# run smr analysis
echo "Running SMR analysis."
threadsPerAnalysis=$(expr ${threads} / 22)

# BrainMeta eQTL
mkdir -p ${targetDir}/smr.eqtl.output
N=22
(
for chr in {1..22}; do # X XY Y MT
	((i=i%N)); ((i++==0)) && wait 
	smr --bfile ${geneticsDir}/chr${chr}/imp_mri_qc/bed/chr${chr}_mri_qc \
	--gwas-summary ${targetDir}/smr.input.txt \
	--beqtl-summary data/smr/BrainMeta/BrainMeta_cis_eqtl_summary/BrainMeta_cis_eQTL_chr${chr} \
	--out ${targetDir}/smr.eqtl.output/chr${chr} \
	--thread-num ${threadsPerAnalysis} &
	done
	wait
)

# BrainMeta sQTL
mkdir -p ${targetDir}/smr.sqtl.output
N=22
(
for chr in {1..22}; do # X XY Y MT
	((i=i%N)); ((i++==0)) && wait 
	smr --bfile ${geneticsDir}/chr${chr}/imp_mri_qc/bed/chr${chr}_mri_qc \
	--gwas-summary ${targetDir}/smr.input.txt \
	--beqtl-summary data/smr/BrainMeta/BrainMeta_cis_sqtl_summary/BrainMeta_cis_sQTL_chr${chr} \
	--out ${targetDir}/smr.sqtl.output/chr${chr} \
	--thread-num ${threadsPerAnalysis} &
	done
	wait
)

# get summary of results
awk 'NR == 1 { print; next} FNR > 1 {sub(23,"X",$6); sub(24,"Y",$6); sub(25,"XY",$6); sub(26,"MT",$6); print $0}' OFS="\t" ${targetDir}/smr.eqtl.output/chr*.smr > ${targetDir}/smr.eqtl.txt
(head -n 1 ${targetDir}/smr.eqtl.txt && tail -n +2 ${targetDir}/smr.eqtl.txt | sort -k19,19g) > ${targetDir}/smr.eqtl.txt.tmp; \mv ${targetDir}/smr.eqtl.txt.tmp ${targetDir}/smr.eqtl.txt
awk 'NR == 1 { print; next} FNR > 1 {sub(23,"X",$6); sub(24,"Y",$6); sub(25,"XY",$6); sub(26,"MT",$6); print $0}' OFS="\t" ${targetDir}/smr.sqtl.output/chr*.smr > ${targetDir}/smr.sqtl.txt
(head -n 1 ${targetDir}/smr.sqtl.txt && tail -n +2 ${targetDir}/smr.sqtl.txt | sort -k19,19g) > ${targetDir}/smr.sqtl.txt.tmp; \mv ${targetDir}/smr.sqtl.txt.tmp ${targetDir}/smr.sqtl.txt

# filter results with p_HEIDI < 0.01 and (FDR < 0.05 OR Bonferroni < 0.05)
scriptDir=$(dirname "$0")
Rscript "${scriptDir}/smr.filter.R" "${targetDir}/smr.eqtl.txt" "${targetDir}/smr.eqtl.filtered.txt" "${smrCorrect}"
Rscript "${scriptDir}/smr.filter.R" "${targetDir}/smr.sqtl.txt" "${targetDir}/smr.sqtl.filtered.txt" "${smrCorrect}"

# =====================================================================
# === assign genes to independent snp-level discoveries by clumping ===
# =====================================================================

# get conditional discoveries
chr=$(awk -F'\t' -v pthresh="${pthresh}" 'NR==1 { next } $2==1 && $17<pthresh && $31<pthresh { print $3 }' ${conditionalFile}) 
chr=($chr)

snp=$(awk -F'\t' -v pthresh="${pthresh}" 'NR==1 { next } $2==1 && $17<pthresh && $31<pthresh { print $5 }' ${conditionalFile}) 
snp=($snp)

# get top snp ID and SMR P values from SMR results
cat ${targetDir}/smr.eqtl.filtered.txt | awk -F'\t' 'BEGIN { print "ID\tP" } NR==1 { next } !seen[$5]++ { print $5, $19}' OFS='\t' > "${targetDir}/smr.eqtl.sumstats4clumping.txt"
cat ${targetDir}/smr.sqtl.filtered.txt | awk -F'\t' 'BEGIN { print "ID\tP" } NR==1 { next } !seen[$5]++ { print $5, $19}' OFS='\t' > "${targetDir}/smr.sqtl.sumstats4clumping.txt"

# get participant IDs (they serve as reference sample for clumping)
awk '{ print $1, $2 }' "${LDsample}" > "${targetDir}/smr.subs.txt"

# do clumping
echo "Assign genes to independent snp-level discoveries."
mkdir -p "${targetDir}/smr.eqtl.clumping"
mkdir -p "${targetDir}/smr.sqtl.clumping"

task () {
for xqtl in eqtl sqtl; do
awk -v snp=${snp[i]} 'NR==1 { print; print snp, 1E-307; next } { print }' OFS='\t' "${targetDir}/smr.${xqtl}.sumstats4clumping.txt" | awk -F'\t' '!seen[$1]++ { print }' > "${targetDir}/smr.${xqtl}.clumping/sumstats4clumping_${snp[i]}.txt"
LD_src="${geneticsDir}/chr${chr[i]}/imp_mri_qc/bed/chr${chr[i]}_mri_qc"
plink --bfile "$LD_src" \
--keep "${targetDir}/smr.subs.txt" \
--clump "${targetDir}/smr.${xqtl}.clumping/sumstats4clumping_${snp[i]}.txt" \
--clump-snp-field "ID" \
--clump-field "P" \
--clump-p1 1E-300 \
--clump-p2 1.0 \
--clump-r2 ${clumpR2} \
--clump-kb 5000 \
--clump-verbose  \
--out "${targetDir}/smr.${xqtl}.clumping/${snp[i]}"
rm -rf "${targetDir}/smr.${xqtl}.clumping/sumstats4clumping_${snp[i]}.txt"
done
}

(
for i in $(seq 0 $((${#snp[@]}-1))); do
task "$i" &
done
wait
)

# Merge clumped files
cat ${targetDir}/smr.eqtl.clumping/*.clumped > ${targetDir}/smr.eqtl.clumping/clumped.tmp
cat ${targetDir}/smr.sqtl.clumping/*.clumped > ${targetDir}/smr.sqtl.clumping/clumped.tmp
awk '$11 ~ /^[0-9]+$/ { indexsnp = $3; print $3, $3, 0, 1, "A1/A2"; next } $1=="" || $1=="KB" || $1=="CHR" || $1=="(INDEX)" || $2 == "not" || $1 == "RANGE:" || $1 == "SPAN:" || $1 ~ /-/ { next } {print indexsnp, $1, $2, $3, $4}' OFS='\t' ${targetDir}/smr.eqtl.clumping/clumped.tmp > ${targetDir}/smr.eqtl.clumping/clumped; rm -f ${targetDir}/smr.eqtl.clumping/clumped.tmp
awk '$11 ~ /^[0-9]+$/ { indexsnp = $3; print $3, $3, 0, 1, "A1/A2"; next } $1=="" || $1=="KB" || $1=="CHR" || $1=="(INDEX)" || $2 == "not" || $1 == "RANGE:" || $1 == "SPAN:" || $1 ~ /-/ { next } {print indexsnp, $1, $2, $3, $4}' OFS='\t' ${targetDir}/smr.sqtl.clumping/clumped.tmp > ${targetDir}/smr.sqtl.clumping/clumped; rm -f ${targetDir}/smr.sqtl.clumping/clumped.tmp

# create smr.assigned file (sort results by LEAD_SNP_P, LEAD_SNP_RSQ, P)
for xqtl in eqtl sqtl; do
	awk -F'\t' 'NR==FNR { snp[$2]=$1"\t"$3"\t"$4"\t"$5; next} $5 in snp { print snp[$5], $0 }' OFS="\t" ${targetDir}/smr.${xqtl}.clumping/clumped ${targetDir}/smr.${xqtl}.filtered.txt > ${targetDir}/smr.${xqtl}.clumping/assigned
	awk -F'\t' '$4 == "A1/A2" { output=$1"\t"$2"\t"$3"\t"$12"/"$13; for(i=5;i<=NF;i++) { output=output"\t"$i }; print output; next } { print }' OFS="\t" ${targetDir}/smr.${xqtl}.clumping/assigned > ${targetDir}/smr.${xqtl}.clumping/assigned.tmp; \mv ${targetDir}/smr.${xqtl}.clumping/assigned.tmp ${targetDir}/smr.${xqtl}.clumping/assigned
	header=$(echo "COND_LOCUS_COUNT"$'\t'"COND_LEAD_SNP_P"$'\t'"COND_LEAD_SNP"$'\t'"COND_LEAD_SNP_KB"$'\t'"COND_LEAD_SNP_RSQ"$'\t'"COND_LEAD_SNP_ALLELES"$'\t'"$(head -1 ${targetDir}/smr.${xqtl}.filtered.txt)")                                                                                     
	awk -F'\t' -v header="${header}" 'BEGIN { print header } NR==1 { next } $2==1 { cond[$5]=$1"\t"$7; next } $1 in cond { print cond[$1], $0 }' OFS='\t' ${conditionalFile} ${targetDir}/smr.${xqtl}.clumping/assigned | sort -k1,1g -k25,25g > ${targetDir}/smr.${xqtl}.filtered.assigned.txt
done

# get summary
scriptDir=$(dirname "$0")
Rscript "${scriptDir}/smr.summarize.R" "${targetDir}/smr.eqtl.filtered.assigned.txt" "${targetDir}/smr.sqtl.filtered.assigned.txt" "${targetDir}/smr.summary.txt"

# clean up
echo "Cleaning up."
rm -f ${targetDir}/smr.subs.txt
rm -f ${targetDir}/smr.input.txt
for xqtl in eqtl sqtl; do
	rm -f ${targetDir}/smr.${xqtl}.sumstats4clumping.txt
	rm -f ${targetDir}/smr.${xqtl}.clumping/*clumped;
	rm -f ${targetDir}/smr.${xqtl}.clumping/*assigned;
	rm -f ${targetDir}/smr.${xqtl}.clumping/*.hh; 
	tar -cvzf ${targetDir}/smr.${xqtl}.clumping.tar.gz --directory=${targetDir} smr.${xqtl}.clumping --remove-files
	tar -cvzf ${targetDir}/smr.${xqtl}.output.tar.gz --directory=${targetDir} smr.${xqtl}.output --remove-files
done
chmod 770 ${targetDir}/*
echo "--- SMR analysis finished. --- "

