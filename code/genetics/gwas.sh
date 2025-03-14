#!/bin/bash

# ========================================
# === Genome-wide association analysis ===
# ========================================

# get arguments
trait="$1" # trait="gap_gm"
traitFile="$2" # traitFile="data/gap_gm/gap_gm.txt"
covsFile="$3" # covsFile="data/gap_gm/covs.txt"
targetDir="$4" # targetDir="results/gap_gm/gwas"
chrFileHandle="$5" # chrFileHandle="data/genetics/chr\$i/imp_mri_qc/chr\${i}_mri_qc"
snpQC="$6" # snpQC="data/genetics/02_r1955_ukb_snp_qc/ukb_snp_qc.txt"
impMFI="$7" # impMFI="data/genetics/05_r1967_ukb_imp_mfi"
maf=$8 # maf=0.01
threads=$9 # threads=50
ancestry="${10}" # ancestry="AFR" | optional: if ancestry is part of chrFileHandle

# echo settings
echo $'\n'"--- GWAS settings ---"
echo "trait: ${trait}"
echo "traitFile: ${traitFile}"
echo "covsFile: ${covsFile}"
echo "targetDir: ${targetDir}"
echo "chrFileHandle: ${chrFileHandle}"
echo "snpQC: ${snpQC}"
echo "impMFI: ${impMFI}"
echo "maf: ${maf}"
echo "threads: ${threads}"$'\n'

# set targetdir and make folder
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# how many covariates
ncovs=$(head -1 "${covsFile}" | mawk '{print NF-2}')

# run gwas
mkdir -p "${targetDir}"/chr
for i in {1..22} X XY Y MT; do
	plink2 \
	--pfile $(eval echo "${chrFileHandle}") \
	--keep <(awk '$3!="NA" { print $1, $2}' ${traitFile}) \
	--glm no-x-sex hide-covar \
	--maf "${maf}" \
	--pheno "${traitFile}" \
	--covar "${covsFile}" \
	--covar-variance-standardize \
	--vif 100000 \
	--freq \
	--threads "${threads}" \
	--out "${targetDir}"/chr/chr"${i}"
done

# merge freq and gwas results
echo "Adding allele frequencies." 
for i in {1..22} X XY Y MT; do (
header="CHR"$'\t'BP$'\t'"ID"$'\t'"A1"$'\t'"A2"$'\t'"A1_FREQ"$'\t'"BETA"$'\t'"SE"$'\t'"P"$'\t'"N"
mawk -v header="$header" 'BEGIN { print header } \
NR==FNR { if ($4==$6) { chr[$3]=$1; pos[$3]=$2; snp[$3]=$3; a1[$3]=$6; a2[$3]=$5; b[$3]=$9; se[$3]=$10; p[$3]=$12; n[$3]=$8; next } else if ($5==$6) { chr[$3]=$1; pos[$3]=$2; snp[$3]=$3; a1[$3]=$6; a2[$3]=$4; b[$3]=$9; se[$3]=$10; p[$3]=$12; n[$3]=$8; next }} \
$2 in chr { if ($4==a1[$2]) { print chr[$2], pos[$2], snp[$2], a1[$2], a2[$2], $5, b[$2], se[$2], p[$2], n[$2] } else { print chr[$2], pos[$2], snp[$2], a1[$2], a2[$2], 1-$5, b[$2], se[$2], p[$2], n[$2] }}' OFS='\t' "${targetDir}"/chr/chr"${i}".*.glm.linear "${targetDir}"/chr/chr"${i}".afreq > "${targetDir}"/chr/chr"${i}"_sumstats.txt
) &
done
wait

# add INFO score (case-sensitive since multiallelic and duplicate variants have not been analyzed)
echo "Adding info scores."
header="CHR"$'\t'BP$'\t'"ID"$'\t'"A1"$'\t'"A2"$'\t'"A1_FREQ"$'\t'"BETA"$'\t'"SE"$'\t'"P"$'\t'"N"$'\t'"INFO"
for i in {1..22} X XY; do (
mawk -v header="$header" 'BEGIN { print header } NR==FNR { info[$2]=$8; next } $3 in info { print $0, info[$3] }' OFS='\t' "${impMFI}"/ukb_mfi_chr"${i}"_v3.txt "${targetDir}"/chr/chr${i}_sumstats.txt > "${targetDir}"/chr/chr"${i}"_sumstats.tmp.txt
mv -f "${targetDir}"/chr/chr"${i}"_sumstats.tmp.txt "${targetDir}"/chr/chr"${i}"_sumstats.txt) &
done
wait

# add TYPED (case-sensitive since multiallelic and duplicate variants have been removed before analyses)
echo "Adding TYPED variable."
header="CHR"$'\t'BP$'\t'"ID"$'\t'"A1"$'\t'"A2"$'\t'"A1_FREQ"$'\t'"BETA"$'\t'"SE"$'\t'"P"$'\t'"N"$'\t'"INFO"$'\t'"TYPED"
for i in {1..22} X XY; do (
mawk -v header="$header" 'BEGIN { print header } NR==FNR { id[$1]; next } FNR==1 { next } $3 in id { print $0, 1; next } { print $0, 0;}' OFS='\t' ${snpQC} ${targetDir}/chr/chr${i}_sumstats.txt > ${targetDir}/chr/chr${i}_sumstats.tmp.txt
mv -f "${targetDir}"/chr/chr"${i}"_sumstats.tmp.txt "${targetDir}"/chr/chr"${i}"_sumstats.txt) &
done
wait

# concatenate plink results
echo "Concatenating plink results."
header="CHROM"$'\t'BP$'\t'"ID"$'\t'"REF"$'\t'"ALT"$'\t'"A1"$'\t'"TEST"$'\t'"OBC_CT"$'\t'"BETA"$'\t'"SE"$'\t'"T_STAT"$'\t'"P"
mawk -v header="$header" 'BEGIN { print header } FNR>1' "${targetDir}"/chr/chr{{1..22},X,XY,Y,MT}.*.glm.linear > "${targetDir}"/sumstats_plink.txt

# concatenate sumstats results
echo "Concatenating custom sumstats."
header="CHR"$'\t'BP$'\t'"ID"$'\t'"A1"$'\t'"A2"$'\t'"A1_FREQ"$'\t'"BETA"$'\t'"SE"$'\t'"P"$'\t'"N"$'\t'"INFO"$'\t'"TYPED"
mawk -v header="$header" 'BEGIN { print header } FNR>1 && ($1=="Y" || $1=="MT") { print $0, 1, 1; next } FNR>1 { print }' OFS='\t' "${targetDir}"/chr/chr{{1..22},X,XY,Y,MT}_sumstats.txt > "${targetDir}"/sumstats.txt
		
# add Z column
echo "Adding Z score column."
header="CHR"$'\t'BP$'\t'"ID"$'\t'"A1"$'\t'"A2"$'\t'"A1_FREQ"$'\t'"BETA"$'\t'"SE"$'\t'"Z"$'\t'"P"$'\t'"N"$'\t'"INFO"$'\t'"TYPED"
mawk -v header="$header" 'BEGIN { print header } FNR>1 { print $1, $2, $3, $4, $5, $6, $7, $8, $7/$8, $9, $10, $11, $12 }' OFS='\t' "${targetDir}"/sumstats.txt > "${targetDir}"/sumstats.tmp.txt
mv -f "${targetDir}"/sumstats.tmp.txt "${targetDir}"/sumstats.txt

# add ETA2 column 
echo "Adding ETA2 column."
header="CHR"$'\t'BP$'\t'"ID"$'\t'"A1"$'\t'"A2"$'\t'"A1_FREQ"$'\t'"BETA"$'\t'"SE"$'\t'"Z"$'\t'"P"$'\t'"N"$'\t'"ETA2"$'\t'"INFO"$'\t'"TYPED"
mawk -v header="$header" -v ncovs="${ncovs}" 'BEGIN { print header } FNR>1 { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $9^2/($9^2+($11-2-ncovs)), $12, $13 }' OFS='\t' ${targetDir}/sumstats.txt > ${targetDir}/sumstats.tmp.txt
mv -f "${targetDir}"/sumstats.tmp.txt "${targetDir}"/sumstats.txt

# remove results with NA (observed 11 snps in low sample-size replication cohort MID)
echo "Remove and log variants with NA."
NaN=$(mawk '$10 == "NA"' OFS='\t' "${targetDir}"/sumstats.txt) 
if [[ -z "$NaN" ]]; then
	echo "No variants with invalid results (NA) found." | tee ${targetDir}/chr/NAvariants.log
else
	mawk '$10 != "NA"' OFS='\t' "${targetDir}"/sumstats.txt > "${targetDir}"/sumstats.tmp.txt; mv -f "${targetDir}"/sumstats.tmp.txt "${targetDir}"/sumstats.txt
	mawk '$12 != "NA"' OFS='\t' "${targetDir}"/sumstats_plink.txt > "${targetDir}"/sumstats_plink.tmp.txt; mv -f "${targetDir}"/sumstats_plink.tmp.txt "${targetDir}"/sumstats_plink.txt
	echo "$(head -1 "${targetDir}/sumstats.txt")"$'\n'"${NaN}" > "${targetDir}"/chr/NAvariants.log
	echo "$(echo "${NaN}" | wc -l) variants with invalid results (NA) removed."
fi

# clean up folder
echo "Cleaning up folder."
rm -rf "${targetDir}"/chr/*.afreq
rm -rf "${targetDir}"/chr/*.linear
rm -rf "${targetDir}"/chr/*.txt
mv "${targetDir}"/chr "${targetDir}"/log
tar -czf "${targetDir}"/sumstats.log.tar.gz --directory="${targetDir}" log --remove-files
echo "G-zipping files."
pigz -f "${targetDir}"/sumstats.txt > "${targetDir}"/sumstats.txt.gz &
pigz -f "${targetDir}"/sumstats_plink.txt > "${targetDir}"/sumstats_plink.txt.gz
wait
chmod 770 "${targetDir}"/*
echo "--- GWAS finished. --- "

