#!/bin/bash

# ============================================================
# === Run GWASs in ancestry-stratified replication samples ===
# ============================================================

# get arguments
trait="${1}" # trait="t1.gap_gm"
targetDir="${2}" # targetDir="results/gap_gm/replicate"
chrFileHandle="${3}" # chrFileHandle="data/genetics/chr\$i/imp_mri_qc_${ancestry}/chr\${i}_mri_qc"
data="${4}" # data="data/replicate/mri_repl_pan.txt"
covs="${5}" # covs="sex,t1.age,t1.age2,t1.ac1,t1.ac2,t1.ac3,t1.TIV,array,PanC{1..20}"
ancestries="${6}" # ancestries="AFR,AMR,CSA,EAS,EUR,MID"
ancestriesCol="${7}" # ancestriesCol="pan"
snpQC="${8}" # snpQC="data/genetics/02_r1955_ukb_snp_qc/ukb_snp_qc.txt"
impMFI="${9}" # impMFI="data/genetics/05_r1967_ukb_imp_mfi"
maf=${10} # maf=0.01
threads=${11} # threads=50

# echo settings
echo $'\n'"--- Replicate GWAS discoveries ---"
echo "trait: ${trait}"
echo "targetDir: ${targetDir}"
echo "chrFileHandle: ${chrFileHandle}"
echo "data: ${data}"
echo "covs: ${covs}"
echo "ancestries: ${ancestries}"
echo "ancestriesCol: ${ancestriesCol}"
echo "snpQC: ${snpQC}"
echo "impMFI: ${impMFI}"
echo "maf: ${maf}"
echo "threads: ${threads}"$'\n'

# covs: replace PanC{1..n} with PanC1,PanC2,..,PanCn
PCshort=$(echo "${covs}" | awk -F',' '$NF ~ /../ { print $NF }')
PClong=$(eval "echo $PCshort" | awk '{ gsub(/ /,","); print }')
tmp=$(echo "${covs}" | sed -e "s/$PCshort/$PClong/g")
covs=${tmp}

# loop over ancestries
ancestries=$(echo "${ancestries}" | sed 's/,/ /g')
for ancestry in ${ancestries}; do

	# create folder
	mkdir -p "${targetDir}/${ancestry}"

	# prepare covsFile and traitFile
	scriptDir=$(dirname "$0")
	Rscript "${scriptDir}"/replicate.R "${trait}" "${targetDir}/${ancestry}" "${data}" "${covs}" "${ancestry}" "${ancestriesCol}"

	# run gwas
	"${scriptDir}"/gwas.sh "${trait}" "${targetDir}/${ancestry}/trait.txt" "${targetDir}/${ancestry}/covs.txt" "${targetDir}/${ancestry}" "${chrFileHandle}" "${snpQC}" "${impMFI}" "${maf}" "${threads}" "${ancestry}"

	# clean destination folder
	rm -f "${targetDir}"/"${ancestry}"/trait.txt
	rm -f "${targetDir}"/"${ancestry}"/covs.txt

# end loop over ancestries
done
echo "--- Replication GWAS finished. --- "

