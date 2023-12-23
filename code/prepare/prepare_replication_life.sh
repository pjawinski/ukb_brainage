#!/bin/bash

# ========================================
# === Harmonize LIFE with UKB sumstats ===
# ========================================

# get arguments
sumstatsLIFE="$1" # sumstatsLIFE="data/life/gwas/s303_1_GWAS_brain_noCovStats.brainage_gap_gm_stack.glm.linear"
ukbSNPs="$2" # ukbSNPs="data/genetics/qc_snplist/ukb_info_08_nodups.txt"
targetDir="$3" # targetDir="${targetDir}/"

# echo settings
echo $'\n'"--- Harmonize LIFE with UKB sumstats ---"
echo "sumstatsLIFE: "${sumstatsLIFE}
echo "ukbSNPs: "${ukbSNPs}
echo "targetDir: "${targetDir}$'\n'

# set targetdir and make folder
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# match by chr_bp_ref_alt or chr_bp_alt_ref
ukbSNPs="data/genetics/qc_snplist/ukb_info_08_nodups.txt"
sumstatsLIFE="data/life/gwas/s303_1_GWAS_brain_noCovStats.brainage_gap_gm_stack.glm.linear"

header=ID$'\t'$(head -1 ${sumstatsLIFE} | sed 's/ID/ID_LIFE/g' | sed 's/[#]//g')$'\t'A2$'\t'Z$'\t'ETA2
mkdir -p ${targetDir}/chr
(for chr in {1..22} X; do (

	# match variant ids
	if [ ${chr} -eq 1 ]; then echo "Matching variants by chr_bp_ref_alt or chr_bp_alt_ref."; fi
	mawk -v chr=${chr} 'FNR==1 { next }
		  NR==FNR && $1==chr { id[$1"_"$2"_"$4"_"$5]=$3; next }
		  NR==FNR { next }
		  $1==chr && $1"_"$2"_"$4"_"$5 in id { print id[$1"_"$2"_"$4"_"$5], $0; next }
		  $1==chr && $1"_"$2"_"$5"_"$4 in id { print id[$1"_"$2"_"$5"_"$4], $0; next }
		  ' OFS='\t' ${ukbSNPs} ${sumstatsLIFE} > ${targetDir}/chr/chr${chr}.txt

	# remove potential duplicates
	if [ ${chr} -eq 1 ]; then echo "Removing potential duplicates."; fi
	mawk 'NR==FNR && seen[$1]++ { dupVar[$1]; next } 
	  	  NR==FNR { next }
	  	  !($1 in dupvar) { print }
	  	  ' OFS='\t' ${targetDir}/chr/chr${chr}.txt ${targetDir}/chr/chr${chr}.txt > ${targetDir}/chr/chr${chr}.tmp.txt
	  	  \mv ${targetDir}/chr/chr${chr}.tmp.txt ${targetDir}/chr/chr${chr}.txt

	# deduce A2 from REF ALT A1
	if [ ${chr} -eq 1 ]; then echo "Deducing A2 allele from REF, ALT, and A1."; fi
	mawk '$5==$7 { print $0, $6; next } 
		  $6==$7 { print $0, $5; next }
		  ' OFS='\t' ${targetDir}/chr/chr${chr}.txt ${targetDir}/chr/chr${chr}.txt > ${targetDir}/chr/chr${chr}.tmp.txt
	  	  \mv ${targetDir}/chr/chr${chr}.tmp.txt ${targetDir}/chr/chr${chr}.txt

	# calculate Z and ETA2 (with sex, age, age2, PC1-4 serving as covariates)
	if [ ${chr} -eq 1 ]; then echo "Calculating Z and ETA2."; fi
	mawk '{ print $0, $12/$13, ($12/$13)^2/(($12/$13)^2+($11-2-7)); next } 
		  ' OFS='\t' ${targetDir}/chr/chr${chr}.txt ${targetDir}/chr/chr${chr}.txt > ${targetDir}/chr/chr${chr}.tmp.txt
	  	  \mv ${targetDir}/chr/chr${chr}.tmp.txt ${targetDir}/chr/chr${chr}.txt
	) & 
done
wait)

# merge chromosome output
echo "Merging chromosome-wise output."
mawk -v header="${header}" 'BEGIN { print header } { print }' OFS='\t' ${targetDir}/chr/chr{{1..22},X}.txt \
	> ${targetDir}/sumstats.txt

# rename columns
echo "Harmonizing column names with UKB sumstats."
header=ID$'\t'CHR$'\t'BP$'\t'ID_LIFE$'\t'REF$'\t'ALT$'\t'A1$'\t'A1_FREQ$'\t'INFO$'\t'TEST$'\t'N$'\t'BETA$'\t'SE$'\t'L95$'\t'U95$'\t'T_STAT$'\t'P$'\t'A2$'\t'Z$'\t'ETA2
mawk -v header="${header}" 'BEGIN { print header } NR==1 { next } { print }' OFS='\t' ${targetDir}/sumstats.txt > ${targetDir}/sumstats.tmp.txt
	\mv ${targetDir}/sumstats.tmp.txt ${targetDir}/sumstats.txt

# only keep columns also present in ukb sumstats
echo "Re-ordering and selecting columns in accordance with UKB sumstats."
mawk '{ print $2,$3,$1,$7,$18,$8,$12,$13,$19,$17,$11,$9 }' OFS='\t' ${targetDir}/sumstats.txt > ${targetDir}/sumstats.tmp.txt
	\mv ${targetDir}/sumstats.tmp.txt ${targetDir}/sumstats.txt

# clean up
echo "Cleaning up."
rm -rf ${targetDir}/chr
pigz -f ${targetDir}/sumstats.txt

