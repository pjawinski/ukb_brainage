#!/bin/bash

# ===========================
# === LD Score Regression ===
# ===========================

# get arguments
trait="$1" # trait="gap_gm"
targetDir="$2" # targetDir="results/gap_gm/ldsc"
sumstats="$3"  # sumstats="results/gap_gm/gwas/sumstats.txt.gz"
LDsnplist="$4" # LDsnplist="data/ldsc/resources/w_hm3.noMHC.snplist"
LDchr="$5" # LDchr="data/ldsc/resources/eur_w_ld_chr/"
LDbaseline="$6" # LDbaseline="data/ldsc/resources/baselineLD/baselineLD."
LDweights="$7" # LDweights="data/ldsc/resources/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
LDfrqfiles="$8" # LDfrqfiles="data/ldsc/resources/1000G_Phase3_frq/1000G.EUR.QC."

# echo settings
echo $'\n'"--- LD Score Regression Settings ---"
echo "trait: "${trait}
echo "targetDir: "${targetDir}
echo "sumstats: "${sumstats}
echo "LDsnplist: "${LDsnplist}
echo "LDchr: "${LDchr}
echo "LDbaseline: "${LDbaseline}
echo "LDweights: "${LDweights}
echo "LDfrqfiles: "${LDfrqfiles}$'\n'

# set targetdir and make folder
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# prepare sumstats
echo "Preparing sumstats."
header="CHR"$'\t'"POS"$'\t'"SNP"$'\t'"A1"$'\t'"A2"$'\t'"BETA"$'\t'"SE"$'\t'"P"$'\t'"N"
awk -v header="$header" 'BEGIN { print header } NR > 1 { print $1, $2, $3, $4, $5, $7, $8, $10, $11}' OFS="\t" <(gzip -dc "${sumstats}") > "${targetDir}/ldsc.${trait}.prep"

# munge sumstats
echo "Munging sumstats."
munge_sumstats.py \
--sumstats "${targetDir}/ldsc.${trait}.prep" \
--merge-alleles "${LDsnplist}" \
--out "${targetDir}/ldsc.${trait}"
chmod -R 770 "${targetDir}"

# ld score regression
echo "Running LD Score regression."
ldsc.py \
--h2 "${targetDir}/ldsc.${trait}.sumstats.gz" \
--ref-ld-chr "${LDchr}" \
--w-ld-chr "${LDchr}" \
--out "${targetDir}/ldsc.h2"

# get summary
snps="$(cat "${targetDir}/ldsc.h2.log" | grep -m1 "After merging with regression SNP LD, " | awk '{ print $7 }')"
h2="$(cat "${targetDir}/ldsc.h2.log" | grep -m1 "Total Observed scale h2: " | awk '{ print $5, $6 }')"
Lambda="$(cat "${targetDir}/ldsc.h2.log" | grep -m1 "Lambda GC: " | awk '{ print $3 }')"
Chi2="$(cat "${targetDir}/ldsc.h2.log" | grep -m1 "Mean Chi^2: " | awk '{ print $3 }')"
Intercept="$(cat "${targetDir}/ldsc.h2.log" | grep -m1 "Intercept: " | awk '{ print $2, $3 }')"
Ratio="$(cat "${targetDir}/ldsc.h2.log" | grep -m1 "Ratio: " | awk '{ print $2, $3 }')"
echo trait$'\t'snp_count$'\t'h2$'\t'lambda_GC$'\t'chi2$'\t'intercept$'\t'ratio > "${targetDir}/ldsc.h2.results"
echo ${trait}$'\t'${snps}$'\t'${h2}$'\t'${Lambda}$'\t'${Chi2}$'\t'${Intercept}$'\t'${Ratio} >> "${targetDir}/ldsc.h2.results"

# run partitioned ld score regression
echo "Running Partitioned LD Score regression."
ldsc.py \
--h2 "${targetDir}/ldsc.${trait}.sumstats.gz" \
--ref-ld-chr "${LDbaseline}" \
--w-ld-chr "${LDweights}" \
--overlap-annot \
--frqfile-chr "${LDfrqfiles}" \
--out "${targetDir}/ldsc.partitioned"

# clean up
rm -f "${targetDir}/ldsc.${trait}.prep"
chmod -R 770 "${targetDir}"

