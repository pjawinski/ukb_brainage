#!/bin/bash

# ==================================
# === run GWAS catalog screening ===
# ==================================

# get arguments
trait="$1" # trait="gap_gm"
targetDir="$2" # targetDir="results/${trait}/catalog"
conditionalFile="$3" # conditionalFile="results/${trait}/conditional/conditional.cleaned.tophits.annovar.txt"
catalogFile="$4" # catalogFile="data/gwas_catalog/gwas_catalog_v1.0-associations_e105_r2021-12-21.tsv"
gwasSig=$5 # gwasSig=5E-8
catalogSig=$6 # catalogSig=5E-8

# echo settings
echo $'\n'"--- GWAS Catalog Settings ---"
echo "trait: "${trait}
echo "targetDir: "${targetDir}
echo "Conditional file: "${conditionalFile}
echo "GWAS catalog file: "${catalogFile}
echo "Significance threshold: "${gwasSig}
echo "GWAS catalog gwasSignificance threshold: "${catalogSig}$'\n'

# create target directory
mkdir -p ${targetDir}

# get conditional tophits and drop loci due to 'non-independence' (pvalue has changed by two orders of magnitude and variation is with 10 MB window of another more gwasSignificant variation)
# awk -F '\t' 'NR==1 { next } NR==FNR { dropped[$2]; next } FNR==1 { print; next } !($6 in dropped) { print }' ${conditionalDropped} <(gzip -dc ${conditionalFile}) > ${targetDir}/conditional.txt

# run catalog screening
header=$(echo "$(head -1 ${conditionalFile})")
header=$(echo "$header"$'\t'"CATALOG_PUBMEDID"$'\t'"CATALOG_FIRST_AUTHOR"$'\t'"CATALOG_DATE"$'\t'"CATALOG_DISEASE/TRAIT"$'\t'"CATALOG_REPORTED_GENE(S)"$'\t'"CATALOG_RISK_ALLELE"$'\t'"CATALOG_RISK_ALLELE_FREQUENCY"$'\t'"CATALOG_OR_or_BETA"$'\t'"CATALOG_P-VALUE")
awk -v header="$header" -v gwasSig="$gwasSig" -v catalogSig="$catalogSig" -F'\t' 'BEGIN { print header }
	NR==FNR && $17 <= gwasSig && $31 <= gwasSig { snp[$5]=$0; snp2[$5]=$0"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"; next}
	NR==FNR { next }
	$22 in snp && $28 < catalogSig { print snp[$22], $2, $3, $4, $8, $14, $21, $27, $31, $28; delete snp2[$22]} END { for (i in snp2) print snp2[i] }' OFS="\t" ${conditionalFile} ${catalogFile} \
	| sort -k1,1n -k2,2n > ${targetDir}/catalog.txt

# finalise output file - drop some columns & remove variations without catalog match (except for index variations)
scriptDir=$(dirname "$0")
Rscript "${scriptDir}/catalog.output.R" "${targetDir}/catalog.txt"

# finish analysis
chmod -R 770 "${targetDir}"
echo "--- GWAS catalog screening finished. ---"
