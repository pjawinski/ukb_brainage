#!/bin/bash

# ====================================
# === Genetic correlation analysis ===
# ====================================

# get arguments
trait="$1" # trait="gap_gm"
targetDir="$2" # targetDir="results/gap_gm/rg"
sumstats="$3" # sumstats="results/gap_gm/ldsc/ldsc.gap_gm.sumstats.gz"
nealeDir="$4" # nealeDir="data/ldsc/resources/neale"
nealeManifest="$5" # nealeManifest="data/ldsc/resources/neale_manifest_dropbox.tsv"
showcaseCategories="$6" # showcaseCategories="data/ldsc/resources/neale_showcase_categories.txt"
LDchr="$7" # LDchr="data/ldsc/resources/eur_w_ld_chr/"
threads=$8 # threads=50

# echo settings
echo $'\n'"--- Genetic Correlation Analysis Settings ---"
echo "trait: "${trait}
echo "targetDir: "${targetDir}
echo "sumstats: "${sumstats}
echo "nealeDir: "${nealeDir}
echo "nealeManifest: "${nealeManifest}
echo "showcaseCategories: "${showcaseCategories}
echo "LDchr: "${LDchr}
echo "threads: "${threads}$'\n'

# set targetdir and make folder
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# run genetic correlation analysis
echo "Running genetic correlations."
mkdir -p "${targetDir}/rg.output/"
nealeFiles=$(ls "${nealeDir}")

for nealeFile in ${nealeFiles}; do
	((i=i%threads)); ((i++==0)) && wait
	ldsc.py \
	--rg "${sumstats}","${nealeDir}/${nealeFile}" \
	--ref-ld-chr "${LDchr}" \
	--w-ld-chr "${LDchr}" \
	--out "${targetDir}/rg.output/rg.${nealeFile}" & 
done
wait

# collect results
echo "Collecting results."
echo p1$'\t'p2$'\t'rg$'\t'se$'\t'z$'\t'p$'\t'h2_obs$'\t'h2_obs_se$'\t'h2_int$'\t'h2_int_se$'\t'gcov_int$'\t'gcov_int_se > "${targetDir}/rg.results"
logFiles=$(ls ${targetDir}/rg.output/*.log)
for i in $logFiles; do
    summary="$(grep -m1 -A 2 "Summary of Genetic Correlation Results" $i | awk 'NR==3 { print }' OFS='\t')"
    echo "$summary" | awk '{sub(".*/","",$1); sub(".*/","",$2); print}' OFS='\t' >> "${targetDir}/rg.results"
done

# add trait description and ukb showcase link
echo "Adding trait description and ukb showcase link."
awk -F'\t' 'NR==1 { next } NR==FNR { trait[$5]=$2"\t"$3; next }
	FNR==1 { sub(" ", "\t"); print "trait_description", "showcase", $0; next}
	$2 in trait { sub(" ", "\t"); print trait[$2], $0 }' OFS="\t" "${nealeManifest}" "${targetDir}/rg.results" > "${targetDir}/rg.results.tmp"
\mv "${targetDir}/rg.results.tmp" "${targetDir}/rg.results"

# add categories 
echo "Adding trait categories (ukb dictionary)."
header=trait_description$'\t'category$'\t'showcase$'\t'p1$'\t'p2$'\t'rg$'\t'se$'\t'z$'\t'p$'\t'h2_obs$'\t'h2_obs_se$'\t'h2_int$'\t'h2_int_se$'\t'gcov_int$'\t'gcov_int_se
awk -F'\t' -v header="$header" 'NR==FNR { categ[$1]=$2; next} FNR==1 { print header; next }
	$2 in categ { print $1, categ[$2], $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14; next }
	{ print $1, "N\A", $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}' OFS="\t" "${showcaseCategories}" "${targetDir}/rg.results" > "${targetDir}/rg.results.tmp"
\mv "${targetDir}/rg.results.tmp" "${targetDir}/rg.results"

# clean up
tar -cvzf "${targetDir}/rg.output.tar.gz" --directory="${targetDir}" rg.output --remove-files
chmod -R 770 "${targetDir}"
