#!/bin/bash

# =====================================
# ======= combine ldsc results ========
# =====================================

# get arguments
ldscFiles="$1" # ldscFiles="results/gap_gm/ldsc/ldsc.h2.results results/gap_wm/ldsc/ldsc.h2.results results/gap_gwm/ldsc/ldsc.h2.results"
outFile="$2" # outFile="results/combined/ldsc.h2.results"

# echo settings
echo $'\n'"--- Combining LDSC output files ---"
echo "ldscFiles: "${ldscFiles}
echo "outFile: "${outFile}

# combine ldsc.h2.results
ldscFiles=($ldscFiles)
echo trait$'\t'snp_count$'\t'h2$'\t'h2_se$'\t'lambda_GC$'\t'chi2$'\t'intercept$'\t'intercept_se$'\t'ratio$'\t'ratio_se > ${outFile}
for (( i=0; i<${#ldscFiles[@]}; i++ )); do
	awk -F '\t' 'FNR == 1 { next } { gsub(/\(/,""); gsub(/\)/,""); gsub(" ","\t"); print }' OFS='\t' "${ldscFiles[i]}" >> ${outFile}
done

# clean up
chmod 770 "${outFile}"
echo "--- Combining LDSC output files completed. ---"
