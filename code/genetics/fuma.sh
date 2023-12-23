#!/bin/bash

# ==============================================
# === Preparing summary statistics for FUMA ===
# ==============================================

# get arguments
trait="$1" # trait="gap_gm"
targetDir="$2" # targetDir="results/gap_gm/fuma"
sumstatsPLINK="$3" # sumstatsPLINK="results/gap_gm/gwas/sumstats_plink.txt.gz"

# echo settings
echo $'\n'"--- Preparing summary statistics for FUMA ---"
echo "trait: "${trait}
echo "targetDir: "${targetDir}
echo "sumstats: "${sumstatsPLINK}$'\n'

# prepare sumstats for ldhub
echo "Removing chromosomes Y, XY, and MT from PLINK sumstats."
mkdir -p "${targetDir}/"
awk '$1!="Y" && $1!="XY" && $1!="MT" { print }' OFS="\t" <(gzip -dc "${sumstatsPLINK}") > "${targetDir}/fuma.${trait}.txt"

# zip file
echo "gzipping file."
gzip -f "${targetDir}/fuma.${trait}.txt"

# clean up
echo "Cleaning up."
chmod -R 770 $targetDir
echo "--- FUMA preparation finished --- "
