#!/bin/bash

# =============================
# === Heritability analysis ===
# =============================

# get arguments
trait="$1" # trait="gap_gm"
traitFile="$2" # traitFile="data/gap_gm/gap_gm.txt"
covsFile="$3" # covsFile="data/gap_gm/covs.txt"
targetDir="$4" # targetDir="results/gap_gm/greml"
covsDiscr="$5" # covsDiscr="sex,ac1,ac2,array"
covsQuant="$6" # covsQuant="age,age2,TIV,PC{1..20}"
grm="$7" # grm="data/grm/ukb_merged"
threads=$8 # threads=50

# echo settings
echo $'\n'"--- Heritability Analysis Settings ---"
echo "trait: "${trait}
echo "traitFile: "${traitFile}
echo "covsFile: "${covsFile}
echo "targetDir: "${targetDir}
echo "covsDiscr: "${covsDiscr}
echo "covsQuant: "${covsQuant}
echo "grm: "${grm}
echo "threads: "${threads}$'\n'

# set targetdir and make folder
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# covsQuant: replace {PC1..n} with PC1,PC2,..,PCn
PCshort=$(echo $covsQuant | awk -F',' '$NF ~ /../ { print $NF }')
PClong=$(eval "echo $PCshort" | awk '{ gsub(/ /,","); print }')
temp=$(echo $covsQuant | sed -e "s/$PCshort/$PClong/g")
covsQuant=$temp

# get covariate files
echo "Creating covariate files."
awk -F '\t' -v covs="${covsQuant}" 'BEGIN { ncovs=split(covs,covnames,",") }
    NR==1 {output=$1"\t"$2; for(i=1;i<=ncovs;++i) { for(j=3;j<=NF;++j) { if($j==covnames[i]) { colnums[i]=j; output=output"\t"$j } } }; print output }
    NR>1 { output=$1"\t"$2; for(i=1;i<=ncovs;++i) { output=output"\t"$colnums[i] }; print output}' "${covsFile}" > "${targetDir}/greml.covsQuant.txt"
awk -F '\t' -v covs="${covsDiscr}" 'BEGIN { ncovs=split(covs,covnames,",") }
    NR==1 {output=$1"\t"$2; for(i=1;i<=ncovs;++i) { for(j=3;j<=NF;++j) { if($j==covnames[i]) { colnums[i]=j; output=output"\t"$j } } }; print output }
    NR>1 { output=$1"\t"$2; for(i=1;i<=ncovs;++i) { output=output"\t"$colnums[i] }; print output}' "${covsFile}" > "${targetDir}/greml.covsDiscr.txt"

# conduct greml analysis
gcta64 --grm "${grm}" \
--pheno "${traitFile}" \
--reml \
--covar "${targetDir}/greml.covsDiscr.txt" \
--qcovar "${targetDir}/greml.covsQuant.txt" \
--out "${targetDir}/greml" \
--thread-num ${threads}

# calculate exact pvalue of heritability estimate
scriptDir=$(dirname "$0")
Rscript "${scriptDir}/greml.pval.R" "${targetDir}"

# clean up
echo "Cleaning up."
head -1 "${targetDir}/greml.covsDiscr.txt" > "${targetDir}/greml.covsDiscr.txt.tmp"; \mv "${targetDir}/greml.covsDiscr.txt.tmp" "${targetDir}/greml.covsDiscr.txt"
head -1 "${targetDir}/greml.covsQuant.txt" > "${targetDir}/greml.covsQuant.txt.tmp"; \mv "${targetDir}/greml.covsQuant.txt.tmp" "${targetDir}/greml.covsQuant.txt"
chmod -R 770 "${targetDir}/"
echo "--- Heritability analysis finished. --- "


