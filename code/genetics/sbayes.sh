#!/bin/bash

# ==================
# === run SBayes ===
# ==================

# get arguments
targetDir="${1}" # targetDir="results/gap_gm/discovery/sbayesrc"
sumstats="${2}" # sumstats="results/gap_gm/discovery/gwas/sumstats.txt.gz"
cols="${3}" # set input column names that correspond to header columns SNP A1 A2 freq b se p N | cols="ID,A1,A2,A1_FREQ,BETA,SE,P,N"
ldmFolder="${4}" # ldmFolder="/fast/software/gctb/resources/ukbEUR_HM3"
annotFile="${5}" # annotFile="/fast/software/gctb/resources/annot_baseline2.2.txt"
threads=${6} # threads=100
sbayes=${7} # sbayes="sbayesrc" # options are 'sbayesrc' or 'sbayesr'
imputation=${8} # imputation=1 either 1 or 0
mcmc=${9} # mcmc=1 # write mcmc samples? (required to compute credible sets)

# echo settings
echo $'\n'"--- Run SBayes | Settings ---"
echo "targetDir: ${targetDir}"
echo "sumstats: ${sumstats}"
echo "cols: ${cols}"
echo "ldmFolder: ${ldmFolder}"
echo "annotFile: ${annotFile}"
echo "threads: ${threads}"
echo "sbayes: ${sbayes}"
echo "imputation: ${imputation}"
echo "mcmc: ${mcmc}"$'\n'
sleep 10

# set targetdir and make folder
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# create input file from sumstats
echo "Creating input file."
awk -v cols="${cols}" '
	BEGIN { ncols=split(cols,colnames,",") } # get input column names
    NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } } } # get input column indices
    { for(i=1;i<=ncols;++i) { length(output)==0 ? output = $colidx[i] : output = output"\t"$colidx[i] }; print output; output=""; next }' OFS='\t' <(gzip -dc "${sumstats}") > "${targetDir}/${sbayes}.input.ma"

# imputation
if [ "$imputation" -eq 1 ]; then
    echo "Imputing missing variants."
    gctb \
    --ldm-eigen "${ldmFolder}" \
    --gwas-summary "${targetDir}/${sbayes}.input.ma" \
    --impute-summary \
    --thread ${threads} \
    --out "${targetDir}/${sbayes}"
    gwasFile="${targetDir}/${sbayes}.imputed.ma"
else
    gwasFile="${targetDir}/${sbayes}.input.ma"
fi

# run sbayes
if [ "$sbayes" = "sbayesrc" ]; then    
    if [ "$mcmc" -eq 1 ]; then
        echo "Running SBayesRC with MCMC output"
        sleep 5
        gctb \
        --sbayes RC \
        --ldm-eigen "${ldmFolder}" \
        --gwas-summary "${gwasFile}" \
        --annot "${annotFile}" \
        --write-mcmc-bin \
        --thread ${threads} \
        --out "${targetDir}/${sbayes}"
    else
        echo "Running SBayesRC without MCMC output"
        sleep 5
        gctb \
        --sbayes RC \
        --ldm-eigen "${ldmFolder}" \
        --gwas-summary "${gwasFile}" \
        --annot "${annotFile}" \
        --thread ${threads} \
        --out "${targetDir}/${sbayes}"
    fi
else
    echo "Running SBayesR"
    sleep 5
    gctb \
    --sbayes R \
    --ldm-eigen "${ldmFolder}" \
    --gwas-summary "${gwasFile}" \
    --thread ${threads} \
    --out "${targetDir}/${sbayes}"
fi

# clean up
echo "Cleaning up."
rm -f "${targetDir}/${sbayes}.input.ma"
pigz -f "${targetDir}/${sbayes}.imputed.ma"
pigz -f "${targetDir}/${sbayes}.snpRes"
chmod 770 "${targetDir}/${sbayes}."*

echo "--- Completed: Running SBayes --- "
