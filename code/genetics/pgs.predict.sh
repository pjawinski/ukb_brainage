#!/bin/bash

# =====================
# === Calculate PGS ===
# =====================

# get arguments
weightFile="${1}" # weightFile="results/gap_gm/discovery/sbayesrc/sbayesrc.snpRes.gz"
pgenFileHandler="${2}" # pgenFileHandler='data/genetics/chr${i}/imp_mri_qc_EUR/chr${i}_mri_qc'
idCol="${3}" # idCol="Name"
a1Col="${4}" # a1Col="A1"
effectCol="${5}" # effectCol="A1Effect"
maf="${6}" # maf=0.01
threads=${7} # threads=100
outFile="${8}" # outFile="data/genetics/prs/imp_mri_qc_EUR/gap_gm.sbayesrc"

# echo settings
echo $'\n'"--- Calculate PGS | Settings ---"
echo "weightFile: ${weightFile}"
echo "pgenFileHandler: ${pgenFileHandler}"
echo "idCol: ${idCol}"
echo "a1Col: ${a1Col}"
echo "effectCol: ${effectCol}"
echo "maf: ${maf}"
echo "threads: ${threads}"
echo "outFile: ${outFile}"$'\n'

# set targetdir and make folder
targetDir="$(dirname "${outFile}")"
mkdir -p "${targetDir}"
mkdir -p "${outFile}.chr"
targetDir="$(readlink -f "${targetDir}")"

# get column numbers
idColnum=$(awk -v idCol="${idCol}" ' NR==1 { for(i=1;i<=NF;++i) { if($i==idCol) { print i; exit 1} } }' <(zcat "${weightFile}"))
a1Colnum=$(awk -v a1Col="${a1Col}" ' NR==1 { for(i=1;i<=NF;++i) { if($i==a1Col) { print i; exit 1} } }' <(zcat "${weightFile}"))
effectColnum=$(awk -v effectCol="${effectCol}" ' NR==1 { for(i=1;i<=NF;++i) { if($i==effectCol) { print i; exit 1} } }' <(zcat "${weightFile}"))

# calculate polygenic score
echo "Calculating scores for individual chromosomes"
for i in {1..22}; do
    plink2 \
        --pfile $(eval echo "${pgenFileHandler}") \
        --maf "${maf}" \
        --score "${weightFile}" "${idColnum}" "${a1Colnum}" "${effectColnum}" header cols=+scoresums \
        --threads "${threads}" \
        --out "${outFile}.chr/chr${i}"
done

awk 'NR==1 { print "FID", "IID", "CHR_COUNT", "NMISS_ALLELE_CT", "SCORE1_SUM"; next }
     NR==FNR { id[$1]=$1; chrcount[$1]=1; nmiss[$1]=$3; score[$1]=$6; next }
     FNR==1 { next }
     { chrcount[$1]++; nmiss[$1]=nmiss[$1]+$3; score[$1]=score[$1]+$6 }
     END { for (i in score) { print id[i], id[i], chrcount[i], nmiss[i], score[i] } }' OFS='\t' "${outFile}.chr/"*.sscore \
     > ${outFile}.score

# clean up
echo "Cleaning up."
tar -cvzf "${outFile}.chr.tar.gz" --directory="${targetDir}" "$(basename ${outFile}.chr)" --remove-files
chmod 770 "${outFile}"*
echo "--- Completed: Calculate PGS --- "


