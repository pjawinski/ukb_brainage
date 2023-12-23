#!/bin/bash

# =====================
# === GSMR analysis ===
# =====================

# get arguments
trait="${1}" # trait="gap_gm"
targetDir="${2}" # targetDir="results/${trait}/gsmr"
sumstats="${3}" # sumstats="results/${trait}/gwas/sumstats.txt.gz"
sumstatsCols="${4}" # sumstatsCols="ID,A1,A2,A1_FREQ,BETA,SE,P,N" | same order as in header variable, see 'prepare sumstats' below
sampleFile="${5}" # sampleFile=data/${trait}/${trait}.txt
mrFiles="${6}" # mrFiles=$(echo $(ls data/sumstats/mr/mr_*.gz) | sed 's/ /,/g')
mrLabels="${7}" # mrLabels=$(echo $(for i in $(echo $mrFiles | sed 's/,/ /g'); do echo $i | sed 's%.*/%%g' | sed 's%\..*%%g'; done) | sed 's/ /,/g')
chrFilehandler="${8}" # chrFilehandler='data/genetics/chr${i}/imp_mri_qc/bed/chr${i}_mri_qc'

# echo settings
echo $'\n'"--- GSMR analysis ---"
echo "trait: "${trait}
echo "targetDir: "${targetDir}
echo "sumstats: "${sumstats}
echo "sumstatsCols: "${sumstatsCols}
echo "sampleFile: "${sumstatsCols}
echo "mrFiles: "${mrFiles}
echo "mrLabels: "${mrFiles}
echo "chrFilehandler: "${chrFilehandler}

# set targetdir and make folder
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# prepare sumstats
echo "Preparing sumstats."
header="SNP"$'\t'"A1"$'\t'"A2"$'\t'"freq"$'\t'"b"$'\t'"se"$'\t'"p"$'\t'"N"
awk -v header="${header}" -v cols="${sumstatsCols}" '
	BEGIN { ncols=split(cols,colnames,","); print header }
    NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } } }
    NR>1 { output=$colidx[1]; for(i=2;i<=ncols;++i) { output=output"\t"$colidx[i] }; print output}
	' OFS="\t" <(gzip -dc "${sumstats}") > "${targetDir}/gsmr.${trait}.prep"

# prepare sampleFile
echo "Preparing sample file."
awk '{ print $1, $2 }' OFS='\t' ${sampleFile} > ${targetDir}/gsmr.sample.txt

# create list of LD reference data files
> ${targetDir}/gsmr.ref.data.txt
for i in {1..22} X Y MT; do
	eval echo $chrFilehandler >> ${targetDir}/gsmr.ref.data.txt
done

# create gsmr exposure file
mrFiles=($(echo $mrFiles | sed 's/,/ /g'))
mrLabels=($(echo $mrLabels | sed 's/,/ /g'))
> ${targetDir}/gsmr.exposure.txt
for i in $(seq 0 $((${#mrFiles[@]}-1))); do
	echo "${mrLabels[i]} ${mrFiles[i]}" >> ${targetDir}/gsmr.exposure.txt
done

# create gsmr outcome file
> ${targetDir}/gsmr.outcome.txt
echo ${trait} ${targetDir}/gsmr.${trait}.prep >> ${targetDir}/gsmr.outcome.txt

# run mendelian randomisation (https://cnsgenomics.com/software/gcta/#Mendelianrandomisation)
gcta64 \
--mbfile ${targetDir}/gsmr.ref.data.txt \
--keep ${targetDir}/gsmr.sample.txt \
--gsmr-file ${targetDir}/gsmr.exposure.txt ${targetDir}/gsmr.outcome.txt \
--gsmr-direction 2 \
--effect-plot \
--gwas-thresh 5e-8 \
--clump-r2 0.05 \
--diff-freq 1 \
--gsmr-snp-min 10 \
--out ${targetDir}/gsmr

# clean up
rm -f "${targetDir}/gsmr.${trait}.prep"
rm -f "${targetDir}/gsmr.sample.txt"
rm -f "${targetDir}/gsmr.exposure.txt"
rm -f "${targetDir}/gsmr.outcome.txt"
chmod -R 770 "${targetDir}"
echo $'\n'"--- Completed: GSMR analysis ---"

