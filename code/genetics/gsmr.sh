#!/bin/bash

# =====================
# === GSMR analysis ===
# =====================

# get arguments
trait="${1}" # trait="gap_gm"
outFile="${2}" # targetDir="results/${trait}/gsmr"
sumstats="${3}" # sumstats="results/${trait}/gwas/sumstats.txt.gz"
sumstatsCols="${4}" # sumstatsCols="ID,A1,A2,A1_FREQ,BETA,SE,P,N" | same order as in header variable, see 'prepare sumstats' below
sampleFile="${5}" # sampleFile=data/${trait}/${trait}.txt
mrFiles="${6}" # mrFiles=$(echo $(ls data/sumstats/mr/mr_*.gz) | sed 's/ /,/g')
mrLabels="${7}" # mrLabels=$(echo $(for i in $(echo $mrFiles | sed 's/,/ /g'); do echo $i | sed 's%.*/%%g' | sed 's%\..*%%g'; done) | sed 's/ /,/g')
chrFilehandler="${8}" # chrFilehandler='data/genetics/chr${i}/imp_mri_qc/bed/chr${i}_mri_qc'
clumpr2="${9}" # default: clumpr2=0.05 (default)
clumpkb="${10}" # default: clumpkb=1000 (in kbp)
gsmr2beta="${11}" # either 1 (TRUE) or 0 (FALSE)

# echo settings
echo $'\n'"--- GSMR analysis ---"
echo "trait: ${trait}"
echo "outFile: ${outFile}"
echo "sumstats: ${sumstats}"
echo "sumstatsCols: ${sumstatsCols}"
echo "sampleFile: ${sampleFile}"
echo "mrFiles: ${mrFiles}"
echo "mrLabels: ${mrFiles}"
echo "chrFilehandler: ${chrFilehandler}"
echo "clumpr2: ${clumpr2}"
echo "clumpkb: ${clumpkb}"
echo "gsmr2beta: ${gsmr2beta}"

# set targetdir and make folder
targetDir=$(dirname "${outFile}")
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# prepare sumstats
echo "Preparing sumstats."
header="SNP"$'\t'"A1"$'\t'"A2"$'\t'"freq"$'\t'"b"$'\t'"se"$'\t'"p"$'\t'"N"
awk -v header="${header}" -v cols="${sumstatsCols}" '
	BEGIN { ncols=split(cols,colnames,","); print header }
    NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } } }
    NR>1 { output=$colidx[1]; for(i=2;i<=ncols;++i) { output=output"\t"$colidx[i] }; print output}
	' OFS="\t" <(gzip -dc "${sumstats}") > "${outFile}"."${trait}".prep

# prepare sampleFile
echo "Preparing sample file."
awk '{ print $1, $2 }' OFS='\t' "${sampleFile}" > "${outFile}".sample.txt

# create list of LD reference data files
> "${outFile}".ref.data.txt
for i in {1..22} X Y MT; do
	eval echo "${chrFilehandler}" >> "${outFile}".ref.data.txt
done

# create gsmr exposure file
mrFiles=($(echo "${mrFiles}" | sed 's/,/ /g'))
mrLabels=($(echo "${mrLabels}" | sed 's/,/ /g'))
> "${outFile}".exposure.txt
for i in $(seq 0 $((${#mrFiles[@]}-1))); do
	echo "${mrLabels[i]} ${mrFiles[i]}" >> "${outFile}".exposure.txt
done

# create gsmr outcome file
> "${outFile}".outcome.txt
echo "${trait}" "${outFile}"."${trait}".prep >> "${outFile}".outcome.txt

# run mendelian randomisation (https://cnsgenomics.com/software/gcta/#Mendelianrandomisation)
(if [[ "${gsmr2beta}" == 1 ]]; then

	# run with HEIDI outlier procedure
	gcta64 \
	--mbfile "${outFile}".ref.data.txt \
	--keep "${outFile}".sample.txt \
	--gsmr-file "${outFile}".exposure.txt "${outFile}".outcome.txt \
	--gsmr2-beta \
	--gsmr-direction 2 \
	--effect-plot \
	--gwas-thresh 5e-8 \
	--heidi-thresh 0.01 \
	--clump-r2 ${clumpr2} \
	--clump-kb ${clumpkb} \
	--diff-freq 0.2 \
	--gsmr-snp-min 10 \
	--out "${outFile}" &

	# run without HEIDI outlier procedure
	gcta64 \
	--mbfile "${outFile}".ref.data.txt \
	--keep "${outFile}".sample.txt \
	--gsmr-file "${outFile}".exposure.txt "${outFile}".outcome.txt \
	--gsmr2-beta \
	--gsmr-direction 2 \
	--effect-plot \
	--gwas-thresh 5e-8 \
	--heidi-thresh 0 \
	--clump-r2 ${clumpr2} \
	--clump-kb ${clumpkb} \
	--diff-freq 0.2 \
	--gsmr-snp-min 10 \
	--out "${outFile}".noHEIDI

else

	# run with HEIDI outlier procedure
	gcta64 \
	--mbfile "${outFile}".ref.data.txt \
	--keep "${outFile}".sample.txt \
	--gsmr-file "${outFile}".exposure.txt "${outFile}".outcome.txt \
	--gsmr-direction 2 \
	--effect-plot \
	--gwas-thresh 5e-8 \
	--heidi-thresh 0.01 \
	--clump-r2 ${clumpr2} \
	--clump-kb ${clumpkb} \
	--diff-freq 0.2 \
	--gsmr-snp-min 10 \
	--out "${outFile}" & 

	# run without HEIDI outlier procedure
	gcta64 \
	--mbfile "${outFile}".ref.data.txt \
	--keep "${outFile}".sample.txt \
	--gsmr-file "${outFile}".exposure.txt "${outFile}".outcome.txt \
	--gsmr-direction 2 \
	--effect-plot \
	--gwas-thresh 5e-8 \
	--heidi-thresh 0 \
	--clump-r2 ${clumpr2} \
	--clump-kb ${clumpkb} \
	--diff-freq 0.2 \
	--gsmr-snp-min 10 \
	--out "${outFile}".noHEIDI

fi
wait)

# clean up
rm -f "${outFile}"."${trait}".prep
rm -f "${outFile}".sample.txt
rm -f "${outFile}".exposure.txt
rm -f "${outFile}".outcome.txt
rm -f "${outFile}".ref.data.txt
chmod -R 770 "${outFile}"*
echo $'\n'"--- Completed: GSMR analysis ---"

