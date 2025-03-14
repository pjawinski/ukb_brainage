#!/bin/bash

# =============================================
# === run phenome-wide association analysis ===
# =============================================

# get arguments
trait="$1" # trait="gap_gm"
traitFile="$2" # traitFile="data/${trait}/${trait}.txt"
covsFile="$3" # covsFile="data/${trait}/covs.txt"
covs="$4" # covs="sex,age,age2,ac1,ac2,TIV"
targetDir="$5" # targetDir="results/${trait}/phesant"
phenoFile="$(readlink -f "$6")" # phenoFile="data/basket/20200409_2007685/data/ukb41573_imaging_wba_phesant.csv"
phesantDir="$7" # phesantDir="/fast/software/PHESANT"
nparts=$8 # nparts=1
standardise="$9" # standardise="FALSE"

# echo settings
echo $'\n'"--- PHESANT analysis Settings ---"
echo "trait: ${trait}"
echo "traitFile: ${traitFile}"
echo "covsFile: ${covsFile}"
echo "covs: ${covs}"
echo "targetDir: ${targetDir}"
echo "phenoFile: ${phenoFile}"
echo "phesantDir: ${phesantDir}"
echo "nparts: ${nparts}"
echo "standardise: ${standardise}"$'\n'

# set targetdir and make folder
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# get exposure file
echo "Creating exposure file."
awk -F '\t' -v trait="${trait}" 'BEGIN { print "eid,"trait } NR > 1 { print $1, $3}' OFS=',' "${traitFile}" > "${targetDir}"/phesant.exposure.csv

# get confounder file
echo "Creating confounder file."
header="eid,${covs}" 
awk -F '\t' -v header="${header}" -v covs="${covs}" 'BEGIN { print header; ncovs=split(covs,covnames,",") }
	NR==1 { for(i=1;i<=ncovs;++i) { for(j=2;j<=NF;++j) { if($j==covnames[i]) { colnums[i]=j } } } }
	NR>1 { output=$1; for(i=1;i<=ncovs;++i) { output=output","$colnums[i] }; print output}' "${covsFile}" > "${targetDir}"/phesant.confounder.csv

# run analysis
echo "Running phenome-wide scan."$'\n'
mkdir -p "${targetDir}"/phesant.output
cd "${phesantDir}"/WAS

for ((i = 1; i <= nparts; i++)); do (
Rscript phenomeScan.r \
--phenofile="${phenoFile}" \
--traitofinterestfile="${targetDir}"/phesant.exposure.csv \
--variablelistfile="${phesantDir}"/variable-info/outcome-info.tsv \
--datacodingfile="${phesantDir}"/variable-info/data-coding-ordinal-info.txt \
--traitofinterest="${trait}" \
--resDir="${targetDir}"/phesant.output/ \
--userId="eid" \
--partIdx="${i}" \
--numParts="${nparts}" \
--confounderfile="${targetDir}"/phesant.confounder.csv \
--standardise="${standardise}" \
--genetic="FALSE" \
) &
done
wait
echo "Finished."

# post phenome scan results processing
cd "${phesantDir}"/resultsProcessing/
Rscript mainCombineResults.r \
--resDir="${targetDir}"/phesant.output/ \
--variablelistfile="${phesantDir}"/variable-info/outcome-info.tsv \
--numParts="${nparts}"

# create visualization
# mkdir -p "${targetDir}/results/web"
# cd ${phesantDir}/PHESANT-viz/bin/
# java -cp .:../jar/json-simple-1.1\ 2.jar ResultsToJSON "${targetDir}/results/results-combined.txt" "../node-positions.csv" "${targetDir}/results/web/java-json.json"

# set permissions
chmod -R 770 "${targetDir}"

