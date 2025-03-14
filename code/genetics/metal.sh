#!/bin/bash

# =================
# === run METAL ===
# =================

# get arguments
trait="${1}" # trait="gap_gm"
targetDir="${2}" # targetDir="results/${trait}/prepare"
cohorts="${3}" # cohorts="AFR,AMR,CSA,EAS,EUR,LIFE,MID" # set output file handler & make sure that the order matches the list of sumstats (below)
sumstats="${4}" # sumstats=$(echo $(for i in $(echo $cohorts | sed 's/,/ /g'); do echo ${targetDir}/${i}/sumstats.txt.gz; done) | sed 's/ /,/g')
type="${5}" # type="ivweight" | type="nweight"
idCol="${6}" # idCol="ID"
carryOver="${7}" # carryOver="CHR,BP"
leaveOneOut="${8}" # leaveOneOut=0

# echo settings
echo $'\n'"--- Running METAL ---"
echo "trait: ${trait}"
echo "targetDir: ${targetDir}"
echo "cohorts: ${cohorts}"
echo "type: ${type}"
echo "idCol: ${idCol}"
echo "carryOver: ${carryOver}"
echo "leaveOneOut ${leaveOneOut}"$'\n'

# get space-delimited list of cohorts and sumstats
cohorts=$(echo "${cohorts}" | sed 's/,/ /g'); cohorts=($cohorts)
sumstats=$(echo "${sumstats}" | sed 's/,/ /g'); sumstats=($sumstats)

# ===== define function for getting metal command ====
createMETAL(){
# set type of analysis
if [ "${type}" == "ivweight" ]; then
	scheme="STDERR"
elif [ "${type}" == "nweight" ]; then
	scheme="SAMPLESIZE"
fi

# main part
read -r -d '' metalCMD <<- eof
# scheme
SCHEME ${scheme}
AVERAGEFREQ ON
MINMAXFREQ ON
EFFECT_PRINT_PRECISION 8
STDERR_PRINT_PRECISION 8

# custom variables
CUSTOMVARIABLE Ntotal
LABEL Ntotal as N

# define header
MARKER ID
ALLELE A1 A2
FREQ A1_FREQ
EFFECT BETA
STDERR SE
PVALUE P
WEIGHT N

# process files
eof

# add files to be processed
for i in ${inputFiles}; do 
metalCMD=$(echo "${metalCMD}"$'\n'PROCESS "${i}")
done

# add output file and analyze commands
metalCMD=$(echo "${metalCMD}"$'\n'$'\n'\# run analysis)
if [ "${type}" == "ivweight" ]; then
	metalCMD=$(echo "${metalCMD}"$'\n'OUTFILE "${targetDir}"/metal.ivweight."${leaveOut}" .txt)
elif [ "${type}" == "nweight" ]; then
	metalCMD=$(echo "${metalCMD}"$'\n'OUTFILE "${targetDir}"/metal.nweight."${leaveOut}" .txt)
fi
metalCMD=$(echo "${metalCMD}"$'\n'ANALYZE HETEROGENEITY$'\n'QUIT)
}
# ===== end function for getting metal command ====
	
# get metal command to run analysis over all samples
leaveOut=""
inputFiles=${sumstats[@]}
createMETAL

# run metal without/with leave-one-out method
if [[ ${leaveOneOut} = 0 ]]; then
	echo "Run METAL over all samples."
	metal <(echo "$metalCMD")
elif [[ ${leaveOneOut} = 1 ]]; then
	echo "Run METAL with leave-one-out method."
	(echo " - over all samples."
	metal <(echo "$metalCMD") &
	for (( j=0; j<${#cohorts[@]}; j++ )); do (
		echo " - w/out ${cohorts[j]}"
		leaveOut="loo.${cohorts[j]}."
		inputFiles=${sumstats[@]}; inputFiles=($inputFiles); unset 'inputFiles[j]'; inputFiles=${inputFiles[@]}
		createMETAL
		metal <(echo "${metalCMD}") ) &
	done
	wait)
fi

# add variables from metal input sumstats
metalInput=$(echo ${sumstats[@]} | sed 's/ /,/g')
if [[ ${leaveOneOut} = 0 ]]; then
	metalOutput="${targetDir}/metal.${type}.1.txt"
elif [[ ${leaveOneOut} = 1 ]]; then
	metalOutput=$(echo "${targetDir}"/metal."${type}".1.txt $(for i in ${cohorts[@]}; do echo "${targetDir}"/metal."${type}".loo."${i}".1.txt; done) | sed 's/ /,/g')
fi
scriptDir=$(dirname "$0")
"${scriptDir}"/metal.carryOver.sh "${targetDir}" "${metalInput}" "${metalOutput}" "${idCol}" "${carryOver}"

# clean up
echo "Cleaning up target directory."
metalOutput=$(echo "${metalOutput}" | sed 's/,/ /g')
rm -f ${metalOutput}
for file in ${metalOutput}; do
    newFileName=$(echo "${file}" | sed 's/.1.txt/.info/g')
    mv "${file}".info "${newFileName}"
    chmod 770 "${newFileName}" 
done
echo "--- Completed: Running Metal ---"
