#!/bin/bash

# ====================================
# === Genetic correlation analysis ===
# ====================================

# get arguments
targetDir="$1" # targetDir="results/PWR_COGA_vs_LIFEadult/rgxy"
sumstatsX="$2" # sumstatsX=$(find results/LIFEadult/PWR*/*/*sumstats.gz)
sumstatsY="$3" # sumstatsY=$(find results/COGA/PWR*/*/*sumstats.gz)
xprefix="$4" # xprefix="LIFE_"
yprefix="$5" # yprefix="COGA_"
LDchr="$6" # LDchr="data/ldsc/resources/eur_w_ld_chr/"
threads=$7 # threads=50

# echo settings
echo $'\n'"--- Genetic Correlation Analysis Settings ---"
echo "targetDir: ${targetDir}"
echo "sumstatsX: ${sumstatsX}"
echo "sumstatsY: ${sumstatsY}"
echo "xprefix: ${xprefix}"
echo "yprefix: ${yprefix}"
echo "LDchr: ${LDchr}"
echo "threads: ${threads}"$'\n'

# set targetdir and make folder
mkdir -p "${targetDir}"
# targetDir="$(readlink -f "${targetDir}")"

# run genetic correlation analysis
mkdir -p "${targetDir}"/rg.output/
sumstatsX=($sumstatsX)
sumstatsY=($sumstatsY)

# get list of pairs (only calculate upper triangle if x == y)
pairs=()
if [ "$(echo ${sumstatsX[@]})" == "$(echo ${sumstatsY[@]})" ] ; then
	for x in ${sumstatsX[@]}; do
		sumstatsY=( "${sumstatsY[@]/$x}" ); #sumstatsY=($sumstatsY)
		for y in ${sumstatsY[@]}; do
			pairs+=("$x $y")
		done
	done
else
	for x in ${sumstatsX[@]}; do
		for y in ${sumstatsY[@]}; do
			pairs+=("$x $y")
		done
	done
fi

# run genetic correlation analysis
echo "Running genetic correlations."
(
for pair in "${pairs[@]}"; do
	((i=i%threads)); ((i++==0)) && wait
		(
		pair=($pair)
		x=${pair[0]}
		y=${pair[1]}
		xname=$(echo "${x}" | sed 's/.*\///g' | sed 's/.sumstats.gz//g' | sed 's/ldsc.//g')
		yname=$(echo "${y}" | sed 's/.*\///g' | sed 's/.sumstats.gz//g' | sed 's/ldsc.//g')
		ldsc.py \
			--rg "${x}","${y}" \
			--ref-ld-chr "${LDchr}" \
			--w-ld-chr "${LDchr}" \
			--out "${targetDir}"/rg.output/rg."${xprefix}${xname}"_VS_"${yprefix}${yname}"
		) &
done
wait
)

# collect results
echo "Collecting results."
echo p1$'\t'p2$'\t'rg$'\t'se$'\t'z$'\t'p$'\t'h2_obs$'\t'h2_obs_se$'\t'h2_int$'\t'h2_int_se$'\t'gcov_int$'\t'gcov_int_se > "${targetDir}"/rg.results
logFiles=$(ls "${targetDir}"/rg.output/*.log)
for i in ${logFiles}; do
    summary="$(grep -m1 -A 2 "Summary of Genetic Correlation Results" "${i}" | awk 'NR==3 { print }' OFS='\t')"
    echo "${summary}" | awk -v xprefix="${xprefix}" -v yprefix="${yprefix}" '{sub(".*/","",$1); sub("ldsc.","",$1); sub(".sumstats.gz","",$1); sub(".*/","",$2); sub("ldsc.","",$2); sub(".sumstats.gz","",$2); $1=xprefix$1; $2=yprefix$2; print}' OFS='\t' >> "${targetDir}"/rg.results
done

# clean up
tar -cvzf "${targetDir}"/rg.output.tar.gz --directory="${targetDir}" rg.output --remove-files
chmod -R 770 "${targetDir}"


