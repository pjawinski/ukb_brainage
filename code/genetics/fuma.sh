#!/bin/bash

# ==============================================
# === Preparing summary statistics for FUMA ===
# ==============================================

# get arguments
trait="${1}" # trait="gap_gm"
targetDir="${2}" # targetDir="results/gap_gm/gwama/eur/fuma"
sumstats="${3}" # sumstats="results/gap_gm/gwama/eur/metal.ivweight.qc.gz"
sumstatsCols="${4}" # sumstatsCols="CHR,BP,ID,A1,A2,BETA,SE,P,N" | first column must be chromosome column due to filtering for Y, XY, and MT
sumstatsColNames="${5}" # sumstatsColNames="CHR,BP,SNP,A1,A2,BETA,SE,P,N" | column names must be in accordance with fuma input https://fuma.ctglab.nl/tutorial#prepare-input-files

# echo settings
echo $'\n'"--- Preparing summary statistics for FUMA ---"
echo "trait: ${trait}"
echo "targetDir: ${targetDir}"
echo "sumstats: ${sumstats}"
echo "sumstatsCols: ${sumstatsCols}"
echo "sumstatsColNames: ${sumstatsColNames}"$'\n'

# set targetdir and make folder
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# prepare sumstats
echo "Preparing sumstats."
header=$(echo "${sumstatsColNames}" | sed 's%,%\t%g')
awk -v header="${header}" -v cols="${sumstatsCols}" '
	BEGIN { ncols=split(cols,colnames,","); print header }
    NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } }; next}
    NR>1 { for(i=1;i<=ncols;++i) { i==1 ? output=$colidx[i] : output=output"\t"$colidx[i] } }
    $colidx[1]=="X" || ($colidx[1]!="Y" && $colidx[1]!="XY" && $colidx[1]!="MT" && $colidx[1] < 24) { print output }
	' OFS="\t" <(gzip -dc "${sumstats}") > "${targetDir}"/fuma."${trait}"

# clean up
echo "Cleaning up."
pigz -f "${targetDir}"/fuma."${trait}"
chmod -R 770  "${targetDir}"/fuma."${trait}".gz
echo "--- FUMA preparation finished --- "
