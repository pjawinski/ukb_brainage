#!/bin/bash

# =====================================
# === Add variables to metal output ===
# =====================================

# get arguments
targetDir="${1}" # set target directory | targetDir="results/pwr_cz_alpha/metal/"
metalInput="${2}" # list of input sumstats used in meta-analysis | metalInput=$(echo $(ls $targetDir/metal.input*) | sed 's/ /,/g')
metalOutput="${3}" # metal output file | metalOutput="${targetDir}/metal.ivweight.1.txt"
idCol="${4}" # column name in input file that contains variant identifiers | idCol="ID"
carryOver="${5}" # variables to be added | carryOver="CHR,BP,KGP_REF,KGP_ALT,KGP_EUR_AF,KGP_RSID"

# echo settings
echo $'\n'"--- Carry over variables from metal input to metal output ---"
echo "metalInput: ${metalInput}"
echo "metalOutput: ${metalOutput}"
echo "idCol: ${idCol}"
echo "carryOver: ${carryOver}"$'\n'

# set targetdir and make folder
targetDir="$(readlink -f "${targetDir}")"

# convert list of metal input sumstats from comma- to space-delimited
metalInput=$(echo "${metalInput}" | sed 's/,/ /g')
metalOutput=$(echo "${metalOutput}" | sed 's/,/ /g')

# get variables from metal input sumstats
echo "Getting information from metal input sumstats."
> "${targetDir}"/metal.carryOver.txt
for file in ${metalInput}; do
    echo " - getting information from ${file}"
    zcat "${file}" | mawk -v cols="${idCol},${carryOver}" '
        BEGIN { ncols=split(cols,colnames,",") }
        NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } } }
        { for(i=1;i<=ncols;++i) { if(i==1) { output = $colidx[i] } else { output = output" "$colidx[i] } }; print output }
        ' >> "${targetDir}"/metal.carryOver.txt
done

# Identify unique rows
echo "Identifying unique lines."
mawk -v targetDir="${targetDir}" '
    !seen[$0]++ { print > targetDir"/metal.carryOver.tmp.txt"; nuniq++ }
    (FNR) % 5000000 == 0 { printf(" - lines processed: %d | unique lines: %d\n", FNR, nuniq) }
    END { printf(" - lines processed: %d | unique lines: %d\n", FNR-1, nuniq) }
    ' "${targetDir}"/metal.carryOver.txt
\mv "${targetDir}"/metal.carryOver.tmp.txt "${targetDir}"/metal.carryOver.txt

# Add variables to metal output
echo "Carrying over variables from metal input to metal output."
(for file in ${metalOutput}; do (
    header=$(echo $(head -1 "${file}") $(head -1 "${targetDir}"/metal.carryOver.txt | awk '{ sub($1" ","",$0); print $0 }'))
    mawk -v header="${header}" '
        BEGIN { print header }
        NR==FNR { for(i=2;i<=NF;++i) { if (i==2) { output = $i } else { output = output" "$i } }; snp[$1]=output; next }
        FNR==1 { next }
        { for(i=1;i<=NF;++i) { if (i==1) { output = $i } else { output = output" "$i } }; print output, snp[$1] }
        ' "${targetDir}"/metal.carryOver.txt "${file}" > "${file}".tmp
    \mv "${file}".tmp "${file}"

    # rename and reorder columns
    echo "Renaming and reordering columns of metal output."
    inputFile="${file}"
    outputFile=$(echo "${inputFile}" | sed 's/.1.txt//g').gz
    scriptDir=$(dirname "$0")
    Rscript "${scriptDir}"/metal.carryOver.R "${inputFile}" "${outputFile}" ) &
done
wait)

# clean up
rm -f "${targetDir}"/metal.carryOver.txt

