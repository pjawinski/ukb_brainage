#!/bin/bash

# ================================================
# === add 1KGP RSID, CHR, BP, REF, ALT, EUR_AF ===
# ================================================

# get arguments
sumstats=$1 # sumstats="data/sumstats/harmonized/03_ad_wrightman_2021.tmp.txt"
kgp=$2 # kgp="data/1kgp/v5a/ALL.chr\${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.meta.vcf.gz"

# echo settings
echo $'\n'"--- Adding 1KGP RSID, CHR, BP, REF, ALT, EUR_AF ---"
echo "sumstats: ${sumstats}"
echo "1kgp: ${kgp}"$'\n'

# use 1kgp to determine REF, ALT, RSID, and EUR_AF
echo "Matching variants with 1kgp entries."
(
for chr in {1..22}; do (
mawk -v sumstats="${sumstats}" -v chr="${chr}" '

    # Skip header line
    NR==1 { next }

    # Cohort file: save content in array
    NR==FNR && $1 == chr {snp[$1":"$2":"$3":"$4]=$0; nvariants++; next} 
    NR==FNR { next }

    # 1kgp file: skip lines starting with hashtag
    /^[#]/ { next } 

    # get EUR_AF from 8th column
    { sub(/.*EUR_AF=/,"",$8); sub(/[;].*/,"",$8); 
      nAF=split($8,EUR_AF,",") }

    # loop over alternate allele[s
    { nALT=split($5,ALT,","); for(i=1;i<=nALT;++i) {

        # find match by CHR:BP:REF:ALT or CHR:BP:REF:ALT
        if ($1":"$2":"$4":"ALT[i] in snp) { 
            print snp[$1":"$2":"$4":"ALT[i]],$4,ALT[i],EUR_AF[i],$3 > sumstats".chr"chr".txt";
            delete snp[$1":"$2":"$4":"ALT[i]]; matchCount++ }

        else if ($1":"$2":"ALT[i]":"$4 in snp) {
            print snp[$1":"$2":"ALT[i]":"$4],$4,ALT[i],EUR_AF[i],$3 > sumstats".chr"chr".txt";
            delete snp[$1":"$2":"ALT[i]":"$4]; matchCount++ }
        }
    }

    # print status
    END { printf(" - chr: %d | variants processed: %d | 1kgp matched: %d\n", chr, nvariants, matchCount) }
    ' OFS='\t' "${sumstats}" <(eval zcat "${kgp}")
) &
done
wait
)

# save matched variants in a single file
echo "Saving matched variants in a single file."
echo CHR$'\t'BP$'\t'A1$'\t'A2$'\t'KGP_REF$'\t'KGP_ALT$'\t'KGP_EUR_AF$'\t'KGP_RSID > "${sumstats}"
cat "${sumstats}".chr{1..22}.txt >> "${sumstats}"

# clean up
echo "Cleaning up target directory."
rm -f "${sumstats}".chr*
echo "--- completed: Adding 1KGP RSID, CHR, BP, REF, ALT, EUR_AF ---"

