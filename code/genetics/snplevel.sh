#!/bin/bash

# ==============================================
# === combine snplevel results across traits ===
# ==============================================

# get arguments
traitlist="${1}" # traitlist="gap_gm gap_wm gap_gwm"
conditionalFiles="${2}" # conditionalFiles="results/gap_gm/conditional/conditional.cleaned.tophits.annovar.txt results/gap_wm/conditional/conditional.cleaned.tophits.annovar.txt results/gap_gwm/conditional/conditional.cleaned.tophits.annovar.txt"
nonsynFiles="${3}" # nonsynFiles="results/gap_gm/conditional/conditional.cleaned.tophits.annovar.nonsynonymous.summary.txt results/gap_wm/conditional/conditional.cleaned.tophits.annovar.nonsynonymous.summary.txt results/gap_gwm/conditional/conditional.cleaned.tophits.annovar.nonsynonymous.summary.txt"
catalogFiles="${4}" # catalogFiles="results/gap_gm/catalog/catalog.by.locus.txt results/gap_wm/catalog/catalog.by.locus.txt results/gap_gwm/catalog/catalog.by.locus.txt"
finemapFiles="${5}" # finemapFiles="results/gap_gm/credibleSet/credibleSet.ls.txt results/gap_wm/credibleSet/credibleSet.ls.txt results/gap_gwm/credibleSet/credibleSet.ls.txt"
smrFiles="${6}" # smrFiles="results/gap_gm/smr/smr.summary.txt results/gap_wm/smr/smr.summary.txt results/gap_gwm/smr/smr.summary.txt"
eqtlFiles="${7}" # eqtlFiles="results/gap_gm/eqtl/eqtl.summary.txt results/gap_wm/eqtl/eqtl.summary.txt results/gap_gwm/eqtl/eqtl.summary.txt"
popsFiles="${8}" # popsFiles="results/gap_gm/pops/pops.prio.txt results/gap_wm/pops/pops.prio.txt results/gap_gwm/pops/pops.prio.txt"
targetDir="${9}" # targetDir="results/combined"
geneticsDir="${10}" # geneticsDir="data/genetics"
LDsample="${11}" # LDsample="data/gap_gm/gap_gm.txt"
pthresh=${12} # pthresh=1E-6

# echo settings
echo $'\n'"--- Settings for combining snplevel results across traits ---"
echo "traitlist: "${traitlist}
echo "conditionalFiles: "${conditionalFiles}
echo "nonsynFiles: "${nonsynFiles}
echo "catalogFiles: "${catalogFiles}
echo "finemapFiles: "${finemapFiles}
echo "smrFiles: "${smrFiles}
echo "eqtlFiles: "${eqtlFiles}
echo "popsFiles: "${popsFiles}
echo "targetDir: "${targetDir}
echo "geneticsDir: "${geneticsDir}
echo "LDsample: "${LDsample}
echo "pthresh: "${pthresh}

# make targetDir
rm -f "${targetDir}/snplevel.clumping"
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# combine conditional results
traitlist=($traitlist)
conditionalFiles=($conditionalFiles)
awk -F'\t' 'NR == 1 { output="TRAIT\t"$3"\t"$4"\t"$5; for(i=11;i<=NF;++i) { output=output"\t"$i}; print output; next}' OFS='\t' "${conditionalFiles[0]}" > "${targetDir}/snplevel.conditional.txt"
for (( i=0; i<${#traitlist[@]}; i++ )); do
	awk -F'\t' -v trait="${traitlist[$i]}" -v pthresh="${pthresh}" 'FNR == 1 { next }
	$2 == 1 && $17 < pthresh && $31 < pthresh { output=trait"\t"$3"\t"$4"\t"$5; for(i=11;i<=NF;++i) { output=output"\t"$i}; print output}
	' OFS='\t' "${conditionalFiles[$i]}" >> ${targetDir}/snplevel.conditional.txt
done

# sort by ID and p-value and remove redundant SNPs
(head -n 1 "${targetDir}/snplevel.conditional.txt" && tail -n +2 "${targetDir}/snplevel.conditional.txt" | sort -k4,4 -k11,11g | awk -F'\t' '!seen[$4]++') > "${targetDir}/snplevel.sumstats.txt"

# get CHR, ID and P values of all variations from summary statistics
cat "${targetDir}/snplevel.sumstats.txt" | awk '{ print $2, $4, $11}' OFS='\t' > "${targetDir}/snplevel.sumstats4clumping.txt"

# get participant IDs of LD reference sample
awk '{ print $1, $2 }' "${LDsample}" > "${targetDir}/snplevel.subs.txt"

# get relevant chromosomes
chr=$(cat ${targetDir}/snplevel.sumstats.txt | awk 'NR==1 { next } { print $2 }' | sort -u ) 
chr=($chr)

# do clumping
mkdir -p "${targetDir}/snplevel.clumping"
for i in ${chr[@]}; do (
	awk -v chr="${i}" 'NR==1 { print; next } $1==chr { print }' OFS='\t' "${targetDir}/snplevel.sumstats4clumping.txt" > "${targetDir}/snplevel.clumping/sumstats4clumping_${i}.txt"
	LD_src="${geneticsDir}/chr$i/imp_mri_qc/bed/chr${i}_mri_qc"
	plink --bfile "$LD_src" \
	--keep "${targetDir}/snplevel.subs.txt" \
	--clump "${targetDir}/snplevel.clumping/sumstats4clumping_${i}.txt" \
	--clump-snp-field "ID" \
	--clump-field "P" \
	--clump-p1 1.000 \
	--clump-p2 1.000 \
	--clump-r2 0.1 \
	--clump-kb 10000 \
	--clump-verbose \
	--out "$targetDir/snplevel.clumping/chr${i}"
	rm -f "${targetDir}/snplevel.clumping/sumstats4clumping_chr${i}.txt"
	) &
done
wait

# Merge clumped files
cat ${targetDir}/snplevel.clumping/*.clumped > ${targetDir}/snplevel.clumping/clumped.tmp
awk '$11 ~ /^[0-9]+$/ { indexsnp = $3; print $3, $3, 0, 1, "A1/A2"; next } $1=="" || $1=="KB" || $1=="CHR" || $1=="(INDEX)" || $2 == "not" || $1 == "RANGE:" || $1 == "SPAN:" || $1 ~ /-/ { next } {print indexsnp, $1, $2, $3, $4}' OFS='\t' ${targetDir}/snplevel.clumping/clumped.tmp > ${targetDir}/snplevel.clumping/clumped 
rm -f ${targetDir}/snplevel.clumping/clumped.tmp

# add lead association (across phenotypes) to conditional analysis results
awk -F'\t' 'NR==FNR { snp[$2]=$1"\t"$3"\t"$4"\t"$5; next} FNR==1 { next } $4 in snp { print $0, snp[$4] }' OFS="\t" "${targetDir}/snplevel.clumping/clumped" "${targetDir}/snplevel.conditional.txt" > "${targetDir}/snplevel.txt"

# replace "A1/A2" (if SNP is Lead SNP) by actual alleles 
awk -F'\t' '$NF == "A1/A2" { output=$1; for(i=2;i<NF;++i) { output=output"\t"$i }; output=output"\t"$5"/"$6; print output; next} { print }' OFS="\t" "${targetDir}/snplevel.txt" > "${targetDir}/snplevel.txt.tmp"
\mv "${targetDir}/snplevel.txt.tmp" "${targetDir}/snplevel.txt"

# add CHR, BP, and GWS status of lead association
header=$(echo "$(head -1 "${targetDir}/snplevel.conditional.txt" | sed 's/LEAD_SNP_/COJO_/g')"$'\t'LEAD_SNP$'\t'LEAD_SNP_KB$'\t'LEAD_SNP_RSQ$'\t'LEAD_SNP_ALLELES$'\t'LEAD_SNP_CHR$'\t'LEAD_SNP_BP$'\t'LEAD_SNP_GWS)  
awk -F'\t' -v header="${header}" '
	BEGIN { print header }
	NR==FNR && $4==$27 && (!($27 in pval) || $25 < pval[$27]) { pval[$27]=$25; $25 < 5E-8 ? gws=1 : gws=2; chrbp[$27]=$2"\t"$3"\t"gws; next }
	NR==FNR { next }
	{ print $0, chrbp[$27] }' OFS="\t" "${targetDir}/snplevel.txt" "${targetDir}/snplevel.txt" > "${targetDir}/snplevel.txt.tmp"
\mv "${targetDir}/snplevel.txt.tmp" "${targetDir}/snplevel.txt"

# sort by LEAD_SNP_GWS, LEAD_SNP_CHR, LEAD_SNP_BP, and P
awk -F'\t' 'NR == 1 { print; next } $31=="X" { $31=23 } $31=="Y" { $31=24 } $31=="XY" { $31=25 } $31=="MT" { $31=26 } { print }' OFS="\t" "${targetDir}/snplevel.txt" | 
	sort -t$'\t' -k33,33g -k31,31g -k32,32g -k11,11g > "${targetDir}/snplevel.txt.tmp"
\mv "${targetDir}/snplevel.txt.tmp" "${targetDir}/snplevel.txt"

# add LOCUS_COUNT and DISCOVERY_COUNT
header=$(echo LOCUS_COUNT$'\t'DISCOV_COUNT$'\t'"$(head -1 "${targetDir}/snplevel.txt")") 
awk -F'\t' -v header="$header" '
	BEGIN { print header; leadsnp="null"; locus_count = 0; snp_count = 1 }
	NR==1 { next } leadsnp != $27 { leadsnp = $27; locus_count++; snp_count = 1; print locus_count, snp_count, $0; next }
	{ snp_count++; print locus_count, snp_count, $0; next }' OFS='\t' "${targetDir}/snplevel.txt" > "${targetDir}/snplevel.txt.tmp"
\mv "${targetDir}/snplevel.txt.tmp" "${targetDir}/snplevel.txt"

# add finemap results
finemapFiles=($finemapFiles)
header=$(echo "$(head -1 "${targetDir}/snplevel.txt")"$'\t'FINEMAP_regionalh2$'\t'FINEMAP_snpCausal$'\t'FINEMAP_bestK$'\t'FINEMAP_credibleSet_size$'\t'FINEMAP_genes) 
rm -f "${targetDir}/snplevel.txt.tmp"
for (( i=0; i<${#finemapFiles[@]}; i++ )); do
	awk -F'\t' -v header="$header" -v trait="${traitlist[$i]}" ' BEGIN { print header } NR == FNR { sub(" .*","",$4); finemap[$1]=$4"\t"$9"\t"$10"\t"$11"\t"$13; next } FNR==1 { next } 
	$6 in finemap && $3 == trait { print $0, finemap[$6]; next } { print $0 }
	' OFS='\t' "${finemapFiles[$i]}" "${targetDir}/snplevel.txt" > "${targetDir}/snplevel.txt.tmp"
	\mv "${targetDir}/snplevel.txt.tmp" "${targetDir}/snplevel.txt"
done
awk -F'\t' 'NR==1 { fieldnum=NF; print; next } NF < fieldnum { print $0, "NA", "NA", "NA", "NA", "NA"; next} { print }' OFS='\t' "${targetDir}/snplevel.txt" > "${targetDir}/snplevel.txt.tmp"
\mv "${targetDir}/snplevel.txt.tmp" "${targetDir}/snplevel.txt"

# add nonsynonymous exonic variations
nonsynFiles=($nonsynFiles)
header=$(echo "$(head -1 "${targetDir}/snplevel.txt")"$'\t'nonsyn_customPthresh$'\t'nonsyn_gwsPthresh) 
rm -f "${targetDir}/snplevel.txt.tmp"
for (( i=0; i<${#nonsynFiles[@]}; i++ )); do
	awk -F'\t' -v header="$header" -v trait="${traitlist[$i]}" ' BEGIN { print header } NR == FNR { nonsyn[$1]=$2"\t"$3; next } FNR==1 { next } 
	$6 in nonsyn && $3 == trait { print $0, nonsyn[$6]; next } { print $0 }
	' OFS='\t' "${nonsynFiles[$i]}" "${targetDir}/snplevel.txt" > "${targetDir}/snplevel.txt.tmp"
	\mv "${targetDir}/snplevel.txt.tmp" "${targetDir}/snplevel.txt"
done
awk -F'\t' 'NR==1 { fieldnum=NF; print; next } NF < fieldnum { print $0, "NA", "NA"; next} { print }' OFS='\t' "${targetDir}/snplevel.txt" > "${targetDir}/snplevel.txt.tmp"
\mv "${targetDir}/snplevel.txt.tmp" "${targetDir}/snplevel.txt"
	
# add SMR results
smrFiles=($smrFiles)
header=$(echo "$(head -1 "${targetDir}/snplevel.txt")"$'\t'SMR_eqtl$'\t'SMR_sqtl) 
rm -f "${targetDir}/snplevel.txt.tmp"
for (( i=0; i<${#smrFiles[@]}; i++ )); do
	awk -F'\t' -v header="$header" -v trait="${traitlist[$i]}" ' BEGIN { print header } NR == FNR { smr[$2]=$3"\t"$4; next } FNR==1 { next } 
	$6 in smr && $3 == trait { print $0, smr[$6]; next } { print $0 }
	' OFS='\t' "${smrFiles[$i]}" "${targetDir}/snplevel.txt" > "${targetDir}/snplevel.txt.tmp"
	\mv "${targetDir}/snplevel.txt.tmp" "${targetDir}/snplevel.txt"
done
awk -F'\t' 'NR==1 { fieldnum=NF; print; next } NF < fieldnum { print $0, "NA", "NA"; next} { print }' OFS='\t' "${targetDir}/snplevel.txt" > "${targetDir}/snplevel.txt.tmp"
\mv "${targetDir}/snplevel.txt.tmp" "${targetDir}/snplevel.txt"

# add GTEx eQTL matches
eqtlFiles=($eqtlFiles)
header=$(echo "$(head -1 "${targetDir}/snplevel.txt")"$'\t'GTEx_singleTissue$'\t'GTEx_multiTissue) 
rm -f "${targetDir}/snplevel.txt.tmp"
for (( i=0; i<${#eqtlFiles[@]}; i++ )); do
	awk -F'\t' -v header="$header" -v trait="${traitlist[$i]}" ' BEGIN { print header } NR == FNR { eqtl[$2]=$3"\t"$4; next } FNR==1 { next } 
	$6 in eqtl && $3 == trait { print $0, eqtl[$6]; next } { print $0 }
	' OFS='\t' "${eqtlFiles[$i]}" "${targetDir}/snplevel.txt" > "${targetDir}/snplevel.txt.tmp"
	\mv "${targetDir}/snplevel.txt.tmp" "${targetDir}/snplevel.txt"
done
awk -F'\t' 'NR==1 { fieldnum=NF; print; next } NF < fieldnum { print $0, "NA", "NA"; next} { print }' OFS='\t' "${targetDir}/snplevel.txt" > "${targetDir}/snplevel.txt.tmp"
\mv "${targetDir}/snplevel.txt.tmp" "${targetDir}/snplevel.txt"

# add polygenic priority scores (PoPS)
popsFiles=($popsFiles)
header=$(echo "$(head -1 "${targetDir}/snplevel.txt")"$'\t'PoPS) 
rm -f "${targetDir}/snplevel.txt.tmp"
for (( i=0; i<${#popsFiles[@]}; i++ )); do
	awk -F'\t' -v header="$header" -v trait="${traitlist[$i]}" ' BEGIN { print header } NR == FNR { pops[$1]=$4; next } FNR==1 { next } 
	$6 in pops && $3 == trait { print $0, pops[$6]; next } { print $0 }
	' OFS='\t' "${popsFiles[$i]}" "${targetDir}/snplevel.txt" > "${targetDir}/snplevel.txt.tmp"
	\mv "${targetDir}/snplevel.txt.tmp" "${targetDir}/snplevel.txt"
done
awk -F'\t' 'NR==1 { fieldnum=NF; print; next } NF < fieldnum { print $0, "NA"; next} { print }' OFS='\t' "${targetDir}/snplevel.txt" > "${targetDir}/snplevel.txt.tmp"
\mv "${targetDir}/snplevel.txt.tmp" "${targetDir}/snplevel.txt"

# add GWAS catalog matches
catalogFiles=($catalogFiles)
header=$(echo "$(head -1 "${targetDir}/snplevel.txt")"$'\t'GWAS_CATALOG) 
rm -f "${targetDir}/snplevel.txt.tmp"
for (( i=0; i<${#catalogFiles[@]}; i++ )); do
	awk -F'\t' -v header="$header" -v trait="${traitlist[$i]}" ' BEGIN { print header } NR == FNR { catalog[$3]=$4; next } FNR==1 { next } 
	$6 in catalog && $3 == trait { print $0, catalog[$6]; next } { print $0 }
	' OFS='\t' "${catalogFiles[$i]}" "${targetDir}/snplevel.txt" > "${targetDir}/snplevel.txt.tmp"
	\mv "${targetDir}/snplevel.txt.tmp" "${targetDir}/snplevel.txt"
done
awk -F'\t' 'NR==1 { fieldnum=NF; print; next } NF < fieldnum { print $0, "NA"; next} { print }' OFS='\t' "${targetDir}/snplevel.txt" > "${targetDir}/snplevel.txt.tmp"
\mv "${targetDir}/snplevel.txt.tmp" "${targetDir}/snplevel.txt"

# create seperate file with gwas hits only
awk -F'\t' 'NR == 1 { print; for(i=1; i<=NF; i++) { if($i=="COJO_pJ") pJ=i; if($i=="P") P=i } }
	$P < 5E-8 && $pJ < 5E-8 { print }' OFS='\t' "${targetDir}/snplevel.txt" > "${targetDir}/snplevel.gws.txt"

# clean up
echo "Cleaning up."
rm -f ${targetDir}/snplevel.conditional.txt
rm -f ${targetDir}/snplevel.sumstats.txt
rm -f ${targetDir}/snplevel.subs.txt
rm -f ${targetDir}/snplevel.sumstats4clumping.txt
rm -f ${targetDir}/snplevel.clumping/sumstats4clumping_*
tar -cvzf ${targetDir}/snplevel.clumping.tar.gz --directory=${targetDir} snplevel.clumping --remove-files
chmod -R 770 ${targetDir}
echo "--- Combining snp-level results completed. --- "


