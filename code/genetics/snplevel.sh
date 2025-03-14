#!/bin/bash

# ==============================================
# === combine snplevel results across traits ===
# ==============================================

# get arguments
traitlist="${1}" # traitlist="gap_gm gap_wm gap_gwm"
conditionalFiles="${2}" # conditionalFiles="results/gap_gm/conditional/conditional.cleaned.tophits.annovar.txt results/gap_wm/conditional/conditional.cleaned.tophits.annovar.txt results/gap_gwm/conditional/conditional.cleaned.tophits.annovar.txt"
nonsynFiles="${3}" # nonsynFiles="results/gap_gm/conditional/conditional.cleaned.tophits.annovar.nonsynonymous.summary.txt results/gap_wm/conditional/conditional.cleaned.tophits.annovar.nonsynonymous.summary.txt results/gap_gwm/conditional/conditional.cleaned.tophits.annovar.nonsynonymous.summary.txt"
catalogFiles="${4}" # catalogFiles="results/gap_gm/catalog/catalog.by.locus.txt results/gap_wm/catalog/catalog.by.locus.txt results/gap_gwm/catalog/catalog.by.locus.txt"
sbayesFiles="${5}" # sbayesFiles="results/gap_gm/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.joined.genes.txt results/gap_wm/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.joined.genes.txt results/gap_gwm/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.joined.genes.txt"
susieRFiles="${6}" # susieRFiles="results/gap_gm/credibleSet/credibleSet.ls.txt results/gap_wm/credibleSet/credibleSet.ls.txt results/gap_gwm/credibleSet/credibleSet.ls.txt"
finemapFiles="${7}" # finemapFiles="results/gap_gm/gwama/eur/finemap/finemap.ls.genes.txt results/gap_wm/gwama/eur/finemap/finemap.ls.genes.txt results/gap_gwm/gwama/eur/finemap/finemap.ls.genes.txt"
smrFiles="${8}" # smrFiles="results/gap_gm/smr/smr.summary.txt results/gap_wm/smr/smr.summary.txt results/gap_gwm/smr/smr.summary.txt"
eqtlFiles="${9}" # eqtlFiles="results/gap_gm/eqtl/eqtl.summary.txt results/gap_wm/eqtl/eqtl.summary.txt results/gap_gwm/eqtl/eqtl.summary.txt"
popsFiles="${10}" # popsFiles="results/gap_gm/pops/pops.prio.txt results/gap_wm/pops/pops.prio.txt results/gap_gwm/pops/pops.prio.txt"
targetDir="${11}" # targetDir="results/combined"
chrFilehandler="${12}" # geneticsDir="data/genetics"
LDsample="${13}" # LDsample="data/gap_gm/gap_gm.txt"
pthresh=${14} # pthresh=1E-6

# echo settings
echo $'\n'"--- Settings for combining snplevel results across traits ---"
echo "traitlist: ${traitlist}"
echo "conditionalFiles: ${conditionalFiles}"
echo "nonsynFiles: ${nonsynFiles}"
echo "catalogFiles: ${catalogFiles}"
echo "sbayesFiles: ${sbayesFiles}"
echo "susieRFiles: ${susieRFiles}"
echo "finemapFiles: ${finemapFiles}"
echo "smrFiles: ${smrFiles}"
echo "eqtlFiles: ${eqtlFiles}"
echo "popsFiles: ${popsFiles}"
echo "targetDir: ${targetDir}"
echo "chrFilehandler: ${chrFilehandler}"
echo "LDsample: ${LDsample}"
echo "pthresh: ${pthresh}"

# make targetDir
rm -f "${targetDir}"/snplevel.clumping
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# combine conditional results
traitlist=($traitlist)
conditionalFiles=($conditionalFiles)
awk -F'\t' 'NR == 1 { output="TRAIT\t"$3"\t"$4"\t"$5; for(i=11;i<=NF;++i) { output=output"\t"$i}; print output; next}' OFS='\t' "${conditionalFiles[0]}" > "${targetDir}"/snplevel.conditional.txt
for (( i=0; i<${#traitlist[@]}; i++ )); do
	awk -F'\t' -v trait="${traitlist[$i]}" -v pthresh="${pthresh}" -v cols="CHR,BP,ID,P,LEAD_SNP_pJ" '
	BEGIN { ncols=split(cols,colnames,",") }
	NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } }; next } 
	$2 == 1 && $colidx[4] < pthresh && $colidx[5] < pthresh { output=trait"\t"$colidx[1]"\t"$colidx[2]"\t"$colidx[3]; for(i=11;i<=NF;++i) { output=output"\t"$i}; print output}
	' OFS='\t' "${conditionalFiles[$i]}" >> "${targetDir}"/snplevel.conditional.txt
done

# sort by ID and p-value and remove redundant SNPs
(head -n 1 "${targetDir}"/snplevel.conditional.txt && tail -n +2 "${targetDir}"/snplevel.conditional.txt | sort -k4,4 -k10,10g | awk -F'\t' '!seen[$4]++') > "${targetDir}"/snplevel.sumstats.txt

# get CHR, ID and P values of all variations from summary statistics
awk '{ print $2, $4, $10}' "${targetDir}"/snplevel.sumstats.txt OFS='\t' > "${targetDir}"/snplevel.sumstats4clumping.txt

# get participant IDs of LD reference sample
awk '{ print $1, $2 }' "${LDsample}" > "${targetDir}"/snplevel.subs.txt

# get relevant chromosomes
chr=$(awk 'NR==1 { next } { print $2 }' "${targetDir}"/snplevel.sumstats.txt | sort -u ) 
chr=($chr)

# do clumping
mkdir -p "${targetDir}"/snplevel.clumping
(for i in ${chr[@]}; do (
	awk -v chr="${i}" 'NR==1 { print; next } $1==chr { print }' OFS='\t' "${targetDir}"/snplevel.sumstats4clumping.txt > "${targetDir}"/snplevel.clumping/sumstats4clumping_"${i}".txt
	plink --bfile $(eval echo "${chrFilehandler}") \
	--keep "${targetDir}/snplevel.subs.txt" \
	--clump "${targetDir}"/snplevel.clumping/sumstats4clumping_"${i}".txt \
	--clump-snp-field "ID" \
	--clump-field "P" \
	--clump-p1 1.000 \
	--clump-p2 1.000 \
	--clump-r2 0.1 \
	--clump-kb 10000 \
	--clump-verbose \
	--out "${targetDir}"/snplevel.clumping/chr"${i}"
	rm -f "${targetDir}"/snplevel.clumping/sumstats4clumping_chr"${i}".txt
	) &
done
wait)

# Merge clumped files
cat "${targetDir}"/snplevel.clumping/*.clumped > "${targetDir}"/snplevel.clumping/clumped.tmp
awk '$11 ~ /^[0-9]+$/ { indexsnp = $3; print $3, $3, 0, 1, "A1/A2"; next } $1=="" || $1=="KB" || $1=="CHR" || $1=="(INDEX)" || $2 == "not" || $1 == "RANGE:" || $1 == "SPAN:" || $1 ~ /-/ { next } {print indexsnp, $1, $2, $3, $4}' OFS='\t' "${targetDir}"/snplevel.clumping/clumped.tmp > "${targetDir}"/snplevel.clumping/clumped 
rm -f "${targetDir}"/snplevel.clumping/clumped.tmp

# add lead association (across phenotypes) to conditional analysis results
header=$(echo "$(head -1 "${targetDir}/snplevel.conditional.txt" | sed 's/LEAD_SNP_/COJO_/g')"$'\t'LEAD_SNP$'\t'LEAD_SNP_KB$'\t'LEAD_SNP_RSQ$'\t'LEAD_SNP_ALLELES)
awk -F'\t' -v header="${header}" 'BEGIN { print header } NR==FNR { snp[$2]=$1"\t"$3"\t"$4"\t"$5; next} FNR==1 { next } $4 in snp { print $0, snp[$4] }' OFS="\t" "${targetDir}"/snplevel.clumping/clumped "${targetDir}"/snplevel.conditional.txt > "${targetDir}"/snplevel.txt

# replace "A1/A2" (if SNP is Lead SNP) by actual alleles 
awk -F'\t' '$NF == "A1/A2" { output=$1; for(i=2;i<NF;++i) { output=output"\t"$i }; output=output"\t"$5"/"$6; print output; next} { print }' OFS="\t" "${targetDir}"/snplevel.txt > "${targetDir}"/snplevel.txt.tmp
\mv "${targetDir}"/snplevel.txt.tmp "${targetDir}"/snplevel.txt

# add CHR, BP, and GWS status of lead association
header=$(echo "$(head -1 "${targetDir}"/snplevel.txt)"$'\t'LEAD_SNP_CHR$'\t'LEAD_SNP_BP$'\t'LEAD_SNP_GWS)  
awk -F'\t' -v header="${header}" -v cols="CHR,BP,ID,LEAD_SNP,P,COJO_pJ" '
	BEGIN { print header; ncols=split(cols,colnames,",") }
	NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } }; next } 
	NR==FNR && $colidx[3]==$colidx[4] && (!($colidx[4] in pval) || $colidx[6] < pval[$colidx[4]]) { pval[$colidx[4]]=$colidx[6]; ($colidx[6] < 5E-8 && $colidx[5]<5E-8) ? gws=1 : gws=2; chrbp[$colidx[4]]=$colidx[1]"\t"$colidx[2]"\t"gws; next }
	NR==FNR || FNR==1 { next }
	{ print $0, chrbp[$colidx[4]] }' OFS="\t" "${targetDir}"/snplevel.txt "${targetDir}"/snplevel.txt > "${targetDir}"/snplevel.txt.tmp
\mv "${targetDir}"/snplevel.txt.tmp "${targetDir}"/snplevel.txt

# sort by LEAD_SNP_GWS, LEAD_SNP_CHR, LEAD_SNP_BP, and P
sortby=($(awk -F'\t' -v cols="LEAD_SNP_GWS,LEAD_SNP_CHR,LEAD_SNP_BP,P" 'BEGIN { ncols=split(cols,colnames,",") } NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } }; print colidx[1],colidx[2],colidx[3],colidx[4] }' "${targetDir}/snplevel.txt"))
awk -F'\t' -v chr="${sortby[1]}" '
	$chr=="X" { $chr=23 } $chr=="Y" { $chr=24 } $chr=="XY" { $chr=25 } $chr=="MT" { $chr=26 } { print }' OFS="\t" "${targetDir}/snplevel.txt" | 
	sort -t$'\t' -k${sortby[0]},${sortby[0]}g -k${sortby[1]},${sortby[1]}g -k${sortby[2]},${sortby[2]}g -k${sortby[3]},${sortby[3]}g > "${targetDir}/snplevel.txt.tmp"
\mv "${targetDir}"/snplevel.txt.tmp "${targetDir}"/snplevel.txt

# add LOCUS_COUNT and DISCOVERY_COUNT
header=$(echo LOCUS_COUNT$'\t'DISCOV_COUNT$'\t'"$(head -1 "${targetDir}/snplevel.txt")") 
awk -F'\t' -v header="$header" -v cols="LEAD_SNP" '
	BEGIN { print header; leadsnp="null"; locus_count = 0; snp_count = 1; ncols=split(cols,colnames,",") }
	NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } }; next } 
	leadsnp != $colidx[1] { leadsnp = $colidx[1]; locus_count++; snp_count = 1; print locus_count, snp_count, $0; next }
	{ snp_count++; print locus_count, snp_count, $0; next }' OFS='\t' "${targetDir}"/snplevel.txt > "${targetDir}"/snplevel.txt.tmp
\mv "${targetDir}"/snplevel.txt.tmp "${targetDir}"/snplevel.txt

# add sbayes fine-mapping results
if [[ ! -z "${sbayesFiles}" ]]; then
	sbayesFiles=($sbayesFiles)
	header=$(echo "$(head -1 "${targetDir}"/snplevel.txt)"$'\t'sbayes_numCausal$'\t'sbayes_pep$'\t'sbayes_cssizes$'\t'sbayes_genes) 
	rm -f "${targetDir}"/snplevel.txt.tmp
	for (( i=0; i<${#sbayesFiles[@]}; i++ )); do
		awk -F'\t' -v header="$header" -v trait="${traitlist[$i]}" ' BEGIN { print header } NR == FNR { sbayes[$1]=$9"\t"$7"\t"$3"\t"$10; next } FNR==1 { next } 
		$6 in sbayes && $3 == trait { print $0, sbayes[$6]; next } { print $0 }
		' OFS='\t' "${sbayesFiles[$i]}" "${targetDir}"/snplevel.txt > "${targetDir}"/snplevel.txt.tmp
		\mv "${targetDir}"/snplevel.txt.tmp "${targetDir}"/snplevel.txt
	done
	awk -F'\t' 'NR==1 { fieldnum=NF; print; next } NF < fieldnum { print $0, "NA", "NA", "NA", "NA"; next} { print }' OFS='\t' "${targetDir}"/snplevel.txt > "${targetDir}"/snplevel.txt.tmp
	\mv "${targetDir}"/snplevel.txt.tmp "${targetDir}"/snplevel.txt
fi

# add SusieR fine-mapping results
if [[ ! -z "${susieRFiles}" ]]; then
	susieRFiles=($susieRFiles)
	header=$(echo "$(head -1 "${targetDir}"/snplevel.txt)"$'\t'susieR_snpCausal$'\t'susieR_minPurity$'\t'susieR_cssizes$'\t'susieR_genes) 
	rm -f "${targetDir}"/snplevel.txt.tmp
	for (( i=0; i<${#susieRFiles[@]}; i++ )); do
		awk -F'\t' -v header="$header" -v trait="${traitlist[$i]}" ' BEGIN { print header } NR == FNR { susieR[$1]=$7"\t"$8"\t"$9"\t"$11; next } FNR==1 { next } 
		$6 in susieR && $3 == trait { print $0, susieR[$6]; next } { print $0 }
		' OFS='\t' "${susieRFiles[$i]}" "${targetDir}"/snplevel.txt > "${targetDir}"/snplevel.txt.tmp
		\mv "${targetDir}"/snplevel.txt.tmp "${targetDir}"/snplevel.txt
	done
	awk -F'\t' 'NR==1 { fieldnum=NF; print; next } NF < fieldnum { print $0, "NA", "NA", "NA", "NA"; next} { print }' OFS='\t' "${targetDir}"/snplevel.txt > "${targetDir}"/snplevel.txt.tmp
	\mv "${targetDir}"/snplevel.txt.tmp "${targetDir}"/snplevel.txt
fi

# add FINEMAP results
if [[ ! -z "${finemapFiles}" ]]; then
	susieRFiles=($finemapFiles)
	header=$(echo "$(head -1 "${targetDir}"/snplevel.txt)"$'\t'finemap_snpCausal$'\t'finemap_bestK$'\t'finemap_cssize$'\t'finemap_genes) 
	rm -f "${targetDir}"/snplevel.txt.tmp
	for (( i=0; i<${#susieRFiles[@]}; i++ )); do
		awk -F'\t' -v header="$header" -v trait="${traitlist[$i]}" ' BEGIN { print header } NR == FNR { susieR[$1]=$9"\t"$10"\t"$11"\t"$13; next } FNR==1 { next } 
		$6 in susieR && $3 == trait { print $0, susieR[$6]; next } { print $0 }
		' OFS='\t' "${susieRFiles[$i]}" "${targetDir}"/snplevel.txt > "${targetDir}"/snplevel.txt.tmp
		\mv "${targetDir}"/snplevel.txt.tmp "${targetDir}"/snplevel.txt
	done
	awk -F'\t' 'NR==1 { fieldnum=NF; print; next } NF < fieldnum { print $0, "NA", "NA", "NA", "NA"; next} { print }' OFS='\t' "${targetDir}"/snplevel.txt > "${targetDir}"/snplevel.txt.tmp
	\mv "${targetDir}"/snplevel.txt.tmp "${targetDir}"/snplevel.txt
fi

# add nonsynonymous exonic variations
if [[ ! -z ${nonsynFiles} ]]; then
	nonsynFiles=($nonsynFiles)
	header=$(echo "$(head -1 "${targetDir}"/snplevel.txt)"$'\t'nonsyn_customPthresh) 
	rm -f "${targetDir}"/snplevel.txt.tmp
	for (( i=0; i<${#nonsynFiles[@]}; i++ )); do
		awk -F'\t' -v header="$header" -v trait="${traitlist[$i]}" ' BEGIN { print header } NR == FNR { nonsyn[$1]=$2; next } FNR==1 { next } 
		$6 in nonsyn && $3 == trait { print $0, nonsyn[$6]; next } { print $0 }
		' OFS='\t' "${nonsynFiles[$i]}" "${targetDir}"/snplevel.txt > "${targetDir}"/snplevel.txt.tmp
		\mv "${targetDir}"/snplevel.txt.tmp "${targetDir}"/snplevel.txt
	done
	awk -F'\t' 'NR==1 { fieldnum=NF; print; next } NF < fieldnum { print $0, "NA"; next} { print }' OFS='\t' "${targetDir}"/snplevel.txt > "${targetDir}"/snplevel.txt.tmp
	\mv "${targetDir}"/snplevel.txt.tmp "${targetDir}"/snplevel.txt
fi
		
# add SMR results
if [[ ! -z ${smrFiles} ]]; then
	smrFiles=($smrFiles)
	header=$(echo "$(head -1 "${targetDir}"/snplevel.txt)"$'\t'SMR_eqtl$'\t'SMR_sqtl) 
	rm -f "${targetDir}"/snplevel.txt.tmp
	for (( i=0; i<${#smrFiles[@]}; i++ )); do
		awk -F'\t' -v header="$header" -v trait="${traitlist[$i]}" ' BEGIN { print header } NR == FNR { smr[$2]=$3"\t"$4; next } FNR==1 { next } 
		$6 in smr && $3 == trait { print $0, smr[$6]; next } { print $0 }
		' OFS='\t' "${smrFiles[$i]}" "${targetDir}"/snplevel.txt > "${targetDir}"/snplevel.txt.tmp
		\mv "${targetDir}"/snplevel.txt.tmp "${targetDir}"/snplevel.txt
	done
	awk -F'\t' 'NR==1 { fieldnum=NF; print; next } NF < fieldnum { print $0, "NA", "NA"; next} { print }' OFS='\t' "${targetDir}"/snplevel.txt > "${targetDir}"/snplevel.txt.tmp
	\mv "${targetDir}"/snplevel.txt.tmp "${targetDir}"/snplevel.txt
fi

# add GTEx eQTL matches
if [[ ! -z ${eqtlFiles} ]]; then
	eqtlFiles=($eqtlFiles)
	header=$(echo "$(head -1 "${targetDir}"/snplevel.txt)"$'\t'GTEx_singleTissue$'\t'GTEx_multiTissue) 
	rm -f "${targetDir}"/snplevel.txt.tmp
	for (( i=0; i<${#eqtlFiles[@]}; i++ )); do
		awk -F'\t' -v header="$header" -v trait="${traitlist[$i]}" ' BEGIN { print header } NR == FNR { eqtl[$2]=$3"\t"$4; next } FNR==1 { next } 
		$6 in eqtl && $3 == trait { print $0, eqtl[$6]; next } { print $0 }
		' OFS='\t' "${eqtlFiles[$i]}" "${targetDir}"/snplevel.txt > "${targetDir}"/snplevel.txt.tmp
		\mv "${targetDir}"/snplevel.txt.tmp "${targetDir}"/snplevel.txt
	done
	awk -F'\t' 'NR==1 { fieldnum=NF; print; next } NF < fieldnum { print $0, "NA", "NA"; next} { print }' OFS='\t' "${targetDir}"/snplevel.txt > "${targetDir}"/snplevel.txt.tmp
	\mv "${targetDir}"/snplevel.txt.tmp "${targetDir}"/snplevel.txt
fi

# add polygenic priority scores (PoPS)
if [[ ! -z ${popsFiles} ]]; then
	popsFiles=($popsFiles)
	header=$(echo "$(head -1 "${targetDir}"/snplevel.txt)"$'\t'PoPS) 
	rm -f "${targetDir}"/snplevel.txt.tmp
	for (( i=0; i<${#popsFiles[@]}; i++ )); do
		awk -F'\t' -v header="$header" -v trait="${traitlist[$i]}" ' BEGIN { print header } NR == FNR { pops[$1]=$4; next } FNR==1 { next } 
		$6 in pops && $3 == trait { print $0, pops[$6]; next } { print $0 }
		' OFS='\t' "${popsFiles[$i]}" "${targetDir}"/snplevel.txt > "${targetDir}"/snplevel.txt.tmp
		\mv "${targetDir}"/snplevel.txt.tmp "${targetDir}"/snplevel.txt
	done
	awk -F'\t' 'NR==1 { fieldnum=NF; print; next } NF < fieldnum { print $0, "NA"; next} { print }' OFS='\t' "${targetDir}"/snplevel.txt > "${targetDir}"/snplevel.txt.tmp
	\mv "${targetDir}"/snplevel.txt.tmp "${targetDir}"/snplevel.txt
fi

# add GWAS catalog matches
if [[ ! -z ${catalogFiles} ]]; then
	catalogFiles=($catalogFiles)
	header=$(echo "$(head -1 "${targetDir}/snplevel.txt")"$'\t'GWAS_CATALOG) 
	rm -f "${targetDir}/snplevel.txt.tmp"
	for (( i=0; i<${#catalogFiles[@]}; i++ )); do
		awk -F'\t' -v header="$header" -v trait="${traitlist[$i]}" ' BEGIN { print header } NR == FNR { catalog[$2]=$3; next } FNR==1 { next } 
		$6 in catalog && $3 == trait { print $0, catalog[$6]; next } { print $0 }
		' OFS='\t' "${catalogFiles[$i]}" "${targetDir}"/snplevel.txt > "${targetDir}"/snplevel.txt.tmp
		\mv "${targetDir}"/snplevel.txt.tmp "${targetDir}"/snplevel.txt
	done
	awk -F'\t' 'NR==1 { fieldnum=NF; print; next } NF < fieldnum { print $0, "NA"; next} { print }' OFS='\t' "${targetDir}"/snplevel.txt > "${targetDir}"/snplevel.txt.tmp
	\mv "${targetDir}"/snplevel.txt.tmp "${targetDir}"/snplevel.txt
fi

# create seperate file with gwas hits only
awk -F'\t' 'NR == 1 { print; for(i=1; i<=NF; i++) { if($i=="COJO_pJ") pJ=i; if($i=="P") P=i } }
	$P < 5E-8 && $pJ < 5E-8 { print }' OFS='\t' "${targetDir}"/snplevel.txt > "${targetDir}"/snplevel.gws.txt

# clean up
echo "Cleaning up."
rm -f "${targetDir}"/snplevel.conditional.txt
rm -f "${targetDir}"/snplevel.sumstats.txt
rm -f "${targetDir}"/snplevel.subs.txt
rm -f "${targetDir}"/snplevel.sumstats4clumping.txt
rm -f "${targetDir}"/snplevel.clumping/sumstats4clumping_*
tar -cvzf "${targetDir}"/snplevel.clumping.tar.gz --directory="${targetDir}" snplevel.clumping --remove-files
chmod -R 770 "${targetDir}"
echo "--- Combining snp-level results completed. --- "


