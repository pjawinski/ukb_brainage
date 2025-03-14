#!/bin/bash

# ============================
# === Conditional analysis ===
# ============================

# get arguments
subsFile="${1}" # subsFile="data/genetics/chr1/imp_mri_qc_EUR/chr1_mri_qc.psam"
targetDir="${2}" # targetDir="results/${trait}/replicate/conditional"
chrFileHandle="${3}" # chrFileHandle="data/genetics/chr\$i/imp_mri_qc_EUR/bed/chr\${i}_mri_qc"
sumstats="${4}" # sumstats="results/${trait}/replicate/metal.ivweight.qc.gz"
cols="${5}" # set input column names that correspond to header column names (see below) | cols="ID,A1,A2,A1_FREQ,BETA,SE,P,N,CHR,BP,A1_FREQ_SE,A1_FREQ_MIN,A1_FREQ_MAX,DIRECTION,HET_ISQ,HetChiSq,HET_DF,HET_P"
sumstatsPLINK="${6}" # sumstatsPLINK="results/${trait}/replicate/EUR/sumstats_plink.txt.gz"
pthresh=${7} # pthresh=1E-6
threads=${8} # threads=100
humandb="${9}" # humandb="data/annovar/humandb"
refseq="${10}" # refseq="data/refseq/GRCh37_latest_genomic.edit.gff.gz"

# echo settings
echo $'\n'"--- Conditional Analysis ---"
echo "subsFile: ${subsFile}"
echo "targetDir: ${targetDir}"
echo "chrFileHandle: ${chrFileHandle}"
echo "sumstats: ${sumstats}"
echo "cols: ${cols}"
echo "sumstatsPLINK: ${sumstatsPLINK}"
echo "pthresh: ${pthresh}"
echo "threads: ${threads}"
echo "humandb: ${humandb}"
echo "refseq: ${refseq}"$'\n'

# set targetdir and make folder
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# create input file from sumstats
echo "Creating input file."
awk -v cols="${cols}" '
	BEGIN { ncols=split(cols,colnames,",") } # get input column names
    NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } } } # get input column indices
    { for(i=1;i<=ncols;++i) { length(output)==0 ? output = $colidx[i] : output = output"\t"$colidx[i] }; print output; output=""; next }' OFS='\t' <(gzip -dc "${sumstats}") > "${targetDir}"/input.txt

# get list of chr with significant hits
echo "Getting list of chromosomes with significant hits."
chr=$(awk -v pthresh="${pthresh}" '$7 < pthresh { print $9 }' "${targetDir}"/input.txt | sort -u)
chr=($chr)

# For gcta software: use alternative .bim file with chromsomes converted to numeric values
for i in X Y XY MT; do
	eval \\cp "${chrFileHandle}"_numeric.bim "${chrFileHandle}".bim
done

# run conditional analysis
echo "Running conditional analysis."
mkdir -p "${targetDir}"/output
chrCount=$(echo "${#chr[@]}")
threadsPerAnalysis=$(expr "${threads}" / "${chrCount}")

(
for i in ${chr[@]}; do (
gcta64 --bfile $(eval echo "${chrFileHandle}") \
--cojo-file "${targetDir}"/input.txt \
--cojo-slct \
--cojo-p "${pthresh}" \
--cojo-wind 10000 \
--cojo-collinear 0.9 \
--out "${targetDir}"/output/cond_chr"${i}" \
--thread-num "${threadsPerAnalysis}"
) &
done
wait
)

# replace numeric .bim files with originals
for i in X Y XY MT; do
	eval \\cp "${chrFileHandle}"_original.bim "${chrFileHandle}".bim
done

# get summary of results
awk 'NR == 1 { print; next} FNR > 1 {sub(23,"X",$1); sub(24,"Y",$1); sub(25,"XY",$1); sub(26,"MT",$1); print $0}' OFS="\t" "${targetDir}"/output/cond_chr*.jma.cojo > "${targetDir}"/conditional.txt

# clean list (variants are not independent if p value increases by two orders of magnitude AND if SNP is within a 10 MB window of a more significantly associated SNP)
sort -k1,1n -k8,8g "${targetDir}"/conditional.txt | awk -F'\t' 'BEGIN {chr = ""} NR == 1 { print; next}
$1 != chr { print; chr = $1; delete bp; k = 1; bp[k] = $3; next }
$8*100<$13 { k = k + 1; bp[k] = $3; min = 10000001; for (i = 1; i < k; i++) { if (sqrt(($3 - bp[i])^2) < min) min = sqrt(($3 - bp[i])^2) }; { if (min <= 10000000) next } }
{ k = k + 1; bp[k] = $3; print }' > "${targetDir}"/conditional.cleaned.txt

sort -k1,1n -k8,8g "${targetDir}"/conditional.txt | awk -F'\t' 'BEGIN {chr = ""; k = 1 } NR == 1 { print; next}
$1 != chr { chr = $1; delete bp; k = 1; bp[k] = $3; next }
$8*100<$13 { k = k + 1; bp[k] = $3; min = 10000001; for (i = 1; i < k; i++) { if (sqrt(($3 - bp[i])^2) < min) min = sqrt(($3 - bp[i])^2) }; { if (min <= 10000000) print } }
{ k = k + 1; bp[k] = $3 }' > "${targetDir}"/conditional.dropped.txt

# ==========================================================================================================
# === get 'tophits', i.e., proxy variants in strong ld with conditional discoveries and with p < pthresh ===
# ==========================================================================================================

# get ID and P values of all variations from summary statistics
awk '{ print $1, $7}' OFS='\t' "${targetDir}"/input.txt > "${targetDir}"/sumstats4clumping.txt

# get participant IDs (they serve as reference sample for clumping)
awk '{ print $1, $2 }' "${subsFile}" > "${targetDir}"/subs.txt

# get conditional discoveries
chr=$(awk 'NR==1 { next } { print $1 }' "${targetDir}"/conditional.cleaned.txt) 
chr=($chr)
snp=$(awk 'NR==1 { next } { print $2 }' "${targetDir}"/conditional.cleaned.txt) 
snp=($snp)

# do clumping
echo "Getting tophits through clumping."
mkdir -p "${targetDir}"/clumping

(for j in $(seq 0 $((${#snp[@]}-1))); do
	awk -v snp="${snp[j]}" 'NR==1 { print; next } $1==snp { print $1, 1E-307; next } { print }' OFS='\t' "${targetDir}"/sumstats4clumping.txt > "${targetDir}"/clumping/sumstats4clumping_"${snp[j]}".txt
	i=${chr[j]}
	LD_src=$(eval echo "${chrFileHandle}")
	plink --bfile "$LD_src" \
		--keep "${targetDir}"/subs.txt \
		--clump "${targetDir}"/clumping/sumstats4clumping_"${snp[j]}".txt \
		--clump-snp-field "ID" \
		--clump-field "P" \
		--clump-p1 1E-300 \
		--clump-p2 "${pthresh}" \
		--clump-r2 0.8 \
		--clump-kb 3000 \
		--clump-verbose  \
		--out "${targetDir}"/clumping/"${snp[j]}"
	rm -rf "${targetDir}"/sumstats4clumping_"${snp[j]}".txt
	done
	wait
)

# Merge clumped files
cat "${targetDir}"/clumping/*.clumped > "${targetDir}"/clumping/clumped.tmp
awk '$11 ~ /^[0-9]+$/ { indexsnp = $3; print $3, $3, 0, 1, "A1/A2"; next } $1=="" || $1=="KB" || $1=="CHR" || $1=="(INDEX)" || $2 == "not" || $1 == "RANGE:" || $1 == "SPAN:" || $1 ~ /-/ { next } {print indexsnp, $1, $2, $3, $4}' OFS='\t' "${targetDir}"/clumping/clumped.tmp > "${targetDir}"/clumping/clumped ; rm -f "${targetDir}"/clumping/clumped.tmp

# create conditional.cleaned.tophits file (sort results by LEAD_SNP_P, LEAD_SNP_RSQ, P)
awk 'NR==FNR { snp[$2]=$1"\t"$3"\t"$4"\t"$5; next } 
	 FNR==1 { ncols = NF }
	 $1 in snp { output=$9"\t"$10"\t"$1"\t"snp[$1]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;
	 			 for(i=11;i<=ncols;++i) { output = output"\t"$i }; print output; output=""; next }
	 ' OFS="\t" "${targetDir}"/clumping/clumped "${targetDir}"/input.txt > "${targetDir}"/clumping/tophits

awk 'NR==1 { ncols = NF }
	 $7 == "A1/A2" { output=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$8"/"$9;
					 for(i=8;i<=ncols;++i) { output = output"\t"$i }; print output; output=""; next }
	 { print }' OFS="\t" "${targetDir}"/clumping/tophits > "${targetDir}"/clumping/tophits.tmp
	 \mv "${targetDir}"/clumping/tophits.tmp "${targetDir}"/clumping/tophits

header="LOCUS_CNT"$'\t'"SNP_CNT"$'\t'"CHR"$'\t'"BP"$'\t'"ID"$'\t'"LEAD_SNP"$'\t'"LEAD_SNP_P"$'\t'"LEAD_SNP_KB"$'\t'"LEAD_SNP_RSQ"$'\t'"LEAD_SNP_ALLELES"$'\t'"A1"$'\t'"A2"$'\t'"A1_FREQ"$'\t'"BETA"$'\t'"SE"$'\t'"P"$'\t'"N"
header=$(echo "${header}"$'\t'"$(awk 'NR==1 { for(i=11;i<=NF;++i) { length(output)==0 ? output = $i : output = output"\t"$i }; print output; exit }' OFS='\t' "${targetDir}"/input.txt)")
awk -v pthresh="${pthresh}" '
	 NR==FNR && $7 <= pthresh { leadsnp_p[$1]=$7; next }
	 NR==FNR { next }
	 FNR==1 { ncols = NF }
	 $4 in leadsnp_p { output=$1"\t"$2"\t"$3"\t"$4"\t"leadsnp_p[$4];
	 	for(i=5;i<=ncols;++i) { output = output"\t"$i }; print output; output=""; next }' OFS="\t" "${targetDir}"/input.txt "${targetDir}"/clumping/tophits |
	 	sort -g -k5,5 -k7,7rn -k14,14 |
	 	awk -v header="${header}" '
	 		BEGIN { print header; leadsnp="null"; locus_cnt = 0; snp_cnt = 0 }
	 		leadsnp != $5 { leadsnp = $5; locus_cnt++; snp_cnt = 1; print locus_cnt, snp_cnt, $0; next }
	 		{ snp_cnt++; print locus_cnt, snp_cnt, $0 } ' OFS="\t" > "${targetDir}"/conditional.cleaned.tophits.txt

# ======================
# === run annotation ===
# ======================

# transform clumped sumstats to vcf format CHR START END REF A
mkdir -p "${targetDir}"/annovar
header=$(echo "$(head -1 "${targetDir}"/conditional.cleaned.tophits.txt)"$'\t'"REF"$'\t'"ALT")
awk -v header="$header" '
	BEGIN { print header } 
	FNR==1 { next }
	NR==FNR { variant[$5]=$0; next }
	$3 in variant { print variant[$3], $5, $4 }' OFS='\t' "${targetDir}"/conditional.cleaned.tophits.txt <(gzip -dc "${sumstatsPLINK}") > "${targetDir}"/annovar/conditional.cleaned.tophits4annovar.txt 

header=$(echo "CHR"$'\t'"START"$'\t'"END"$'\t'"REF"$'\t'"ALT"$'\t'"$(head -1 "${targetDir}"/conditional.cleaned.tophits.txt)")
awk -v header="$header" '
	BEGIN { print header }
	FNR==1 { next }
	{ CHR=$3 } CHR=="XY" { CHR="X" } {
		START = $4 - (length($(NF-1)) - 1);
		output=CHR"\t"START"\t"$4"\t"$(NF-1)"\t"$NF;
		for(i=1;i<=NF-2;++i) { output = output"\t"$i }; print output; output=""; next }' OFS='\t' "${targetDir}"/annovar/conditional.cleaned.tophits4annovar.txt | 
		sort -n -k6,6g -k7,7g > "${targetDir}"/annovar/conditional.cleaned.tophits4annovar.tmp.txt
		\mv "${targetDir}"/annovar/conditional.cleaned.tophits4annovar.tmp.txt "${targetDir}"/annovar/conditional.cleaned.tophits4annovar.txt

# make annotation
echo "Running annotation."
mkdir -p "${targetDir}"/annovar/refGene
annotate_variation.pl -geneanno -dbtype refGene -buildver hg19 -out "${targetDir}"/annovar/refGene/conditional.cleaned.tophits.annovar "${targetDir}"/annovar/conditional.cleaned.tophits4annovar.txt "${humandb}"

mkdir -p "${targetDir}"/annovar/cytoBand
annotate_variation.pl -regionanno -dbtype cytoBand -buildver hg19 -out "${targetDir}"/annovar/cytoBand/conditional.cleaned.tophits.annovar "${targetDir}"/annovar/conditional.cleaned.tophits4annovar.txt "${humandb}"

mkdir -p "${targetDir}"/annovar/multianno/
table_annovar.pl "${targetDir}"/annovar/conditional.cleaned.tophits4annovar.txt "${humandb}" -buildver hg19 -out "${targetDir}"/annovar/multianno/conditional.cleaned.tophits.annovar -remove -protocol refGene,cytoBand,dbnsfp35a -operation gx,r,f -nastring NA -polish

# identify nearest gene and create clumped_sparse_annovar.txt file
header=$(echo "$(head -1 "${targetDir}"/conditional.cleaned.tophits.txt)"$'\t'"REGION"$'\t'"NEAREST_GENE"$'\t'"DISTANCE"$'\t'"GENES")
awk -v header="$header" '
	BEGIN { print header } {
		gene_compl=$2; gene1=$2; gene2=$2; dist1=$2; dist2=$2;
		gsub(/[(].*/,"",gene1); gsub(/.*[,]/,"",gene2); gsub(/[(].*/,"",gene2);
		gsub(/[)].*/,"",dist1); gsub(/.*[=]/,"",dist1); gsub(/.*=/,"",dist2); gsub(/[)].*/,"",dist2) }
		dist1 - dist2 < 0 { nearest_gene = gene1; distance = dist1/1 } dist1 - dist2 > 0 { nearest_gene = gene2; distance = dist2/1  } dist1 - dist2 == 0 { nearest_gene = gene1; distance = 0 }
		nearest_gene == "NONE" && dist1 == "NONE" { nearest_gene = gene2; distance = dist2/1 }
		nearest_gene == "NONE" && dist2 == "NONE" { nearest_gene = gene1; distance = dist1/1 }
		nearest_gene == "NONE" { nearest_gene = "Chr"$3":"$5; distance = 0 }
	{ for(i=8;i<=NF;++i) { length(output)==0 ? output = $i : output = output"\t"$i }; output = output"\t"$1"\t"nearest_gene"\t"distance"\t"$2; print output; output=""; next }
	' OFS='\t' "${targetDir}"/annovar/refGene/conditional.cleaned.tophits.annovar.variant_function > "${targetDir}"/conditional.cleaned.tophits.annovar.txt

# add cytoband
echo "Adding cytoband."
header=$(echo "$(head -1 "${targetDir}"/conditional.cleaned.tophits.annovar.txt)"$'\t'"CYTOBAND")
awk -v header="$header" '
	BEGIN { print header } NR==FNR { snp[$8]=$2; next }
	FNR==1 { next }
	{ print $0, snp[$1] }' OFS='\t' "${targetDir}/"annovar/cytoBand/conditional.cleaned.tophits.annovar.hg19_cytoBand "${targetDir}"/conditional.cleaned.tophits.annovar.txt > "${targetDir}"/conditional.cleaned.tophits.annovar.tmp.txt
	\mv "${targetDir}"/conditional.cleaned.tophits.annovar.tmp.txt "${targetDir}"/conditional.cleaned.tophits.annovar.txt

# add biotype and description
echo "Adding nearest gene biotype and description."
header=$(echo "$(head -1 "${targetDir}"/conditional.cleaned.tophits.annovar.txt)"$'\t'"NEAREST_GENE_BIOTYPE"$'\t'"NEAREST_GENE_DESCRIPTION")
awk -F'\t' -v header="${header}" '
	BEGIN { print header }
	NR==FNR { gene[$10]=$11"\t"$12; next }
	FNR==1 { next }
	{ current_gene = $(NF-3); gsub(/[,].*/,"",current_gene) }
	current_gene in gene { print $0, gene[current_gene]; next }
	{ print $0, "NA", "NA"; next }' OFS='\t' <(gzip -dc "${refseq}") "${targetDir}"/conditional.cleaned.tophits.annovar.txt > "${targetDir}"/conditional.cleaned.tophits.annovar.tmp.txt
	\mv "${targetDir}"/conditional.cleaned.tophits.annovar.tmp.txt "${targetDir}"/conditional.cleaned.tophits.annovar.txt

	# no match? get refseq rows that contain id as synonym
	#echo "(2/4) Finding matches in RefSeq synonyms..."
	missings=$(awk -F'\t' '$NF == "NA" { gene = $(NF-5); gsub(/[,].*/,"",gene); print gene}' "${targetDir}"/conditional.cleaned.tophits.annovar.txt)

	if [ ! -z "$missings" ]; then

		# get rownum of identified synonyms
		awk -F'\t|,' 'NR==FNR { missings[$1]; next}
			{ for(k=1; k<=NF; k++) {if($k !~ /^[0-9]+$/ && $k in missings) { print FNR, $k; delete missings[$k] } } }
			' OFS='\t' <(echo "$missings") <(gzip -dc "${refseq}") > "${targetDir}"/annovar/synonyms.rownum.txt

		# for identified synonyms get refseq id, biotype, and description
		awk -F'\t' 'NR==FNR { id[$1]=$2; next } FNR in id { print id[FNR], $10, $11, $12, $13}' OFS='\t' "${targetDir}"/annovar/synonyms.rownum.txt <(gzip -dc "${refseq}") > "${targetDir}"/annovar/synonyms.annotated.txt

		# merge with annotation list
		header=$(echo "$(head -1 "${targetDir}"/conditional.cleaned.tophits.annovar.txt)")
		awk -F'\t' -v header="$header" 'BEGIN { print header }
			NR==FNR { id[$1]=$3"\t"$4; next }
			FNR==1 { next } 
			$(NF-5) in id { for(i=1;i<=(NF-2);++i) { length(output)==0 ? output = $i : output = output"\t"$i }; output = output"\t"id[$(NF-5)]; print output; output=""; next }
			{ print $0}' OFS='\t' "${targetDir}"/annovar/synonyms.annotated.txt "${targetDir}/"conditional.cleaned.tophits.annovar.txt > "${targetDir}"/conditional.cleaned.tophits.annovar.tmp.txt
			\mv "${targetDir}"/conditional.cleaned.tophits.annovar.tmp.txt "${targetDir}"/conditional.cleaned.tophits.annovar.txt
	fi

# add conditional analysis results (pvalue etc.)
header=$(echo "$(head -1 "${targetDir}"/conditional.cleaned.tophits.annovar.txt)"$'\t'"LEAD_SNP_bJ"$'\t'"LEAD_SNP_bJ_se"$'\t'"LEAD_SNP_pJ"$'\t'"LEAD_SNP_LD_r")
awk -F'\t' -v header="$header" '
	BEGIN { print header } NR==FNR { snp[$2]=$11"\t"$12"\t"$13"\t"$14; next }
	$6 in snp { print $0, snp[$6] }' OFS='\t' "${targetDir}"/conditional.cleaned.txt "${targetDir}"/conditional.cleaned.tophits.annovar.txt > "${targetDir}"/conditional.cleaned.tophits.annovar.tmp.txt
	\mv "${targetDir}"/conditional.cleaned.tophits.annovar.tmp.txt "${targetDir}"/conditional.cleaned.tophits.annovar.txt

# create results for nonsynonymous exonic variants
header=$(echo "$(head -1 "${targetDir}"/conditional.cleaned.tophits.annovar.txt)"$'\t'EXONIC_FUNCTION$'\t'TRANSCRIPT_CONSEQUENCE$'\t'CADD_PHRED$'\t'CADD_RAW$'\t'CADD_RANK$'\t'DANN$'\t'DANN_RANK$'\t'REVEL$'\t'REVEL_RANK)
awk -F'\t' -v header="${header}" 'BEGIN { print header } NR==FNR { id[$3":"$4":"$11":"$12]=$0; next}
	$1":"$2":"$4":"$5 in id {print id[$1":"$2":"$4":"$5],$9,$10,$53,$51,$52,$54,$55,$47,$48; next}
	$1":"$2":"$5":"$4 in id {print id[$1":"$2":"$5":"$4],$9,$10,$53,$51,$52,$54,$55,$47,$48; next}' \
	OFS='\t' "${targetDir}"/conditional.cleaned.tophits.annovar.txt <(grep 'nonsynonymous\|frameshift\|stopgain\|stoploss' "${targetDir}"/annovar/multianno/conditional.cleaned.tophits.annovar.hg19_multianno.txt) \
	> "${targetDir}"/conditional.cleaned.tophits.annovar.nonsynonymous.txt

# summarize results for nonsynonymous exonic variants
scriptDir=$(dirname "$0")
Rscript "${scriptDir}"/conditional.exonic.R "${targetDir}"/conditional.cleaned.tophits.annovar.nonsynonymous.txt "${pthresh}" "${targetDir}"/conditional.cleaned.tophits.annovar.nonsynonymous.summary.txt

# clean up
echo "Cleaning up."
mv "${targetDir}"/output/ "${targetDir}"/conditional.log/
rm -f "${targetDir}"/sumstats4clumping.txt
rm -f "${targetDir}"/input.txt
rm -f "${targetDir}"/subs.txt
rm -rf "${targetDir}"/conditional.log/*cma.cojo
rm -f "${targetDir}"/clumping/*clumped
rm -f "${targetDir}"/clumping/sumstats4clumping_*
rm -f "${targetDir}"/clumping/tophits
tar -cvzf "${targetDir}"/conditional.clumping.tar.gz --directory="${targetDir}" clumping --remove-files
tar -cvzf "${targetDir}"/conditional.annovar.tar.gz --directory="${targetDir}" annovar --remove-files
tar -cvzf "${targetDir}"/conditional.log.tar.gz --directory="${targetDir}" conditional.log --remove-files
chmod -R 770 "${targetDir}"/*
echo "--- Conditional analysis finished. --- "

