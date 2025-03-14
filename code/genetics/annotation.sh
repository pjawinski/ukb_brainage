#!/bin/bash

# ==========================
# === variant annotation ===
# ==========================

# get arguments
inputFile="${1}" # inputFile="results/gap_gm/credibleSet/credibleSet.df.txt"
outputFile="${2}" # outputFile="results/gap_gm/credibleSet/credibleSet.df.annot.txt"
chr_bp_ref_alt="${3}" # chr_bp_ref_alt="chromosome,bp,allele1,allele2"
humandb="${4}" # humandb="data/annovar/humandb"
refseq="${5}" # refseq="data/refseq/GRCh37_latest_genomic.edit.gff.gz"

# echo settings
echo $'\n'"--- variant annotation | settings ---"
echo "inputFile: ${inputFile}"
echo "outputFile: ${outputFile}"
echo "chr_bp_ref_alt: ${chr_bp_ref_alt}"
echo "humandb: ${humandb}"
echo "refseq: ${refseq}"$'\n'

# set targetdir and make folder
outputPrefix=$(echo "${outputFile}" | sed 's/.txt//')
targetDir=$(echo "${outputFile}" | sed 's/\/[^/]*$/\//')
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# transform file to vcf format CHR START END REF ALT
echo "Creating input file."
awk -F'\t' -v cols="${chr_bp_ref_alt}" '
	BEGIN { ncols=split(cols,colnames,",") } # get input column names
    NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } } } # get input column indices
    NR==1 { print "CHR", "START", "END", "REF", "ALT", $0; next }
    $colidx[1]=="XY" { $colidx[1]="X" }
    { print $colidx[1], $colidx[2]-(length($colidx[3])-1), $colidx[2], $colidx[3], $colidx[4], $0 }' OFS='\t' "${inputFile}" > "${targetDir}"/input4annovar.txt

# make annotation
echo "Running annotation."
mkdir -p "${targetDir}"/annovar/refGene
annotate_variation.pl -geneanno -dbtype refGene -buildver hg19 -out "${targetDir}"/annovar/refGene/output.annovar "${targetDir}"/input4annovar.txt "${humandb}"
mkdir -p "${targetDir}"/annovar/cytoBand
annotate_variation.pl -regionanno -dbtype cytoBand -buildver hg19 -out "${targetDir}"/annovar/cytoBand/output.annovar "${targetDir}"/input4annovar.txt "${humandb}"
mkdir -p "${targetDir}"/annovar/multianno/
table_annovar.pl "${targetDir}"/input4annovar.txt "${humandb}" -buildver hg19 -out "${targetDir}"/annovar/multianno/output.annovar -remove -protocol refGene,cytoBand,dbnsfp35a -operation gx,r,f -nastring NA -polish

# identify nearest gene and create clumped_sparse_annovar.txt file
header=$(echo "$(head -1 "${inputFile}")"$'\t'"REGION"$'\t'"NEAREST_GENE"$'\t'"DISTANCE"$'\t'"GENES")
awk -F'\t' -v header="$header" '
	BEGIN { print header } {
		gene_compl=$2; gene1=$2; gene2=$2; dist1=$2; dist2=$2;
		gsub(/[(].*/,"",gene1); gsub(/[,].*/,"",gene1); gsub(/.*[,]/,"",gene2); gsub(/[(].*/,"",gene2);
		gsub(/[)].*/,"",dist1); gsub(/.*[=]/,"",dist1); gsub(/.*=/,"",dist2); gsub(/[)].*/,"",dist2) }
		dist1 - dist2 < 0 { nearest_gene = gene1; distance = dist1/1 } dist1 - dist2 > 0 { nearest_gene = gene2; distance = dist2/1  } dist1 - dist2 == 0 { nearest_gene = gene1; distance = 0 }
		nearest_gene == "NONE" && dist1 == "NONE" { nearest_gene = gene2; distance = dist2/1 }
		nearest_gene == "NONE" && dist2 == "NONE" { nearest_gene = gene1; distance = dist1/1 }
		nearest_gene == "NONE" { nearest_gene = "Chr"$3":"$5; distance = 0 }
	{ for(i=8;i<=NF;++i) { length(output)==0 ? output = $i : output = output"\t"$i }; output = output"\t"$1"\t"nearest_gene"\t"distance"\t"$2; print output; output=""; next }
	' OFS='\t' "${targetDir}"/annovar/refGene/output.annovar.variant_function > "${outputFile}"

# add cytoband
echo "Adding cytoband."
header=$(echo "$(head -1 "${outputFile}")"$'\t'"CYTOBAND")
awk -F'\t' -v cols="${chr_bp_ref_alt}" -v header="$header" '
	BEGIN { print header; ncols=split(cols,colnames,",") }  # get input column names
	NR==FNR { id[$3":"$5":"$6":"$7]=$2; next }
	FNR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } }; next} # get input column indices
	{ print $0, id[$colidx[1]":"$colidx[2]":"$colidx[3]":"$colidx[4]] }' OFS='\t' "${targetDir}"/annovar/cytoBand/output.annovar.hg19_cytoBand "${outputFile}" > "${outputFile}".tmp
	\mv "${outputFile}.tmp" "${outputFile}"

# add biotype and description
echo "Adding nearest gene biotype and description."
header=$(echo "$(head -1 "${outputFile}")"$'\t'"NEAREST_GENE_BIOTYPE"$'\t'"NEAREST_GENE_DESCRIPTION")
awk -F'\t' -v header="${header}" '
	BEGIN { print header }
	NR==FNR { gene[$10]=$11"\t"$12; next }
	FNR==1 { next }
	{ current_gene = $(NF-3); gsub(/[,].*/,"",current_gene) }
	current_gene in gene { print $0, gene[current_gene]; next }
	{ print $0, "NA", "NA"; next }' OFS='\t' <(gzip -dc "${refseq}") "${outputFile}" > "${outputFile}".tmp
	\mv "${outputFile}".tmp "${outputFile}"

	# no match? get refseq rows that contain id as synonym
	#echo "(2/4) Finding matches in RefSeq synonyms..."
	missings=$(awk -F'\t' '$NF == "NA" { gene = $(NF-5); gsub(/[,].*/,"",gene); print gene}' "${outputFile}")

	if [ ! -z "$missings" ]; then

		# get rownum of identified synonyms
		awk -F'\t|,' 'NR==FNR { missings[$1]; next}
			{ for(k=1; k<=NF; k++) {if($k !~ /^[0-9]+$/ && $k in missings) { print FNR, $k; delete missings[$k] } } }
			' OFS='\t' <(echo "$missings") <(gzip -dc "${refseq}") > "${targetDir}/annovar/synonyms.rownum.txt"

		# for identified synonyms get refseq id, biotype, and description
		nrows=$(awk 'END { print NF } ' "${targetDir}/annovar/synonyms.rownum.txt")
 		if [[ $nrows != 0 ]]; then
			awk -F'\t' 'NR==FNR { id[$1]=$2; next } FNR in id { print id[FNR], $10, $11, $12, $13}' OFS='\t' "${targetDir}/annovar/synonyms.rownum.txt" <(gzip -dc "${refseq}") > "${targetDir}"/annovar/synonyms.annotated.txt

		# merge with annotation list
		header=$(echo "$(head -1 "${outputFile}")")
		awk -F'\t' -v header="$header" 'BEGIN { print header }
			NR==FNR { id[$1]=$3"\t"$4; next }
			FNR==1 { next } 
			$(NF-5) in id { for(i=1;i<=NF-2;++i) { length(output)==0 ? output = $i : output = output"\t"$i }; output = output"\t"id[$(NF-5)]; print output; output=""; next }
			{ print $0}' OFS='\t' "${targetDir}/annovar/synonyms.annotated.txt" "${outputFile}" > "${outputFile}".tmp
		\mv "${outputFile}".tmp "${outputFile}"
		fi
	fi

# create results for nonsynonymous exonic variants
header=$(echo "$(head -1 "${outputFile}")"$'\t'EXONIC_FUNCTION$'\t'TRANSCRIPT_CONSEQUENCE$'\t'CADD_PHRED$'\t'CADD_RAW$'\t'CADD_RANK$'\t'DANN$'\t'DANN_RANK$'\t'REVEL$'\t'REVEL_RANK)
awk -F'\t' -v cols="${chr_bp_ref_alt}" -v header="${header}" '
	BEGIN { print header; ncols=split(cols,colnames,",") }  # get input column names 
	NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } }; next} # get input column indices
	NR==FNR { id[$colidx[1]":"$colidx[2]":"$colidx[3]":"$colidx[4]]=$0; next}
	$1":"$3":"$4":"$5 in id {print id[$1":"$3":"$4":"$5],$9,$10,$53,$51,$52,$54,$55,$47,$48; next}' \
	OFS='\t' "${outputFile}" <(grep 'nonsynonymous\|frameshift\|stopgain\|stoploss' "${targetDir}"/annovar/multianno/output.annovar.hg19_multianno.txt) \
	> "${outputPrefix}".nonsynonymous.txt

# clean up
echo "Cleaning up."
rm -f "${targetDir}"/input4annovar.txt
tar -cvzf "${outputFile}".tar.gz --directory="${targetDir}" annovar --remove-files
chmod -R 770 "${outputPrefix}"*
echo "--- Completed: variant annotation. --- "

