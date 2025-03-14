#!/bin/bash

# ========================
# === MAGMA annotation ===
# ========================

# get arguments
inputFile="${1}" # inputFile="results/${trait}/magma/magma.genes.out"
outputFile="${2}" # outputFile="results/${trait}/magma/magma.genes.out.annot"
annotationFile="${3}" # annotationFile="/fast/software/pops/data/utils/gene_annot_jun10.txt"
conditionalFile="${4}" # conditionalFile="results/${trait}/conditional/conditional.cleaned.tophits.annovar.txt"
pthreshMapping=${5} # pthreshMapping=5E-8
glist_hg19="${6}" # glist_hg19="data/glist_hg19/glist-hg19-refseq.txt"
threads=${7} # threads=100
window=${8} # window=0
humandb="${9}" # humandb="data/annovar/humandb"
refseq="${10}" # refseq="data/refseq/GRCh37_latest_genomic.edit.gff.gz"
clumpingWindow=${11} # clumpingWindow=3000

# echo settings
echo $'\n'"--- MAGMA gene-based annotation settings ---"
echo "inputFile: ${inputFile}"
echo "outputFile: ${outputFile}"
echo "annotationFile: ${annotationFile}"
echo "conditionalFile: ${conditionalFile}"
echo "pthreshMapping: ${pthreshMapping}"
echo "glist_hg19: ${glist_hg19}"
echo "threads: ${threads}"
echo "window: ${window}"
echo "humandb: ${humandb}"
echo "refseq: ${refseq}"
echo "clumpingWindow: ${clumpingWindow}"$'\n'

# sort results by p value
sort -g -k 9 "${inputFile}" | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' > "${outputFile}".tmp

# add gene symbol
header="SYMBOL"$'\t'"GENE"$'\t'"CHR"$'\t'"START"$'\t'"STOP"$'\t'"NSNPS"$'\t'"NPARAM"$'\t'"N"$'\t'"ZSTAT"$'\t'"P"
awk -v header="${header}" 'NR==1 { next } NR==FNR { symbol[$1]=$2; next }
     FNR==1 { print header; next }
     $1 in symbol { print symbol[$1], $0 }' OFS='\t' "${annotationFile}" "${outputFile}".tmp \
     > "${outputFile}"

# prepare for annovar
echo "Running annotation."
mkdir -p "${outputFile}".dir
header="CHR"$'\t'"START"$'\t'"END"$'\t'"REF"$'\t'"ALT"$'\t'"SYMBOL"$'\t'"GENE"$'\t'"CHR"$'\t'"START"$'\t'"STOP"$'\t'"NSNPS"$'\t'"NPARAM"$'\t'"N"$'\t'"ZSTAT"$'\t'"P"
awk -v header="$header" 'BEGIN { print header } FNR==1 { next } { CHR=$3 } CHR=="23" { CHR="X" } CHR=="24" { CHR="Y" } CHR=="25" { CHR="X" } CHR=="26" { CHR="MT" } { START=int($4+($5-$4)/2); print CHR, START, START, "A", "G", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10 }' OFS='\t' "${outputFile}" > "${outputFile}".dir/magma4annovar

# add cytoband
mkdir -p "${outputFile}".dir/cytoBand
annotate_variation.pl -regionanno -dbtype cytoBand -buildver hg19 -out "${outputFile}".dir/cytoBand/magma.genes.out.annovar "${outputFile}".dir/magma4annovar "${humandb}"
header=$(echo "$(head -1 "${outputFile}")"$'\t'"CYTOBAND")
awk -v header="$header" 'BEGIN { print header } NR==FNR { cyto[$8]=$2; next } FNR == 1 { next } !($1 in cyto) { print $0, "NA"; next} { print $0, cyto[$1] }' OFS='\t' "${outputFile}".dir/cytoBand/magma.genes.out.annovar.hg19_cytoBand "${outputFile}" > "${outputFile}".tmp
\mv "${outputFile}".tmp "${outputFile}"

# add biotype and description
header=$(echo "$(head -1 "${outputFile}")"$'\t'"GENE_DESCRIPTION"$'\t'"GENE_BIOTYPE")
awk -F'\t' -v header="$header" 'BEGIN { print header } NR==FNR { gene[$10]=$12"\t"$11; next } FNR==1 { next } $1 in gene { print $0, gene[$1]; next} { print $0, "NA", "NA"; next }' OFS='\t' <(gzip -dc "${refseq}") "${outputFile}" > "${outputFile}".tmp
\mv "${outputFile}".tmp "${outputFile}"

  # no match? get refseq rows that contain id as synonym
  #echo "(2/4) Finding matches in RefSeq synonyms..."
  missings=$(awk -F'\t' '$NF == "NA" { print $1}' "${outputFile}")

  if [ ! -z "$missings" ]; then

    awk -F'\t|,' 'NR==FNR { missings[$1]; next}
      { for(k=1; k<=NF; k++) {if($k !~ /^[0-9]+$/ && $k in missings) { print FNR, $k; delete missings[$k] } } }
      ' OFS='\t' <(echo "${missings}") <(gzip -dc "${refseq}") > "${outputFile}".dir/synonyms.rownum.txt

    # for identified synonyms get refseq id, biotype, and description
    awk -F'\t' 'NR==FNR { id[$1]=$2; next } FNR in id { print id[FNR], $10, $11, $12, $13}' OFS='\t' "${outputFile}".dir/synonyms.rownum.txt <(gzip -dc "${refseq}") > "${outputFile}".dir/synonyms.annotated.txt

    # merge with annotation list
    awk -F'\t' -v header="$header" 'BEGIN { print header }
      NR==FNR { id[$1]=$3"\t"$4; next }
      FNR==1 { next } 
      $1 in id { output=$1; { for(i=2; i<=NF-2; i++) { output=output"\t"$i } }; print output, id[$1]; next }
      { print $0}' OFS='\t' "${outputFile}".dir/synonyms.annotated.txt "${outputFile}" > "${outputFile}".tmp
    \mv "${outputFile}".tmp "${outputFile}"

  fi

# add ENTREZ
echo "Adding ENTREZ identifier."
header=$(echo "$(head -1 "${outputFile}")"$'\t'"ENTREZ_ID")
awk -F'\t|,|;' -v header="$header" 'BEGIN { print header }
  NR==FNR { gsub(/Dbxref=/,""); for(k=1; k<=NF; k++) { if($k ~ /^GeneID:/) { gene[$10]=$k } }; next }
  FNR==1 { next } $1 in gene { print $0, gene[$1]; next }
  { print $0, "NA"}' OFS='\t' <(gzip -dc "${refseq}") "${outputFile}" > "${outputFile}".tmp
\mv "${outputFile}".tmp "${outputFile}"

# add HGNC
echo "Adding HGNC identifier."
header=$(echo "$(head -1 "${outputFile}")"$'\t'"HGNC_ID")
awk -F'\t|,|;' -v header="$header" 'BEGIN { print header }
  NR==FNR { gsub(/Dbxref=/,""); gsub(/HGNC:HGNC:/,"HGNC:"); for(k=1; k<=NF; k++) { if($k ~ /^HGNC:/) { gene[$10]=$k } }; next }
  FNR==1 { next } $1 in gene { print $0, gene[$1]; next }
  { print $0, "NA"}' OFS='\t' <(gzip -dc "${refseq}") "${outputFile}" > "${outputFile}".tmp
\mv "${outputFile}".tmp "${outputFile}"

# add ".XY" to pseudoautosomal gene-based results
echo "Adding .XY to pseudoautosomal genes."
awk -F'\t' '$3=="25" { $1=$1".XY" } { print }' OFS='\t' "${outputFile}" > "${outputFile}.tmp"
\mv "${outputFile}".tmp "${outputFile}"

# position-based clumping of gene-based results
echo "Clumping gene-based results."
scriptDir=$(dirname "$0") # scriptDir=$(echo $(pwd)/code/genetics)
Rscript "${scriptDir}"/magma.clumping.R "${outputFile}" "${clumpingWindow}"

# map gene-level to snp-level results
echo "Mapping gene-level to snp-level results."
scriptDir=$(dirname "$0") # scriptDir=$(echo $(pwd)/code/genetics)
Rscript "${scriptDir}"/magma.mapping.R  "${conditionalFile}" "${outputFile}".clumped "${pthreshMapping}" "${clumpingWindow}" "${outputFile}".mapping

# clean up
echo "Cleaning up."
outDir=$(echo "${outputFile}" | sed 's!/[^/]*$!/!')
rm -rf "${outputFile}".dir
chmod 770 "${outputFile}"*

echo "--- Completed: MAGMA gene-based annotation --- "
