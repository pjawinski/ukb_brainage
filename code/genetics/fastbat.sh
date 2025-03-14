#!/bin/bash

# =========================================
# === run gene-based analysis (fastBAT) ===
# =========================================

# get arguments
targetDir="${1}" # targetDir="results/${trait}/fastbat"
chrFilehandler="${2}" # chrFilehandler='data/genetics/chr${i}/imp_mri_qc_EURjoined/bed/chr${i}_mri_qc'
sampleFile="${3}" # sampleFile="data/genetics/chr1/imp_mri_qc_EURjoined/imp_mri_qc.psam"
sumstats="${4}" # sumstats="results/${trait}/gwas/sumstats.txt.gz"
id="${5}" # id="ID"
p="${6}" # p="P"
conditionalFile="${7}" # conditionalFile="results/${trait}/conditional/conditional.cleaned.tophits.annovar.txt"
pthreshMapping=${8} # pthreshMapping=5E-8
glist_hg19="${9}" # glist_hg19="data/glist_hg19/glist-hg19-refseq.txt"
threads=${10} # threads=52
window=${11} # window=0
humandb="${12}" # humandb="data/annovar/humandb"
refseq="${13}" # refseq="data/refseq/GRCh37_latest_genomic.edit.gff.gz"
clumpingWindow=${14} # clumpingWindow=10000

# echo settings
echo $'\n'"--- Gene-based Analysis Settings (fastBAT) ---"
echo "targetDir: ${targetDir}"
echo "chrFilehandler: ${chrFilehandler}"
echo "sampleFile: ${sampleFile}"
echo "sumstats: ${sumstats}"
echo "id: ${id}"
echo "p: ${p}"
echo "conditionalFile: ${conditionalFile}"
echo "pthreshMapping: ${pthreshMapping}"
echo "glist_hg19: ${glist_hg19}"
echo "threads: ${threads}"
echo "window: ${window}"
echo "humandb: ${humandb}"
echo "refseq: ${refseq}"
echo "clumpingWindow: ${clumpingWindow}"$'\n'

# set targetdir and make folder
rm -rf "${targetDir}"
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# get reference sample (subjects with available trait data)
awk '{ print $1, $2 }' OFS="\t" "${sampleFile}" > "${targetDir}"/fastbat.sample.txt

# get SNPs and P-values from summary statistics
awk -v cols="${id},${p}" '
  BEGIN { ncols=split(cols,colnames,","); print "SNP\tP" } 
  NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } } }
  NR > 1 { print $colidx[1], $colidx[2] }' OFS="\t" <(gzip -dc "${sumstats}") > "${targetDir}"/fastbat.assoc.txt

# # For gcta software: change chromsomes to numeric values
# \cp ${geneticsDir}/chrX/imp_mri_qc/bed/chrX_mri_qc_numeric.bim ${geneticsDir}/chrX/imp_mri_qc/bed/chrX_mri_qc.bim 
# \cp ${geneticsDir}/chrY/imp_mri_qc/bed/chrY_mri_qc_numeric.bim ${geneticsDir}/chrY/imp_mri_qc/bed/chrY_mri_qc.bim
# \cp ${geneticsDir}/chrXY/imp_mri_qc/bed/chrXY_mri_qc_numeric_23.bim ${geneticsDir}/chrXY/imp_mri_qc/bed/chrXY_mri_qc.bim
# \cp ${geneticsDir}/chrMT/imp_mri_qc/bed/chrMT_mri_qc_numeric.bim ${geneticsDir}/chrMT/imp_mri_qc/bed/chrMT_mri_qc.bim

# run gene-based analysis (http://cnsgenomics.com/software/gcta/#Gene-basedtest)
echo "Running gene-based analysis."
mkdir -p "${targetDir}"/fastbat.log
chrlist=$(echo {1..22} X Y XY MT); chrlist=($chrlist)
chrCount=$(echo "${#chrlist[@]}")
threadsPerAnalysis=$(expr "${threads}" / "${chrCount}")
threadsPerAnalysis=$(( threadsPerAnalysis > 1 ? threadsPerAnalysis : 1 ))

for i in ${chrlist[@]}; do (
gcta64 \
--bfile "$(eval echo "${chrFilehandler}")" \
--keep "${targetDir}"/fastbat.sample.txt \
--fastBAT "${targetDir}"/fastbat.assoc.txt \
--fastBAT-gene-list "${glist_hg19}" \
--out "${targetDir}"/fastbat.log/fastbat.chr"${i}"."${window}"kb \
--fastBAT-wind "${window}" \
--thread-num "${threadsPerAnalysis}" \
--fastBAT-write-snpset
) &
done
wait

# # replace .bim files with originals
# \cp ${geneticsDir}/chrX/imp_mri_qc/bed/chrX_mri_qc_original.bim ${geneticsDir}/chrX/imp_mri_qc/bed/chrX_mri_qc.bim
# \cp ${geneticsDir}/chrY/imp_mri_qc/bed/chrY_mri_qc_original.bim ${geneticsDir}/chrY/imp_mri_qc/bed/chrY_mri_qc.bim
# \cp ${geneticsDir}/chrXY/imp_mri_qc/bed/chrXY_mri_qc_original.bim ${geneticsDir}/chrXY/imp_mri_qc/bed/chrXY_mri_qc.bim
# \cp ${geneticsDir}/chrMT/imp_mri_qc/bed/chrMT_mri_qc_original.bim ${geneticsDir}/chrMT/imp_mri_qc/bed/chrMT_mri_qc.bim

# concatenate gene.fastbat files
echo "Concatenating gene.fastbat files."
awk 'NR==1 { print; next } { $2="25"; print }' OFS='\t' "${targetDir}"/fastbat.log/fastbat.chrXY.0kb.gene.fastbat > "${targetDir}"/fastbat.log/fastbat.chrXY.0kb.gene.fastbat.tmp
\mv "${targetDir}"/fastbat.log/fastbat.chrXY.0kb.gene.fastbat.tmp "${targetDir}"/fastbat.log/fastbat.chrXY.0kb.gene.fastbat
awk 'NR==1 { print; next } FNR>1' "${targetDir}"/fastbat.log/fastbat.chr{{1..22},X,Y,XY,MT}.*.gene.fastbat > "${targetDir}"/fastbat."${window}"kb.gene.fastbat

# sort results by p value
sort -g -k 9 "${targetDir}"/fastbat."${window}"kb.gene.fastbat > "${targetDir}"/fastbat."${window}"kb.gene.fastbat.tmp
\mv "${targetDir}"/fastbat."${window}"kb.gene.fastbat.tmp "${targetDir}"/fastbat."${window}"kb.gene.fastbat

# make annotation
echo "Running annotation."
mkdir -p "${targetDir}"/annovar
header="CHR"$'\t'"START"$'\t'"END"$'\t'"REF"$'\t'"ALT"$'\t'"Gene"$'\t'"Chr"$'\t'"Start"$'\t'"End"$'\t'"No.SNPs"$'\t'"SNP_start"$'\t'"SNP_end"$'\t'"Chisq(Obs)"$'\t'"Pvalue"$'\t'"TopSNP.Pvalue"$'\t'"TopSNP"
awk -v header="$header" 'BEGIN { print header } FNR==1 { next } { CHR=$2 } CHR=="23" { CHR="X" } CHR=="24" { CHR="Y" } CHR=="25" { CHR="X" } CHR=="26" { CHR="MT" } { START=int($3+($4-$3)/2); print CHR, START, START, "A", "G", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11 }' OFS='\t' "${targetDir}"/fastbat."${window}"kb.gene.fastbat > "${targetDir}"/annovar/fastbat."${window}"kb.gene.fastbat4annovar

mkdir -p "${targetDir}"/annovar/cytoBand
annotate_variation.pl -regionanno -dbtype cytoBand -buildver hg19 -out "${targetDir}"/annovar/cytoBand/fastbat."${window}"kb.gene.fastbat.annovar "${targetDir}"/annovar/fastbat."${window}"kb.gene.fastbat4annovar "${humandb}"

# add cytoband
header="Gene"$'\t'"Chr"$'\t'"Start"$'\t'"End"$'\t'"No.SNPs"$'\t'"SNP_start"$'\t'"SNP_end"$'\t'"Chisq(Obs)"$'\t'"Pvalue"$'\t'"TopSNP.Pvalue"$'\t'"TopSNP"$'\t'"CYTOBAND"
awk -v header="$header" 'BEGIN { print header } NR==FNR { cyto[$8]=$2; next } FNR == 1 { next } !($1 in cyto) { print $0, "NA"; next} { print $0, cyto[$1] }' OFS='\t' "${targetDir}"/annovar/cytoBand/fastbat."${window}"kb.gene.fastbat.annovar.hg19_cytoBand "${targetDir}"/fastbat."${window}"kb.gene.fastbat > "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar

# add biotype and description
header="Gene"$'\t'"Chr"$'\t'"Start"$'\t'"End"$'\t'"No.SNPs"$'\t'"SNP_start"$'\t'"SNP_end"$'\t'"Chisq(Obs)"$'\t'"Pvalue"$'\t'"TopSNP.Pvalue"$'\t'"TopSNP"$'\t'"CYTOBAND"$'\t'"GENE_DESCRIPTION"$'\t'"GENE_BIOTYPE"
awk -F'\t' -v header="$header" 'BEGIN { print header } NR==FNR { gene[$10]=$12"\t"$11; next } FNR==1 { next } $1 in gene { print $0, gene[$1]; next} { print $0, "NA", "NA"; next }' OFS='\t' <(gzip -dc "${refseq}") "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar > "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar.tmp.txt
\mv "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar.tmp.txt "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar

  # no match? get refseq rows that contain id as synonym
  #echo "(2/4) Finding matches in RefSeq synonyms..."
  missings=$(awk -F'\t' '$NF == "NA" { print $1}' "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar)

  if [ ! -z "$missings" ]; then

    awk -F'\t|,' 'NR==FNR { missings[$1]; next}
      { for(k=1; k<=NF; k++) {if($k !~ /^[0-9]+$/ && $k in missings) { print FNR, $k; delete missings[$k] } } }
      ' OFS='\t' <(echo "$missings") <(gzip -dc "${refseq}") > "${targetDir}"/annovar/synonyms.rownum.txt

    # for identified synonyms get refseq id, biotype, and description
    awk -F'\t' 'NR==FNR { id[$1]=$2; next } FNR in id { print id[FNR], $10, $11, $12, $13}' OFS='\t' "${targetDir}"/annovar/synonyms.rownum.txt <(gzip -dc "${refseq}") > "${targetDir}"/annovar/synonyms.annotated.txt

    # merge with annotation list
    header="Gene"$'\t'"Chr"$'\t'"Start"$'\t'"End"$'\t'"No.SNPs"$'\t'"SNP_start"$'\t'"SNP_end"$'\t'"Chisq(Obs)"$'\t'"Pvalue"$'\t'"TopSNP.Pvalue"$'\t'"TopSNP"$'\t'"CYTOBAND"$'\t'"GENE_DESCRIPTION"$'\t'"GENE_BIOTYPE"
    awk -F'\t' -v header="$header" 'BEGIN { print header }
      NR==FNR { id[$1]=$3"\t"$4; next }
      FNR==1 { next } 
      $1 in id { output=$1; { for(i=2; i<=NF-2; i++) { output=output"\t"$i } }; print output, id[$1]; next }
      { print $0}' OFS='\t' "${targetDir}"/annovar/synonyms.annotated.txt "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar > "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar.tmp.txt
    \mv "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar.tmp.txt "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar

  fi

# add ENTREZ
echo "Adding ENTREZ identifier."
header=$(awk 'NR==1 { print }' "${targetDir}/fastbat.${window}kb.gene.fastbat.annovar")$'\t'"ENTREZ_ID"
awk -F'\t|,|;' -v header="$header" 'BEGIN { print header }
  NR==FNR { gsub(/Dbxref=/,""); for(k=1; k<=NF; k++) { if($k ~ /^GeneID:/) { gene[$10]=$k } }; next }
  FNR==1 { next } $1 in gene { print $0, gene[$1]; next }
  { print $0, "NA"}' OFS='\t' <(gzip -dc "${refseq}") "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar > "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar.tmp
\mv "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar.tmp "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar

# add HGNC
echo "Adding HGNC identifier."
header=$(awk 'NR==1 { print }' "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar)$'\t'"HGNC_ID"
awk -F'\t|,|;' -v header="$header" 'BEGIN { print header }
  NR==FNR { gsub(/Dbxref=/,""); gsub(/HGNC:HGNC:/,"HGNC:"); for(k=1; k<=NF; k++) { if($k ~ /^HGNC:/) { gene[$10]=$k } }; next }
  FNR==1 { next } $1 in gene { print $0, gene[$1]; next }
  { print $0, "NA"}' OFS='\t' <(gzip -dc "${refseq}") "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar > "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar.tmp
\mv "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar.tmp "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar

# add ".XY" to pseudoautosomal gene-based results
echo "Adding .XY to pseudoautosomal genes."
awk -F'\t' '$2=="25" { $1=$1".XY" } { print }' OFS='\t' "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar > "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar.tmp.txt
\mv "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar.tmp.txt "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar

# position-based clumping of gene-based results
echo "Clumping gene-based results."
scriptDir=$(dirname "$0")
Rscript "${scriptDir}"/fastbat.clumping.R "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar "${clumpingWindow}"

# map gene-level to snp-level results
echo "Mapping gene-level to snp-level results."
scriptDir=$(dirname "$0")
Rscript "${scriptDir}"/fastbat.mapping.R  "${conditionalFile}" "${targetDir}"/fastbat."${window}"kb.gene.fastbat.annovar.clumped "${pthreshMapping}" "${clumpingWindow}" "${targetDir}"/fastbat."${window}"kb.mapping

# clean up
echo "Cleaning up."
rm -f "${targetDir}"/fastbat.assoc.txt
rm -f "${targetDir}"/fastbat.sample.txt
mv "${targetDir}"/annovar "${targetDir}"/fastbat.annovar
tar -cvzf "${targetDir}"/fastbat.annovar.tar.gz --directory="${targetDir}" fastbat.annovar --remove-files
tar -cvzf $"${targetDir}"/fastbat.log.tar.gz --directory="${targetDir}" fastbat.log --remove-files
chmod -R 770 "${targetDir}"/*

echo "--- Gene-based analysis finished. --- "
