#!/bin/bash

# ===============================
# === run GTEx eQTL screening ===
# ===============================

# get arguments
trait="$1" # trait="gap_gm" 
targetDir="$2" # targetDir="results/${trait}/eqtl"
conditionalFile="$3" # conditionalFile="results/${trait}/conditional/conditional.cleaned.tophits.annovar.txt"
gtexDir="$4" # gtexDir="data/gtex"
pthresh=$5 # pthresh=5E-8

# echo settings
echo $'\n'"--- Settings for GTEx eQTL screening ---"
echo "trait: "${trait}
echo "targetDir: "${targetDir}
echo "conditionalFile: "${conditionalFile}
echo "gtexDir: "${gtexDir}
echo "pthresh: "${pthresh}$'\n'

# create target directory
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# get tophits and convert rsid to variant_id_b37
echo "Converting rsid to variant_id_b37."
variants=$(awk -F'\t' -v pthresh=${pthresh} '$17 < pthresh && $31 < pthresh { print $5, $3"_"$4"_"$12"_"$11"_b37", $3"_"$4"_"$11"_"$12"_b37", $1, $2, $6, $9}' OFS='\t' ${conditionalFile})

# convert variant_id_b37 to variant_id_b38 (try to match by CHR_BP_A2_A1_b37, CHR_BP_A1_A2_b37, or rsid)
echo "Converting variant_id_b37 to variant_id_b38."
header=$(echo id$'\t'LOCUS_COUNT$'\t'SNP_COUNT$'\t'LEAD_SNP$'\t'LEAD_SNP_R2$'\t'matchedBy$'\t'"$(zcat "${gtexDir}/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz" | head -1)")
awk -F'\t' -v header="${header}" 'BEGIN { print header } 
	NR==FNR { id1[$1]=$1"\t"$4"\t"$5"\t"$6"\t"$7"\tchr_bp_a1_a2_b37";
			  id2[$2]=$1"\t"$4"\t"$5"\t"$6"\t"$7"\tchr_bp_a2_a1_b37";
			  id3[$3]=$1"\t"$4"\t"$5"\t"$6"\t"$7"\trsid"; next }
	$8 in id2 { print id2[$8], $0; next }
	$8 in id3 { print id3[$8], $0; next } 
	$7 in id1 { print id1[$7], $0; next }
	' OFS='\t' <(echo "${variants}") <(gzip -dc "${gtexDir}/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz") > "${targetDir}/eqtl.translation"

# find single-tissue variant-gene associations
echo "Finding single-tissue variant-gene associations."
tissueFiles="$(ls -l ${gtexDir}/GTEx_Analysis_v8_eQTL/*.signif_variant_gene_pairs.txt.gz | awk '{ print $NF }')"
tissue="$(echo "${tissueFiles}" | sed 's/[.].*//g' | sed 's%.*/%%g')"
tissueFiles=($tissueFiles)
tissue=($tissue)

mkdir -p "${targetDir}/eqtl.singleTissue"
(
for (( i=0; i<${#tissueFiles[@]}; i++ )); do (
	awk -F'\t' -v tissue="${tissue[$i]}" 'NR==1 { next } NR==FNR { variants[$7]=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6 } FNR == 1 { next }
		$1 in variants { print variants[$1], tissue, $0 }' OFS='\t' "${targetDir}/eqtl.translation" <(gzip -dc "${tissueFiles[$i]}") > "${targetDir}/eqtl.singleTissue/${tissue[$i]}.txt"
	) &
done
wait
)

	# combine singleTissue results into one file
	header=$(echo id$'\t'LOCUS_COUNT$'\t'SNP_COUNT$'\t'LEAD_SNP$'\t'LEAD_SNP_R2$'\t'matchedBy$'\t'tissue$'\t'"$(zcat "${tissueFiles[0]}" | head -1)")
	awk -F'\t' -v header="${header}" 'BEGIN { print header } { print }' OFS='\t' ${targetDir}/eqtl.singleTissue/*.txt > "${targetDir}/eqtl.singleTissue.txt"
	(head -n 1 "${targetDir}/eqtl.singleTissue.txt" && tail -n +2 "${targetDir}/eqtl.singleTissue.txt" | sort -k2,2g -k3,3g) > "${targetDir}/eqtl.singleTissue.txt.tmp"
	\mv "${targetDir}/eqtl.singleTissue.txt.tmp" "${targetDir}/eqtl.singleTissue.txt"

	# only keep one `variant-gene-tissue` association per locus (i.e., keep results of the variant that is in strongest ld with lead variant)
	awk -F'\t' 'NR==1 { print; next} !seen[$2$7$9]++' "${targetDir}/eqtl.singleTissue.txt" > "${targetDir}/eqtl.singleTissue.txt.tmp"
	\mv "${targetDir}/eqtl.singleTissue.txt.tmp" "${targetDir}/eqtl.singleTissue.txt"

# find multi-tissue variant-gene associations
echo "Finding multi-tissue variant-gene associations."
header=$(echo id$'\t'LOCUS_COUNT$'\t'SNP_COUNT$'\t'LEAD_SNP$'\t'LEAD_SNP_R2$'\t'matchedBy$'\t'gene_id$'\t'"$(zcat "${gtexDir}/GTEx_Analysis_v8.metasoft.txt.gz" | head -1)")
awk -F'\t' -v header="${header}" 'BEGIN { print header } NR==1 { next } NR==FNR { variants[$7]=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6 } FNR == 1 { next }
	{ id=$1; sub(/[,].*/,"",id); gene=$1; sub(/.*[,]/,"",gene) }
	id in variants { print variants[id], gene, $0 }' OFS='\t' "${targetDir}/eqtl.translation" <(gzip -dc "${gtexDir}/GTEx_Analysis_v8.metasoft.selected.txt.gz") > "${targetDir}/eqtl.multiTissue.txt"

	# sort by LOCUS_COUNT and SNP_COUNT
	(head -n 1 "${targetDir}/eqtl.multiTissue.txt" && tail -n +2 "${targetDir}/eqtl.multiTissue.txt" | sort -k2,2g -k3,3g) > "${targetDir}/eqtl.multiTissue.txt.tmp"
	\mv "${targetDir}/eqtl.multiTissue.txt.tmp" "${targetDir}/eqtl.multiTissue.txt"

	# only keep one `variant-gene` association per locus (i.e., keep results of the variant that is in strongest ld with lead variant)
	awk -F'\t' 'NR==1 { print; next} !seen[$2$7]++' "${targetDir}/eqtl.multiTissue.txt" > "${targetDir}/eqtl.multiTissue.txt.tmp"
	\mv "${targetDir}/eqtl.multiTissue.txt.tmp" "${targetDir}/eqtl.multiTissue.txt"

# convert ensembl gene id to hgnc id and symbol
scriptDir=$(dirname "$0")
Rscript ${scriptDir}/eqtl.ensemblTranslate.R "${targetDir}/eqtl.singleTissue.txt" "${targetDir}/eqtl.singleTissue.txt" "gene_id"
Rscript ${scriptDir}/eqtl.ensemblTranslate.R "${targetDir}/eqtl.multiTissue.txt" "${targetDir}/eqtl.multiTissue.txt" "gene_id"

# get summary
Rscript ${scriptDir}/eqtl.summarize.R "${targetDir}/eqtl.singleTissue.txt" "${targetDir}/eqtl.multiTissue.txt" "${targetDir}/eqtl.summary.txt"

# finish analysis
rm -rf ${targetDir}/eqtl.singleTissue
rm -f ${targetDir}/eqtl.translation
chmod -R 770 "${targetDir}"
echo "--- GTEx eQTL screening finished. ---"
