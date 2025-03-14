#!/bin/bash

# ============================
# === run GTEx eQTL lookup ===
# ============================

# get arguments
inFile="${1}" # inFile="results/pleiofdr/combined/pleio.crosstrait.txt"
outFile="${2}" # outFile="results/pleiofdr/combined/pleio.crosstrait.eqtl"
rsid="${3}" # rsid="snpid"
chr="${4}" # chr="chrnum"
bp="${5}" # bp="chrpos"
ref="${6}" # ref="a2"
alt="${7}" # alt="a1"
locus="${8}" # locus="discovnum"
assoc="${9}" # assoc="assocnum"
leadsnp="${10}" # leadsnp="leadsnp"
gtexDir="${11}" # gtexDir="data/gtex"

# echo settings
echo $'\n'"--- GTEx eQTL lookup: settings ---"
echo "inFile: ${inFile}"
echo "outFile: ${outFile}"
echo "rsid: ${rsid}"
echo "chr: ${chr}"
echo "bp: ${bp}"
echo "ref: ${ref}"
echo "alt: ${alt}"
echo "locus: ${locus}"
echo "assoc: ${assoc}"
echo "leadsnp: ${leadsnp}"
echo "gtexDir: ${gtexDir}"$'\n'
sleep 5

# get tophits and convert rsid to variant_id_b37
echo "Converting rsid to variant_id_b37."
variants=$(awk -F'\t' -v cols="${rsid},${chr},${bp},${ref},${alt},${locus},${assoc},${leadsnp}" 'BEGIN { ncols=split(cols,colnames,",") }
	NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } }; next }
{ print $colidx[1], $colidx[2]"_"$colidx[3]"_"$colidx[4]"_"$colidx[5]"_b37", $colidx[2]"_"$colidx[3]"_"$colidx[5]"_"$colidx[4]"_b37", $colidx[6], $colidx[7], $colidx[8]}' OFS='\t' "${inFile}")

# convert variant_id_b37 to variant_id_b38 (try to match by CHR_BP_A2_A1_b37, CHR_BP_A1_A2_b37, or rsid)
echo "Converting variant_id_b37 to variant_id_b38."
header=$(echo id$'\t'locusnum$'\t'assocnum$'\t'leadsnp$'\t'matchedBy$'\t'"$(zcat "${gtexDir}/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz" | head -1)")
awk -F'\t' -v header="${header}" 'BEGIN { print header } 
	NR==FNR { id1[$1]=$1"\t"$4"\t"$5"\t"$6"\trsid";
			  id2[$2]=$1"\t"$4"\t"$5"\t"$6"\tchr_bp_ref_alt_b37";
			  id3[$3]=$1"\t"$4"\t"$5"\t"$6"\tchr_bp_alt_ref_b37"; next }
	$8 in id2 { print id2[$8], $0; next }
	$8 in id3 { print id3[$8], $0; next } 
	$7 in id1 { print id1[$7], $0; next }
	' OFS='\t' <(echo "${variants}") <(gzip -dc "${gtexDir}"/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz) > "${outFile}.translation"

# find single-tissue variant-gene associations
echo "Finding single-tissue variant-gene associations."
tissueFiles="$(ls -l "${gtexDir}"/GTEx_Analysis_v8_eQTL/*.signif_variant_gene_pairs.txt.gz | awk '{ print $NF }')"
tissue="$(echo "${tissueFiles}" | sed 's/[.].*//g' | sed 's%.*/%%g')"
tissueFiles=($tissueFiles)
tissue=($tissue)

mkdir -p "${outFile}.singleTissue"
(for (( i=0; i<${#tissueFiles[@]}; i++ )); do (
	awk -F'\t' -v tissue="${tissue[$i]}" 'NR==1 { next } NR==FNR { variants[$6]=$1"\t"$2"\t"$3"\t"$4"\t"$5 } FNR == 1 { next }
		$1 in variants { print variants[$1], tissue, $0 }' OFS='\t' "${outFile}.translation" <(gzip -dc "${tissueFiles[$i]}") > "${outFile}.singleTissue/${tissue[$i]}.txt"
	) &
done
wait)

	# combine singleTissue results into one file
	header=$(echo id$'\t'locusnum$'\t'assocnum$'\t'leadsnp$'\t'matchedBy$'\t'tissue$'\t'"$(zcat "${tissueFiles[0]}" | head -1)")
	awk -F'\t' -v header="${header}" 'BEGIN { print header } { print }' OFS='\t' "${outFile}.singleTissue/"*.txt > "${outFile}.singleTissue.txt"
	(head -n 1 "${outFile}.singleTissue.txt" && tail -n +2 "${outFile}.singleTissue.txt" | sort -k2,2g -k3,3g) > "${outFile}.singleTissue.txt.tmp"
	\mv "${outFile}.singleTissue.txt.tmp" "${outFile}.singleTissue.txt"

	# only keep one `variant-gene-tissue` association per locus (i.e., keep results of the variant that is in strongest ld with lead variant)
	awk -F'\t' 'NR==1 { print; next} !seen[$2$6$8]++' "${outFile}.singleTissue.txt" > "${outFile}.singleTissue.txt.tmp"
	\mv "${outFile}.singleTissue.txt.tmp" "${outFile}.singleTissue.txt"

# find multi-tissue variant-gene associations
echo "Finding multi-tissue variant-gene associations."
header=$(echo id$'\t'locusnum$'\t'assocnum$'\t'leadsnp$'\t'matchedBy$'\t'gene_id$'\t'"$(zcat "${gtexDir}"/GTEx_Analysis_v8.metasoft.txt.gz | head -1)")
awk -F'\t' -v header="${header}" 'BEGIN { print header } NR==1 { next } NR==FNR { variants[$6]=$1"\t"$2"\t"$3"\t"$4"\t"$5 } FNR == 1 { next }
	{ id=$1; sub(/[,].*/,"",id); gene=$1; sub(/.*[,]/,"",gene) }
	id in variants { print variants[id], gene, $0 }' OFS='\t' "${outFile}.translation" <(gzip -dc "${gtexDir}"/GTEx_Analysis_v8.metasoft.selected.txt.gz) > "${outFile}.multiTissue.txt"

	# sort by locsnum and assocnum
	(head -n 1 "${outFile}.multiTissue.txt" && tail -n +2 "${outFile}.multiTissue.txt" | sort -k2,2g -k3,3g) > "${outFile}.multiTissue.txt.tmp"
	\mv "${outFile}.multiTissue.txt.tmp" "${outFile}.multiTissue.txt"

	# only keep one `variant-gene` association per locus (i.e., keep results of the variant that is in strongest ld with lead variant)
	awk -F'\t' 'NR==1 { print; next} !seen[$2$6]++' "${outFile}.multiTissue.txt" > "${outFile}.multiTissue.txt.tmp"
	\mv "${outFile}.multiTissue.txt.tmp" "${outFile}.multiTissue.txt"

# convert ensembl gene id to hgnc id and symbol
scriptDir=$(dirname "$0")
Rscript "${scriptDir}/eqtl.ensemblTranslate.R" "${outFile}.singleTissue.txt" "${outFile}.singleTissue.txt" "gene_id"
Rscript "${scriptDir}/eqtl.ensemblTranslate.R" "${outFile}.multiTissue.txt" "${outFile}.multiTissue.txt" "gene_id"

# get summary
Rscript "${scriptDir}/eqtl.summarize.R" "${outFile}.singleTissue.txt" "${outFile}.multiTissue.txt" "locusnum" "leadsnp" "${outFile}.summary.txt"

# finish analysis
rm -rf "${outFile}.singleTissue"
rm -f "${outFile}.translation"
chmod -R 770 "${outFile}"*
echo "--- Completed: GTEx eQTL lookup. ---"
