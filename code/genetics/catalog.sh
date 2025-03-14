#!/bin/bash

# ===============================
# === run GWAS catalog lookup ===
# ===============================

# get arguments
inFile="${1}" # inFile="results/gap_gm/gwama/eur/conditional/conditional.cleaned.tophits.annovar.txt"
outFile="${2}" # outFile="results/gap_gm/gwama/eur/catalog/catalog.txt"
cols2keep="${3}" # cols2keep="LOCUS_CNT,SNP_CNT,CHR,BP,ID,CYTOBAND,LEAD_SNP,LEAD_SNP_P,LEAD_SNP_KB,LEAD_SNP_RSQ,LEAD_SNP_ALLELES,REGION,NEAREST_GENE,NEAREST_GENE_DESCRIPTION,NEAREST_GENE_BIOTYPE,DISTANCE,A1,A2,A1_FREQ,BETA,SE,P,N"
snpCol="${4}" # snpCol="ID"
catalogFile="${5}" # catalogFile="data/gwas_catalog/gwas_catalog_v1.0-associations_e112_r2024-09-13.tsv"
catalogThresh="${6}" # catalogThresh="5E-8"
chrFileHandler="${7}" # chrFileHandler='/fast/software/pleiofdr/bfile/chr${i}' # this argument can be omitted when performing direct variant lookups only
rsq="${8}" # rsq="0.8" # this argument can be omitted when performing direct variant lookups only
kb="${9}" # kb="500" # this argument can be omitted when performing direct variant lookups only

# echo settings
echo $'\n'"--- run GWAS catalog lookup: settings ---"
echo "inFile: ${inFile}"
echo "outFile: ${outFile}"
echo "cols2keep: ${cols2keep}"
echo "snpCol: ${snpCol}"
echo "catalogFile: ${catalogFile}"
echo "catalogThresh: ${catalogThresh}"
echo "chrFileHandler: ${chrFileHandler}"
echo "rsq: ${rsq}"
echo "kb: ${kb}"$'\n'
sleep 5

# get list of variants in ld and add snps to list for catalog lookup
if [[ "${chrFileHandler}" != "" ]]; then
	snps=$(awk -v snpCol="${snpCol}" -F'\t' 'NR==1 { for(i=1;i<=NF;++i) { if($i==snpCol) { snpcolidx=i } }; next } { print $snpcolidx }' "${inFile}")
	> "${outFile}".ldall
	> "${outFile}".logall
	for snp in ${snps}; do
		i="$(cat $(i=""; eval echo "${chrFileHandler}"*.bim) | grep -w "${snp}" | awk 'NR==1 {print $1}')"
		if [[ "${i}" != "" ]]; then
			plink --bfile "$(eval echo ${chrFileHandler})" \
			--r2 dprime \
			--ld-snp "${snp}" \
			--ld-window-kb "${kb}" \
			--ld-window 99999 \
			--ld-window-r2 "${rsq}" \
			--out "${outFile}.${snp}"
		awk '(NR!=1 && FNR==1) || $3 == $6 { next } { print $1,$2,$3,$4,$5,$6,$7 }' OFS='\t' "${outFile}.${snp}.ld" >> "${outFile}".ldall
		awk '{ print }' OFS='\t' "${outFile}.${snp}.log" >> "${outFile}".logall
		fi
	done

	# merge .ld output files, remove lines where tag-snp is ld-snp, save tab-delimited
	

	# add tag-snps to input file
	header=$(echo "$(echo "${cols2keep}" | sed 's/,/\t/g')")
	awk -v header="${header}" -v cols="${cols2keep}" -v snpCol="${snpCol}" -F'\t' '
		BEGIN { ncols=split(cols,colnames,","); print header }
		NR>1 && NR==FNR { snpcount[$3]++; i=snpcount[$3]; snpinfo[$3"_"i]=$3" ("$6";"$7")"; next }
		FNR==1 { for(i=1;i<=NF;++i) { if($i==snpCol) { snpcolidx=i } } }
		FNR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } }; next }
		{ for(i=1;i<=ncols;++i) { if(i>1) { output = output"\t" };
								  		   if(colidx[i]==snpcolidx) { output = output""$snpcolidx }
								  		   else { output = output""$colidx[i] } }; 
		  print output; output = "" }
		$snpcolidx in snpcount { for(j=1;j<=snpcount[$snpcolidx];++j) { for(i=1;i<=ncols;++i)
								 	{ if(i>1) { output = output"\t" };
								  	  if(colidx[i]==snpcolidx) { output = output""snpinfo[$snpcolidx"_"j] }
								 	  else { output = output""$colidx[i] } };
									print output; output = "" } }
		' OFS="\t" "${outFile}.ldall" "${inFile}"   \
		> "${outFile}.tmp"
		inFile="${outFile}.tmp"
		echo ""
fi

# run catalog lookup
echo " - running GWAS catalog lookup."
header=$(echo "$(echo "${cols2keep}" | sed 's/,/\t/g')"$'\t'"catalog_pubmedid"$'\t'"catalog_author"$'\t'"catalog_date"$'\t'"catalog_trait"$'\t'"catalog_genes"$'\t'"catalog_a1"$'\t'"catalog_a1freq"$'\t'"catalog_orbeta"$'\t'"catalog_pval")
awk -v header="${header}" -v cols="${cols2keep}" -v snpCol="${snpCol}" -v catalogThresh="${catalogThresh}" -F'\t' '
	BEGIN { ncols=split(cols,colnames,","); print header }
	FNR==1 { file++ }
	NR==1 { for(i=1;i<=NF;++i) { if($i==snpCol) { snpcolidx=i } } }
	NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } }; next }
	NR==FNR { id = $snpcolidx; sub(".*\(","",id); sub("\;.*","",id); for(i=1;i<=ncols;++i) { length(output)==0 ? output = $colidx[i] : output = output"\t"$colidx[i] }; snp[id] = output; output=""; next}
	file==2 && $22 in snp && $28 < catalogThresh { length(snp2nd[$22])==0 ? snp2nd[$22]=snp[$22]"\t"$2"\t"$3"\t"$4"\t"$8"\t"$14"\t"$21"\t"$27"\t"$31"\t"$28 : snp2nd[$22]=snp2nd[$22]"\n"snp[$22]"\t"$2"\t"$3"\t"$4"\t"$8"\t"$14"\t"$21"\t"$27"\t"$31"\t"$28; next}
	file==3 { id = $snpcolidx; sub(".*\(","",id); sub("\;.*","",id); if(id in snp2nd) { print snp2nd[id] } }' OFS="\t" "${inFile}" "${catalogFile}" "${inFile}" \
	> "${outFile}"

# clean up
echo " - cleaning up."
if [[ "${chrFileHandler}" != "" ]]; then
	rm -f "${outFile}.tmp"
	rm -f "${outFile}"*.ld
	rm -f "${outFile}"*.log
	\mv "${outFile}.ldall" "$(echo "${outFile}" | sed 's/.txt$/.ld/g')"
	\mv "${outFile}.logall" "$(echo "${outFile}" | sed 's/.txt$/.log/g')"
	chmod 770 "$(echo "${outFile}" | sed 's/.txt$/.ld/g')"
	chmod 770 "$(echo "${outFile}" | sed 's/.txt$/.log/g')"
fi
chmod 770 "${outFile}"*

echo $'\n'"--- Completed: GWAS catalog lookup ---"
