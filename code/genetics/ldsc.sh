#!/bin/bash

# ===========================
# === LD Score Regression ===
# ===========================

# get arguments
trait="${1}" # trait="pwr_cz_alpha"
targetDir="${2}" # targetDir="results/${trait}/ldsc/"
sumstats="${3}" # sumstats="results/${trait}/metal/metal.ivweight.gz"
sumstatsCols="${4}" # sumstatsCols="CHR,BP,KGP_RSID,A1,A2,BETA,SE,P,N" | columns to be used for ld score regression
sumstatsColNames="${5}" # sumstatsColNames="CHR,POS,SNP,A1,A2,BETA,SE,P,N" | column names for ld score regression (see ldsc presets)
LDsnplist="${6}" # LDsnplist="data/ldsc/resources/w_hm3.noMHC.snplist"
LDchr="${7}" # LDchr="data/ldsc/resources/eur_w_ld_chr/"
LDbaselineH2="${8}" # LDbaselineH2="data/ldsc/resources/1000G_Phase3_baselineLD_v2.2/baselineLD."
LDbaselineTau="${9}" # LDbaselineTau="data/ldsc/resources/1000G_Phase3_baseline_v1.2/baseline."
LDweights="${10}" # LDweights="data/ldsc/resources/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
LDfrqfiles="${11}" # LDfrqfiles="data/ldsc/resources/1000G_Phase3_frq/1000G.EUR.QC."
partitioned="${12}" # partitioned=0
celltypes="${13}" # celltypes="data/ldsc/resources/Multi_tissue_gene_expr_fullpaths.ldcts"
celltypegroups="${14}" # celltypegroups=$(echo /fast/software/ldsc/resources/1000G_Phase3_cell_type_groups/cell_type_group.{1..10}. | sed 's/ /,/g')
ctglabels="${15}" # ctglabels="adrenal.pancreas,cardiovascular,cns,connective.bone,gi,hematopoietic,kidney,liver,other,skeletalmuscle"

# echo settings
echo $'\n'"--- LD Score Regression Settings ---"
echo "trait: "${trait}
echo "targetDir: "${targetDir}
echo "sumstats: "${sumstats}
echo "sumstatsCols: "${sumstatsCols}
echo "sumstatsColNames: "${sumstatsColNames}
echo "LDsnplist: "${LDsnplist}
echo "LDchr: "${LDchr}
echo "LDbaselineH2: "${LDbaselineH2}
echo "LDbaselineTau: "${LDbaselineTau}
echo "LDweights: "${LDweights}
echo "LDfrqfiles: "${LDfrqfiles}
echo "partitioned: "${partitioned}
echo "celltypes: "${celltypes}
echo "celltypegroups: "${celltypegroups}
echo "ctglabels: "${ctglabels}$'\n'
sleep 5 

# set targetdir and make folder
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# prepare sumstats
echo "Preparing sumstats."
header=$(echo "${sumstatsColNames}" | sed 's%,%\t%g')
awk -v header="${header}" -v cols="${sumstatsCols}" '
	BEGIN { ncols=split(cols,colnames,","); print header }
    NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } } }
    NR>1 { output=$colidx[1]; for(i=2;i<=ncols;++i) { output=output"\t"$colidx[i] }; print output}
	' OFS="\t" <(gzip -dc "${sumstats}") > "${targetDir}/ldsc.${trait}.prep"

# munge sumstats
echo "Munging sumstats."
munge_sumstats.py \
--sumstats "${targetDir}/ldsc.${trait}.prep" \
--merge-alleles "${LDsnplist}" \
--out "${targetDir}/ldsc.${trait}"
chmod -R 770 "${targetDir}"

# ld score regression
echo "Running LD Score regression."
ldsc.py \
--h2 "${targetDir}/ldsc.${trait}.sumstats.gz" \
--ref-ld-chr "${LDchr}" \
--w-ld-chr "${LDchr}" \
--out "${targetDir}/ldsc.h2"

# get summary
snps="$(cat "${targetDir}/ldsc.h2.log" | grep -m1 "After merging with regression SNP LD, " | awk '{ print $7 }')"
h2="$(cat "${targetDir}/ldsc.h2.log" | grep -m1 "Total Observed scale h2: " | awk '{ print $5, $6 }')"
Lambda="$(cat "${targetDir}/ldsc.h2.log" | grep -m1 "Lambda GC: " | awk '{ print $3 }')"
Chi2="$(cat "${targetDir}/ldsc.h2.log" | grep -m1 "Mean Chi^2: " | awk '{ print $3 }')"
Intercept="$(cat "${targetDir}/ldsc.h2.log" | grep -m1 "Intercept: " | awk '{ print $2, $3 }')"
Ratio="$(cat "${targetDir}/ldsc.h2.log" | grep -m1 "Ratio: " | awk '{ print $2, $3 }')"
echo trait$'\t'snp_count$'\t'h2$'\t'lambda_GC$'\t'chi2$'\t'intercept$'\t'ratio > "${targetDir}/ldsc.h2.results"
echo ${trait}$'\t'${snps}$'\t'${h2}$'\t'${Lambda}$'\t'${Chi2}$'\t'${Intercept}$'\t'${Ratio} >> "${targetDir}/ldsc.h2.results"

# run partitioned ld score regression
if [[ ${partitioned} == 1 ]]; then
	echo "Running Partitioned LD Score regression."
	ldsc.py \
	--h2 "${targetDir}/ldsc.${trait}.sumstats.gz" \
	--ref-ld-chr "${LDbaselineH2}" \
	--w-ld-chr "${LDweights}" \
	--overlap-annot \
	--frqfile-chr "${LDfrqfiles}" \
	--print-coefficients \
	--out "${targetDir}/ldsc.partitioned"
fi

# run cell-type group analysis
if [[ ! -z ${celltypegroups} ]]; then
	echo "Running cell type group analysis."
	IFS=',' read -r -a celltypegroups <<< "${celltypegroups}"
	IFS=',' read -r -a ctglabels <<< "${ctglabels}"
	(for (( i=0; i<${#celltypegroups[@]}; i++ )); do (
		ldsc.py \
	    --h2 "${targetDir}/ldsc.${trait}.sumstats.gz" \
	    --ref-ld-chr ${celltypegroups[i]},"${LDbaselineH2}" \
	    --w-ld-chr "${LDweights}" \
	    --overlap-annot \
	    --frqfile-chr "${LDfrqfiles}" \
		--print-coefficients \
	    --out "${targetDir}/ldsc.ctg.${ctglabels[i]}"
	    ) &
	done
	wait)

	# summarise results on one file
	for (( i=0; i<${#celltypegroups[@]}; i++ )); do
		if [[ $i -eq 0 ]]; then head -1 "${targetDir}/ldsc.ctg.${ctglabels[i]}.results" > "${targetDir}/ldsc.ctg.results"; fi
		awk -F'\t' -v label=${ctglabels[i]} 'NR==2 { $1=label; print $0 }' OFS='\t' "${targetDir}/ldsc.ctg.${ctglabels[i]}.results" >> "${targetDir}/ldsc.ctg.results"
	done
fi

# run cell-type specific analysis
if [[ ! -z "${celltypes}" ]]; then
	echo "Running cell type specific analysis."
	ldsc.py \
    --h2-cts "${targetDir}/ldsc.${trait}.sumstats.gz" \
    --ref-ld-chr "${LDbaselineTau}" \
    --ref-ld-chr-cts "${celltypes}" \
    --w-ld-chr "${LDweights}" \
    --out "${targetDir}/ldsc.cts"
fi

# clean up
rm -f "${targetDir}/ldsc.${trait}.prep"
chmod -R 770 "${targetDir}/ldsc."*

