#!/bin/bash

# =========================================
# === Create regional association plots ===
# =========================================

# get arguments
trait="${1}" # trait="gap_wm"
traitDescription="${2}" # traitDescription="grey matter"
targetDir="${3}" # targetDir="results/gap_gm/locuszoom/"
chrFilehandler="${4}" # chrFilehandler='data/genetics/chr${i}/imp_mri_qc_EURjoined/vcf/chr${i}_mri_qc.vcf.gz'
sumstats="${5}" # sumstats="results/gap_gm/gwas/sumstats.txt.gz"
id="${6}" # id="ID"
p="${7}" # p="P"
conditionalFile="${8}" # conditionalFile="results/gap_gm/conditional/conditional.cleaned.tophits.annovar.txt"
pthresh=${9} # pthresh=5E-8
flank="${10}" # flank="500kb"
snpsManual="${11}" # snpsManual="rs796226228"
flankManual="${12}" # flankManual="1000kb"
chrManual="${13}" # chrManual="1"
startManual="${14}" # startManual="214800000"
endManual="${15}" # endManual="216200000"

# echo settings
echo $'\n'"--- Locuszoom Settings ---"
echo "trait: ${trait}"
echo "traitDescription: ${traitDescription}"
echo "targetDir: ${targetDir}"
echo "chrFilehandler: ${chrFilehandler}"
echo "sumstats: ${sumstats}"
echo "id: ${id}"
echo "p: ${p}"
echo "conditionalFile: ${conditionalFile}"
echo "pthresh: ${pthresh}"
echo "flank: ${flank}"
echo "snpsManual: ${snpsManual}"
echo "flankManual: ${flankManual}"
echo "chrManual: ${chrManual}"
echo "startManual: ${startManual}"
echo "endManual: ${endManual}"$'\n'
if [ "${startManual}" != "" ]; then
	echo "Note: startManual has been set, so flankManual will be ignored."$'\n'
fi
sleep 5

# set targetdir and make folder
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# Get Markername and P-value
echo "Getting Markername and P-value."
awk -v cols="${id},${p}" '
  BEGIN { ncols=split(cols,colnames,","); print "MarkerName\tP-value" } 
  NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } } }
  NR > 1 { print $colidx[1], $colidx[2] }' OFS="\t" <(gzip -dc "${sumstats}") > "${targetDir}"/lz.metal.txt

# get list of index snps
echo "Getting index snps and respective chromosomes."
snps=$(awk -F'\t' -v pthresh="${pthresh}" -v cols="ID,LEAD_SNP,CHR,P,LEAD_SNP_pJ" 'BEGIN { ncols=split(cols,colnames,",") } NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } }; next } $colidx[1] == $colidx[2] && $colidx[4] < pthresh && $colidx[5] < pthresh { print $colidx[1] }' "${conditionalFile}"); snps=($snps)
chr=$(awk -F'\t' -v pthresh="${pthresh}" -v cols="ID,LEAD_SNP,CHR,P,LEAD_SNP_pJ" 'BEGIN { ncols=split(cols,colnames,",") } NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } }; next } $colidx[1] == $colidx[2] && $colidx[4] < pthresh && $colidx[5] < pthresh { print $colidx[3] }' "${conditionalFile}"); chr=($chr)

# set flanking region for each snp
echo "Setting flanking regions."
flnk=$(yes "${flank}" | head -${#snps[@]}); flnk=(${flnk})

# adding/replacing settings of manually set SNPs
if [ "${startManual}" != "" ]; then
	echo "Adding start and end positions of manually set SNPs."
	snpsManual=(${snpsManual})
	chrManual=(${chrManual})
	startManual=(${startManual})
	endManual=(${endManual})

	# remove snpsManual from snps
	for i in "${!snpsManual[@]}"; do
		for j in "${!snps[@]}"; do
			if [[ "${snpsManual[$i]}" = "${snps[$j]}" ]]; then
				unset snps["${j}"]
				unset chr["${j}"]
				unset flnk["${j}"]
			fi
		done
	done

	snps=(${snps[@]})
	chr=(${chr[@]})
	flnk=(${flnk[@]})

elif [ "${flankManual}" != "" ]; then
	echo "Adding/replacing flanking regions of manually set SNPs."
	snpsManual=(${snpsManual})
	chrManual=(${chrManual})
	flankManual=(${flankManual})
	unset startManual

	# replace / add snpManual in snps
	for i in "${!snpsManual[@]}"; do
		for j in "${!snps[@]}"; do
			if [[ "${snpsManual[$i]}" = "${snps[$j]}" ]]; then
				flnk[$j]=${flankManual[$i]}
				replaced=TRUE
			fi
		done
		if [ "${replaced}" != "TRUE" ]; then
			snps+=(${snpsManual[$i]})
			chr+=(${chrManual[$i]})
			flnk+=(${flankManual[$i]})
		fi
		unset replaced
	done
fi

# draw locuszoom plot
echo "Drawing locuszoom plots."
rm -rf "${targetDir}"/lz.output/
mkdir -p "${targetDir}"/lz.output/

for (( i=0; i<${#snps[@]}; i++ )); do (
	locuszoom \
	--metal "${targetDir}"/lz.metal.txt \
	--refsnp "${snps[$i]}" \
	--ignore-vcf-filter \
	--flank "${flnk[$i]}" \
	--build hg19 \
	--ld-vcf $(eval echo "$(echo "${chrFilehandler}" | sed 's/${i}/${chr[$i]}/g')") \
	--no-date \
	--prefix "${targetDir}"/lz.output/lz.chr"${chr[$i]}" \
	--plotonly \
	width=8 \
	height=7
) &
done

if [ "${startManual}" != "" ]; then
	for (( i=0; i<${#snpsManual[@]}; i++ )); do (
		locuszoom \
		--metal "${targetDir}"/lz.metal.txt \
		--refsnp "${snpsManual[$i]}" \
		--chr "${chrManual[$i]}" \
		--start "${startManual[$i]}" \
		--end "${endManual[$i]}" \
		--ignore-vcf-filter \
		--build hg19 \
		--ld-vcf $(eval echo "$(echo "${chrFilehandler}" | sed 's/${i}/${chrManual[$i]}/g')") \
		--no-date \
		--prefix "${targetDir}"/lz.output/lz.chr"${chrManual[$i]}" \
		--plotonly \
		width=8 \
		height=7
	) &
	done
fi
wait

# draw six locuszoom plots on a single page (2x3 plots)
scriptDir=$(dirname "$0")
Rscript "${scriptDir}"/locuszoom.plot.R "${traitDescription}" "${targetDir}" "${conditionalFile}"

# clean up
echo "Cleaning up."
rm -f "${targetDir}"/lz.metal.txt
echo "--- Completed: Locuszoom. --- "
