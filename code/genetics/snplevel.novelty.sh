#!/bin/bash

# ==========================================================================
# === identify novel discoveries by comparison with previous discoveries ===
# ==========================================================================

# get arguments
ownDiscoveries="$1" # ownDiscoveries="results/combined/snplevel.gws.txt"
prevDiscoveries="$2" # prevDiscoveries="data/prevDiscoveries/Jonsson_2019_munged.txt data/prevDiscoveries/Kaufmann_2019_munged.txt data/prevDiscoveries/Ning_2020_munged.txt data/prevDiscoveries/Smith_2020_munged.txt data/prevDiscoveries/Leonardsen_2023_munged.txt"
targetDir="$3" # targetDir="results/combined"
geneticsDir="$4" # geneticsDir="data/genetics"
LDsample="$5" # LDsample="data/gap_gm/gap_gm.txt"

# echo settings
echo $'\n'"--- identify novel discoveries ---"
echo "ownDiscoveries: "${ownDiscoveries}
echo "prevDiscoveries: "${prevDiscoveries}
echo "targetDir: "${targetDir}
echo "geneticsDir: "${geneticsDir}
echo "LDsample: "${LDsample}

# make targetDir
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# merge prevDiscoveries
awk '{ print }' ${prevDiscoveries} > ${targetDir}/snplevel.novelty.prevDiscoveries.txt

# identify duplicate (multiallelic) variants in .pvar files and remove them
variants=$(awk -F'\t' '{ print $3 }' ${targetDir}/snplevel.novelty.prevDiscoveries.txt | sort -u)
rm -f ${targetDir}/snplevel.novelty.prevDiscoveries.multiallelic.txt
for i in {1..22} X Y XY MT; do (
awk 'NR==FNR { snp[$1]; next } $3 in snp { print $3 }' OFS='\t' <(echo "$variants") ${geneticsDir}/chr${i}/imp_mri/chr${i}_mri.pvar | sort | uniq -c | awk '$1 > 1 { print $2}' >> ${targetDir}/snplevel.novelty.prevDiscoveries.multiallelic.txt
) &
done
wait
awk 'NR==FNR { multi[$1]; next } !($3 in multi) { print }' ${targetDir}/snplevel.novelty.prevDiscoveries.multiallelic.txt ${targetDir}/snplevel.novelty.prevDiscoveries.txt > ${targetDir}/snplevel.novelty.prevDiscoveries.cleaned.txt

# munge ownDiscoveries
awk -F'\t' 'NR==1 { next } { print "ownDiscoveries", $3, $6, $13 }' OFS='\t' ${ownDiscoveries} > ${targetDir}/snplevel.novelty.ownDiscoveries.txt

# merge together
awk '{ print }' ${targetDir}/snplevel.novelty.ownDiscoveries.txt ${targetDir}/snplevel.novelty.prevDiscoveries.cleaned.txt > ${targetDir}/snplevel.novelty.sumstats.txt

# sort by ID and p-value and remove redundant SNPs
header="STUDY"$'\t'"ID"$'\t'"P"
sort -t$'\t' -k3,3 -k4,4g ${targetDir}/snplevel.novelty.sumstats.txt | awk -F'\t' -v header="${header}" 'BEGIN { print header } !seen[$3]++ { print $1, $3, $4}' OFS='\t' > ${targetDir}/snplevel.novelty.sumstats4clumping.txt

# get participant IDs of LD reference sample
awk '{ print $1, $2 }' "${LDsample}" > "${targetDir}/snplevel.novelty.subs.txt"

# get snplist 
awk 'NR > 1 { print $2 }' ${targetDir}/snplevel.novelty.sumstats4clumping.txt > ${targetDir}/snplevel.novelty.sumstats4clumping.snplist

# do clumping
mkdir -p "${targetDir}/snplevel.novelty.clumping"
for i in {1..22} X XY Y MT; do (
	plink2 --pfile "${geneticsDir}/chr${i}/imp_mri/chr${i}_mri" \
	--keep "${targetDir}/snplevel.novelty.subs.txt" \
	--extract "${targetDir}/snplevel.novelty.sumstats4clumping.snplist" \
	--make-bed \
	--out "${targetDir}/snplevel.novelty.clumping/chr${i}_mri"

	LD_src="${targetDir}/snplevel.novelty.clumping/chr${i}_mri"
	plink --bfile "$LD_src" \
	--keep "${targetDir}/snplevel.novelty.subs.txt" \
	--clump "${targetDir}/snplevel.novelty.sumstats4clumping.txt" \
	--clump-snp-field "ID" \
	--clump-field "P" \
	--clump-p1 1.000 \
	--clump-p2 1.000 \
	--clump-r2 0.1 \
	--clump-kb 10000 \
	--clump-verbose \
	--out "$targetDir/snplevel.novelty.clumping/chr${i}"
	) &
done
wait

# Merge clumped files
echo $'\n'"Merging clumped results."
cat ${targetDir}/snplevel.novelty.clumping/chr*.clumped > ${targetDir}/snplevel.novelty.clumping/clumped
awk '$11 ~ /^[0-9]+$/ { indexsnp = $3; print $3, $3, 0, 1, "A1/A2"; next } $1=="" || $1=="KB" || $1=="CHR" || $1=="(INDEX)" || $2 == "not" || $1 == "RANGE:" || $1 == "SPAN:" || $1 ~ /-/ { next } {print indexsnp, $1, $2, $3, $4}' OFS='\t' ${targetDir}/snplevel.novelty.clumping/clumped > ${targetDir}/snplevel.novelty.clumping/clumped.tmp
\mv ${targetDir}/snplevel.novelty.clumping/clumped.tmp ${targetDir}/snplevel.novelty.clumping/clumped

# merge clumping results and initial file
awk -F'\t' 'NR==FNR { snp[$2]=$1"\t"$3"\t"$4"\t"$5; next} $3 in snp { print $0, snp[$3] }' OFS="\t" ${targetDir}/snplevel.novelty.clumping/clumped ${targetDir}/snplevel.novelty.sumstats.txt > ${targetDir}/snplevel.novelty.clumping/clumped.tmp
\mv ${targetDir}/snplevel.novelty.clumping/clumped.tmp ${targetDir}/snplevel.novelty.clumping/clumped

# get chr:pos of variants
echo "Adding chr and bp of variants."
variants=$(awk -F'\t' '{ print $3 }' ${targetDir}/snplevel.novelty.clumping/clumped | sort -u)
rm -f ${targetDir}/snplevel.novelty.variants_chr_pos.txt
for i in {1..22} X Y XY MT; do (
awk 'NR==FNR { snp[$3]=$1"\t"$2; next } $1 in snp { print $1, snp[$1] }' OFS='\t' ${geneticsDir}/chr${i}/imp_mri/chr${i}_mri.pvar <(echo "$variants") >> ${targetDir}/snplevel.novelty.variants_chr_pos.txt
) &
done
wait

# add chr:pos 
awk -F'\t' 'NR==FNR { snp[$1]=$2"\t"$3; next } $3 in snp { print $0, snp[$3] }' OFS='\t' ${targetDir}/snplevel.novelty.variants_chr_pos.txt ${targetDir}/snplevel.novelty.clumping/clumped > ${targetDir}/snplevel.novelty.clumping/clumped.tmp
\mv ${targetDir}/snplevel.novelty.clumping/clumped.tmp ${targetDir}/snplevel.novelty.clumping/clumped

# add LEAD_SNP chr:pos 
awk -F'\t' 'NR==FNR { snp[$1]=$2"\t"$3; next } $5 in snp { print $0, snp[$5] }' OFS='\t' ${targetDir}/snplevel.novelty.variants_chr_pos.txt ${targetDir}/snplevel.novelty.clumping/clumped > ${targetDir}/snplevel.novelty.clumping/clumped.tmp
\mv ${targetDir}/snplevel.novelty.clumping/clumped.tmp ${targetDir}/snplevel.novelty.clumping/clumped

# sort by LEAD_SNP chr:pos and LD
awk -F'\t' '$11=="X" { $11=23 } $11=="Y" { $11=24 } $11=="XY" { $11=25 } $11=="MT" { $11=26 } { print }' OFS='\t' ${targetDir}/snplevel.novelty.clumping/clumped | sort -t$'\t' -k11,11g -k12,12g -k7,7gr -k4,4g | awk -F'\t' '$11==23 { $11="X" } $11==24 { $11="Y" } $11==25 { $11="XY" } $11==26 { $11="MT" } { print }' OFS='\t' > ${targetDir}/snplevel.novelty.clumping/clumped.tmp
\mv ${targetDir}/snplevel.novelty.clumping/clumped.tmp ${targetDir}/snplevel.novelty.clumping/clumped

# add LOCUS_COUNT and DISCOVERY_COUNT
header=LOCUS_COUNT$'\t'DISCOV_COUNT$'\t'STUDY$'\t'TRAIT$'\t'ID$'\t'P$'\t'LEAD_SNP$'\t'LEAD_SNP_KB$'\t'LEAD_SNP_RSQ$'\t'LEAD_SNP_ALLELES$'\t'CHR$'\t'BP$'\t'LEAD_SNP_CHR$'\t'LEAD_SNP_BP
awk -F'\t' -v header="$header" 'BEGIN { print header; leadsnp="null"; locus_cnt = 0; snp_cnt = 1 } leadsnp != $5 { leadsnp = $5; locus_cnt++; snp_cnt = 1; print locus_cnt, snp_cnt, $0; next } { snp_cnt++; print locus_cnt, snp_cnt, $0; next }' OFS='\t' ${targetDir}/snplevel.novelty.clumping/clumped > ${targetDir}/snplevel.novelty.clumping/clumped.tmp
\mv ${targetDir}/snplevel.novelty.clumping/clumped.tmp ${targetDir}/snplevel.novelty.clumping/clumped

# test if any multiallelic variant is near ownDiscoveries
if [ -e ${targetDir}/snplevel.novelty.prevDiscoveries.multiallelic.txt ]; then
echo "Clumping multiallelic variants by position."

# 1) get chr:pos of variants
	variants=$(awk -F'\t' '{ print $1 }' ${targetDir}/snplevel.novelty.prevDiscoveries.multiallelic.txt | sort -u)
	rm -f ${targetDir}/snplevel.novelty.prevDiscoveries.multiallelic_chr_pos.txt
	for i in {1..22} X Y XY MT; do (
	awk 'NR==FNR { snp[$3]=$1"\t"$2; next } $1 in snp { print $1, snp[$1] }' OFS='\t' ${geneticsDir}/chr${i}/imp_mri/chr${i}_mri.pvar <(echo "$variants") >> ${targetDir}/snplevel.novelty.prevDiscoveries.multiallelic_chr_pos.txt
	) &
	done
	wait

# 2) determine smallest distance to own discoveries
	header=$(echo ID_MULTI$'\t'CHR_MULTI$'\t'BP_MULTI$'\t'DISTANCE_MULTI$'\t'"$(head -1 ${targetDir}/snplevel.novelty.clumping/clumped)")
	awk -F'\t' -v header="$header" 'BEGIN { print header }
				NR==FNR {chr[$1]=$2; pos[$1]=$3; dist[$1]=999999999; info[$1]="NA"; next }
				FNR > 1 && $3 == "ownDiscoveries" { for (key in chr) if ($11 == chr[key] && sqrt((pos[key]-$12)^2) < dist[key]) { dist[key] = sqrt((pos[key]-$12)^2); info[key]=$0 } }
				END { for (key in chr) { print key, chr[key], pos[key], dist[key], info[key] } }' OFS='\t' ${targetDir}/snplevel.novelty.prevDiscoveries.multiallelic_chr_pos.txt ${targetDir}/snplevel.novelty.clumping/clumped > ${targetDir}/snplevel.novelty.prevDiscoveries.multiallelic.clumped.txt

# 3) add to locus if distance is smaller than 10 MB
	awk -F'\t' 'BEGIN { locus=1 } 
				FNR == 1 { filenum++ }
				filenum == 1 && NR > 1 { id[$1]=$1; locus_count[$1]=$5; lead_snp[$1]=$9; lead_snp_kb[$1]=$4/1000; lead_snp_rsq[$1]="NA"; lead_snp_alleles[$1]="NA"; chr[$1]=$2; bp[$1]=$3; lead_snp_chr[$1]=$15; lead_snp_bp[$1]=$16; next}
				filenum == 2 && $3 in id && !seen[$3]++ { study[$3]=$1; trait[$3]=$2; p[$3]=$4; next }
				filenum == 3 && FNR == 1 { print; next }
				filenum == 3 { print }
				filenum == 3 && $1 != locus { locus=$1; for (key in id) if (locus_count[key]==locus && lead_snp_kb[key] <= 10000) { print locus_count[key], 9999, study[key], trait[key], id[key], p[key], lead_snp[key], lead_snp_kb[key], lead_snp_rsq[key], lead_snp_alleles[key], chr[key], bp[key], lead_snp_chr[key], lead_snp_bp[key]; delete id[key]}}
				END { for (key in id) { print 9999, 9999, study[key], trait[key], id[key], p[key], lead_snp[key], lead_snp_kb[key], lead_snp_rsq[key], lead_snp_alleles[key], chr[key], bp[key], lead_snp_chr[key], lead_snp_bp[key]} }' OFS='\t' \
				${targetDir}/snplevel.novelty.prevDiscoveries.multiallelic.clumped.txt ${targetDir}/snplevel.novelty.prevDiscoveries.txt ${targetDir}/snplevel.novelty.clumping/clumped > ${targetDir}/snplevel.novelty.clumping/clumped.tmp
	\mv ${targetDir}/snplevel.novelty.clumping/clumped.tmp ${targetDir}/snplevel.novelty.clumping/clumped
fi

# add literature matches to snplevel results
scriptDir=$(dirname "$0")
Rscript ${scriptDir}/snplevel.novelty.R "${ownDiscoveries}" "${targetDir}/snplevel.novelty.clumping/clumped" "${targetDir}/snplevel.novelty.txt"

# how many prevDiscoveries have been found by id in current genetics dataset?
awk -F'\t' 'FNR == 1 { filenum++ }
			filenum == 1 && !seen1[$3]++ { snp[$3]; count1++; next }
			filenum == 2 && $5 in snp && !seen2[$5]++ { count2++; next }
 			END { print "Previously discovered variants:\t"count1"\nPreviously discovered variants sucessfully matched by id and used for clumping:\t"count2}' OFS='\t' ${targetDir}/snplevel.novelty.prevDiscoveries.txt ${targetDir}/snplevel.novelty.clumping/clumped \
 			> ${targetDir}/snplevel.novelty.log

# how many multiallelic variants 
if [ -e ${targetDir}/snplevel.novelty.prevDiscoveries.multiallelic.txt ]; then
awk -F'\t' '{ count++ } END { print "Multiallelic variants clumped by position:\t" count }' OFS='\t' ${targetDir}/snplevel.novelty.prevDiscoveries.multiallelic.txt >> ${targetDir}/snplevel.novelty.log
else
echo "Multiallelic variants clumped by position only (< 10 Mbp):"$'\t'0 >> ${targetDir}/snplevel.novelty.log
fi

# count novel loci - 13
awk -F'\t' 'NR==1 { next } NR==FNR && $3 != "ownDiscoveries" { locus[$1]; next}
	NR==FNR { next } 
	$3 == "ownDiscoveries" && !($1 in locus) && !($1 in seen) { countNovel++; countTotal++; seen[$1]; next}
	$3 == "ownDiscoveries" && ($1 in locus) && !($1 in seen) { countTotal++; seen[$1]; next}
	END { print "Total count of own independent discoveries:\t"countTotal"\nCount of novel independent discoveries:\t"countNovel }' ${targetDir}/snplevel.novelty.clumping/clumped ${targetDir}/snplevel.novelty.clumping/clumped \
	>> ${targetDir}/snplevel.novelty.log

# clean up
mkdir -p ${targetDir}/snplevel.novelty
rm -f ${targetDir}/snplevel.novelty.clumping/*.bed
rm -f ${targetDir}/snplevel.novelty.clumping/*.bim
rm -f ${targetDir}/snplevel.novelty.clumping/*.fam
mv ${targetDir}/snplevel.novelty.clumping ${targetDir}/snplevel.novelty/
mv ${targetDir}/snplevel.novelty.sumstats* ${targetDir}/snplevel.novelty/
mv ${targetDir}/snplevel.novelty.subs.txt ${targetDir}/snplevel.novelty/
mv ${targetDir}/snplevel.novelty.ownDiscoveries.txt ${targetDir}/snplevel.novelty/
mv ${targetDir}/snplevel.novelty.prevDiscoveries* ${targetDir}/snplevel.novelty/
mv ${targetDir}/snplevel.novelty.variants* ${targetDir}/snplevel.novelty/
\mv $targetDir/snplevel.novelty/snplevel.novelty.clumping/clumped ${targetDir}/snplevel.novelty.clumped.txt
tar -cvzf ${targetDir}/snplevel.novelty.tar.gz --directory=${targetDir} snplevel.novelty --remove-files
chmod -R 770 ${targetDir}/snplevel.novelty*

# echo summary
echo $'\n'"--- novelty analysis summary ---"
cat ${targetDir}/snplevel.novelty.log
echo "--- Novelty analysis finished. --- "

