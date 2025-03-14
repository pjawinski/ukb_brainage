#!/bin/bash

# ============================================
# === Get SBayes credible sets of variants ===
# ============================================

# get arguments
out="${1}" # out="results/gap_gm/gwama/eur/sbayes/sbayesrc"
sbayesRes="${2}" # sbayesRes="results/gap_gwm/gwama/eur/sbayes/sbayesrc"
ldFile="${3}" # ldFile="/fast/software/gctb/resources/ukbEUR_32k/ldm/pairwise0.5.ld.txt"
pip=${4} # pip=0.95
pep=${5} # pep=0.50
threads=${6} # threads=100
conditionalFile="${7}" # conditionalFile="results/gap_gwm/gwama/eur/conditional/conditional.cleaned.tophits.txt"
pthresh=${8} # pthresh=5e-8

# echo settings
echo $'\n'"--- Get SBayes credible sets of variants | Settings ---"
echo "out: ${out}"
echo "sbayesRes: ${sbayesRes}"
echo "ldFile: ${ldFile}"
echo "pip: ${pip}"
echo "pip: ${pep}"
echo "threads: ${threads}"
echo "conditionalFile: ${conditionalFile}"
echo "pthresh: ${pthresh}"$'\n'
sleep 10

# set targetdir and make folder
targetDir=$(basename "${out}")
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# get credible sets
pigz -dc "${sbayesRes}.snpRes.gz" > "${sbayesRes}.snpRes"
gctb \
--cs \
--ld-file "${ldFile}" \
--pip "${pip}" \
--pep "${pep}" \
--mcmc-samples "${sbayesRes}" \
--thread "${threads}" \
--out "${out}"

# assign leadsnps based on ld data (tophits with p < 1e-6 and R2 > 0.8)
awk -v cols="ID,LEAD_SNP,LEAD_SNP_P" -v pthresh="${pthresh}" '
    BEGIN { ncols=split(cols,colnames,",") }
    NR==1 { for(i=1;i<=ncols;++i) { for(j=1;j<=NF;++j) { if($j==colnames[i]) { colidx[i]=j } } } }
    NR==1 { id[$colidx[1]]=$colidx[2]; for(i=3;i<=ncols;++i) { id[$colidx[1]]=id[$colidx[1]]"\t"$colidx[i] }; next}
    NR==FNR && $colidx[3] < pthresh { id[$colidx[1]]=$colidx[2]; for(i=3;i<=ncols;++i) { id[$colidx[1]]=id[$colidx[1]]"\t"$colidx[i] }; next}
    NR==FNR { next }
    FNR==1 { print id[colnames[1]], $2, $3, $4, $5, $6, $7; next } 
    { nsnps=split($7,snps,",") }
    { for(i=1;i<=nsnps;++i) { if(snps[i] in id) { print id[snps[i]], $2, $3, $4, $5, $6, $7; break } }; next }
    ' OFS='\t' "${conditionalFile}" "${sbayesRes}.lcs" \
    > "${out}.lcs.leadsnps"

# if PEP < 0.7, only keep lcs with largest PEP 
PEPidx=$(awk -v col="PEP" 'NR==1 { for(i=1;i<=NF;++i) { if($i==col) { print i } } }' "${out}.lcs.leadsnps")
(head -n 1 "${out}.lcs.leadsnps" && tail -n +2 "${out}.lcs.leadsnps" | sort -k1,1 -k${PEPidx},${PEPidx}r) > "${out}.lcs.leadsnps.tmp"
awk -v PEPidx="${PEPidx}" 'NR==1 { print; next } (!seen[$1]++ || $PEPidx >= 0.7) { print }
    ' "${out}.lcs.leadsnps.tmp" > "${out}.lcs.leadsnps"
rm -f "${out}.lcs.leadsnps.tmp"

# report results of multiple signals per locus in one line
awk -F'\t' 'NR==1 { print $0, "numCausal" }
    NR==FNR && !seen[$1] { id[$1]=$1"\t"$2; size[$1]=$3; pip[$1]=sprintf("%.2f",$4); pgv[$1]=$5; pgvenrich[$1]=$6; pep[$1]=sprintf("%.2f",$7); snp[$1]=$8; seen[$1]++; next}
    NR==FNR && seen[$1]++ { size[$1]=size[$1]" | "$3; pip[$1]=pip[$1]" | "sprintf("%.2f",$4); pgv[$1]=pgv[$1]" | "$5; pgvenrich[$1]=pgvenrich[$1]" | "$6; pep[$1]=pep[$1]" | "sprintf("%.2f",$7); snp[$1]=snp[$1]" | "$8; next}
    FNR==1 { next }
    !seen2[$1]++ { print id[$1], size[$1], pip[$1], pgv[$1], pgvenrich[$1], pep[$1], snp[$1], seen[$1] }
    ' OFS='\t' "${out}.lcs.leadsnps" "${out}.lcs.leadsnps" > "${out}.lcs.leadsnps.joined"

# get details of variants in credible sets
SNPidx=$(awk -v col="SNP" 'NR==1 { for(i=1;i<=NF;++i) { if($i==col) { print i } } }' "${out}.lcs.leadsnps")
awk -v col="${SNPidx}" '
    NR==FNR { id[$2]=$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16; next }
    FNR==1 { print "LEAD_SNP","csCount","csSize","csPIP","csPEP","csPGV","csPGVenrich",id["Name"]; next }
    { count[$1]++; nsnps=split($col,snps,","); for(i=1;i<=nsnps;++i) { print $1,count[$1],$3,$4,$7,$5,$6,id[snps[i]] } }
    ' OFS='\t' "${sbayesRes}.snpRes" "${out}.lcs.leadsnps" > "${out}.lcs.leadsnps.snpRes"

# add cumulative PIP
awk 'NR==1 { print $0, "cumPIP"; next}
    !seen[$1":"$2]++ { cumPIP=0 } { print $0, cumPIP+$NF; cumPIP=cumPIP+$NF }
    ' OFS='\t' "${out}.lcs.leadsnps.snpRes" > "${out}.lcs.leadsnps.snpRes.tmp"
\mv "${out}.lcs.leadsnps.snpRes.tmp" "${out}.lcs.leadsnps.snpRes"

# clean up
echo "Cleaning up."
rm -f "${sbayesRes}.snpRes"
chmod 770 "${out}."*
echo "--- Completed: Get SBayes credible sets  --- "
