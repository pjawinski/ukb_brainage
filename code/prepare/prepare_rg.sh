#!/bin/bash

# ====================================
# === prepare genetic correlations ===
# ====================================
# see prepare_ldsc.sh for general ldsc preparations

# set working directory
cd /slow/projects/ukb_brainage/data/ldsc/resources/

# get dropbox links (see https://nealelab.github.io/UKBB_ldsc/downloads.html)
wget -O neale_manifest_dropbox.tsv "https://docs.google.com/spreadsheets/d/1EmwlCYYqkoVKqAS71nNDKoN18PyKwysYUppKTSMvKiM/export?exportFormat=tsv"

# remove ^M at the end of each row
cat -v neale_manifest_dropbox.tsv | sed -e "s/\^M//g" | sed -e "s/GeneId://g" > neale_manifest_dropbox.tsv.tmp; \mv neale_manifest_dropbox.tsv.tmp neale_manifest_dropbox.tsv; 

# extract download links of heritable traits
awk -F'\t' '$18=="z4" || $18 == "z7" || $18 == "nominal" { print }' neale_manifest_dropbox.tsv | head
awk -F '\t' '$18=="z4" || $18 == "z7" || $18 == "nominal" { print $7 }' neale_manifest_dropbox.tsv > neale_manifest_getFiles.sh

# download files (nparallel at a time)
mkdir -p neale; cd neale
nparallel=50
count=0
while read line; do
    ${line} &
    (( count ++ ))        
    if (( count == nparallel )); then
        wait
        count=0
    fi
done < ../neale_manifest_getFiles.sh

# get category hierarchy from ukb showcase
cd ../
showcaseLinks=$(awk -F '\t' '$3 != "N/A" && ($18=="z4" || $18 == "z7" || $18 == "nominal") { print $3; next}' neale_manifest_dropbox.tsv | sort -u)
rm -rf neale_showcase_categories.txt
for i in $showcaseLinks; do
	showcaseContent="`wget -qO- $i`"
	category=$(echo "${showcaseContent}" | grep Category: | sed -e 's/<[^>]*>//g' | sed -e 's/Category://g')
	echo "$i"$'\t'"$category" >> neale_showcase_categories.txt
done


# ===================================================
# === download sumstats of major mental disorders ===
# ===================================================

# make directory
mkdir -p data/ldsc/resources/sumstatsCollection
cd data/ldsc/resources/sumstatsCollection

# get base files
rsync -auP jawinskp@cluster1.psychologie.hu-berlin.de:/home/groups/markett/UK_Biobank/05_GWAS_for_PRS/00_base_gz/ .

# set working directory and activate ldsc conda environment
cd /slow/projects/ukb_brainage
conda activate envs/ldsc

# settings
LDsnplist=data/ldsc/resources/w_hm3.noMHC.snplist
targetDir=data/ldsc/resources/sumstatsCollection/munged/
files=$(ls data/ldsc/resources/sumstatsCollection/*txt.gz)

# munge sumstats
echo "Munging sumstats."
mkdir -p ${targetDir}
for file in $files; do
    outFile=$(echo $file | sed 's/.*\///g' | sed 's/.txt.gz//g')
    munge_sumstats.py \
    --sumstats "${file}" \
    --merge-alleles "${LDsnplist}" \
    --out "${targetDir}/${outFile}"
    chmod -R 770 ${targetDir}/${outFile}*
done

# ===============================================================
# === download Freesurfer GWAS sumstats from UK Biobank BIG40 ===
# ===============================================================

# list of phenotypes: https://open.win.ox.ac.uk/ukbiobank/big40/BIG40-IDPs_v4/IDPs.html
# download: curl -O -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/0001.txt.gz
# 0189-0221 # aseg_[lh,rh]_volume (33)
# 0344-0409 # aparc-Desikan_[lh,rh]_volume (66)
# 0648-0715 # aparc-Desikan_[lh,rh]_area (68)
# 1020-1087 # aparc-Desikan_[lh,rh]_thickness (68)

LDsnplist=/fast/software/ldsc/resources/w_hm3.noMHC.snplist
mkdir -p /fast/software/ldsc/resources/brainGWAS/munged/
cd /fast/software/ldsc/resources/brainGWAS
targetDir=munged

N=50
(for i in {189..221} {344..409} {648..715}; do
    ((j=j%N)); ((j++==0)) && wait
    (curl -O -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/0${i}.txt.gz
    zcat 0${i}.txt.gz | awk 'NR==1 { sub(/pval\(-log10\)/,"pval"); print } 
        NR > 1 { sub("^0*","",$1); print $1, $2, $3, $4, $5, $6, $7, 10^(-$8) }' > "0${i}.prep"
    N=$(cat brainGWAS.txt | awk -F'\t' -v pheno="0${i}" '$1==pheno { print $16 }')
    outFile=$(cat brainGWAS.txt | awk -F'\t' -v pheno="0${i}" '$1==pheno { print $3 }')
    munge_sumstats.py \
    --sumstats 0${i}.prep \
    --merge-alleles "${LDsnplist}" \
    --N ${N} \
    --out "${targetDir}/${outFile}"
    chmod -R 770 ${targetDir}/${outFile}*
    rm -f "0${i}.prep") &
done
wait)

N=50
(for i in {1020..1087}; do
    ((j=j%N)); ((j++==0)) && wait
    (curl -O -L -C - https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/${i}.txt.gz
    zcat ${i}.txt.gz | awk 'NR==1 { sub(/pval\(-log10\)/,"pval"); print } 
        NR > 1 { sub("^0*","",$1); print $1, $2, $3, $4, $5, $6, $7, 10^(-$8) }' > "${i}.prep"
    N=$(cat brainGWAS.txt | awk -F'\t' -v pheno="0${i}" '$1==pheno { print $16 }')
    outFile=$(cat brainGWAS.txt | awk -F'\t' -v pheno="0${i}" '$1==pheno { print $3 }')
    munge_sumstats.py \
    --sumstats ${i}.prep \
    --merge-alleles "${LDsnplist}" \
    --N ${N} \
    --out "${targetDir}/${outFile}"
    chmod -R 770 ${targetDir}/${outFile}*
    rm -f "${i}.prep") &
done
wait)

# ===============================================
# === download ENIGMA3 cortical GWAS sumstats ===
# ===============================================

mkdir -p /fast/software/ldsc/resources/ENIGMA3/munged/
cd /fast/software/ldsc/resources/ENIGMA3/

# get list of files
wget https://enigma.ini.usc.edu/downloads/ENIGMA3_Global/
files=$(cat index.html | grep .txt.gz | sed 's%.*ENIGMA3%https://enigma.ini.usc.edu/downloads/ENIGMA3_Global/ENIGMA3%g' | sed 's/.txt.gz.*/.txt.gz/g')
rm -f index.html
wget https://enigma.ini.usc.edu/downloads/ENIGMA3_withoutGlobal/
files=$(echo $files $(cat index.html| grep .txt.gz | sed 's%.*ENIGMA3%https://enigma.ini.usc.edu/downloads/ENIGMA3_withoutGlobal/ENIGMA3%g' | sed 's/.txt.gz.*/.txt.gz/g'))
rm -f index.html

# download files in parallel and munge them (calculate Z, because BETA1 shows large values and check.median exceeds tolerance from 0)
LDsnplist=/fast/software/ldsc/resources/w_hm3.noMHC.snplist
targetDir=munged
for i in ${files}; do (
    curl -O -C - ${i}
    iDownl=$(echo ${i} | sed 's%.*/%%g')
    chmod 750 ${iDownl}
    outFile=$(echo ${iDownl} | sed 's/.txt.gz//g')
    zcat ${iDownl} | awk 'NR==1 { print $0, "Z" } NR>1 { print $0, $5/$6 }' > ${outFile}.prep
    munge_sumstats.py \
        --sumstats ${outFile}.prep \
        --signed-sumstats Z,0 \
        --merge-alleles "${LDsnplist}" \
        --out "${targetDir}/${outFile}"
    chmod -R 750 ${targetDir}/${outFile}*
    rm -f ${outFile}.prep
    ) &
done
