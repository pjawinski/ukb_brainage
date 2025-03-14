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