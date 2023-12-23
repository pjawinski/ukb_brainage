#!/bin/bash

# ====================================================================
# === prepare 'novelty' (mapping with previous literature results) ===
# ====================================================================

# set working directory
cd /slow/projects/ukb_brainage

# set targetdir
targetDir="data/prevDiscoveries"

# download association files of Smith et al. (2020) and  Ning et al. (2020)
mkdir -p ${targetDir}
wget -O ${targetDir}/Smith_2020.txt https://www.fmrib.ox.ac.uk/ukbiobank/BrainAgingModes/GWAS/GWAS_2019_09_02_1/Table1.txt
wget -O ${targetDir}/Ning_2020.xlsx https://zenodo.org/record/3496206/files/Supplementary_Table_3_SNP_p_value.xlsx?download=1

# add association results of Jonsson et al. (2019) shown in main article (Table 4)
echo "ID"$'\t'"P"$'\n'\
"rs2435204"$'\t'"1.4e-12"$'\n'\
"rs1452628"$'\t'"2.3e-09"$'\n'\
"rs2790099"$'\t'"8.9e-06"$'\n'\
"rs6437412"$'\t'"6.8e-06"$'\n'\
"rs2184968"$'\t'"7.5e-05" > ${targetDir}/Jonsson_2019.txt

# add association results of Leonardsen et al. (2023) shown in main article
echo "ID"$'\t'"P"$'\n'\
"rs73185796"$'\t'"2.53e-08"$'\n'\
"rs13132853"$'\t'"2.34e-18"$'\n'\
"rs79107704"$'\t'"1.65e-08"$'\n'\
"rs2790102"$'\t'"8.92e-09"$'\n'\
"rs7461069"$'\t'"1.57e-08"$'\n'\
"rs4880424"$'\t'"3.69e-08"$'\n'\
"rs17203398"$'\t'"1.42e-10"$'\n'\
"rs2106786"$'\t'"1.87e-23" > ${targetDir}/Leonardsen_2023.txt

# munge sumstats of Smith et al., Jonsson et al., and Leonardsen et al.
awk 'NR==1 { next } { print "Smith_2020", $2, $1, 10^(-$9) }' OFS='\t' ${targetDir}/Smith_2020.txt > ${targetDir}/Smith_2020_munged.txt
awk 'NR==1 { next } { print "Jonsson_2019", "jgwm", $1, $2 }' OFS='\t' ${targetDir}/Jonsson_2019.txt > ${targetDir}/Jonsson_2019_munged.txt
awk 'NR==1 { next } { print "Leonardsen_2023", "t1w", $1, $2 }' OFS='\t' ${targetDir}/Leonardsen_2023.txt > ${targetDir}/Leonardsen_2023_munged.txt

# munge sumstats of Ning et al.
xlsx2csv -s 1 -d "\t" --floatformat %e ${targetDir}/Ning_2020.xlsx | awk 'NR==1 { print } $4 < 5E-8 { print }' OSF='\t' > ${targetDir}/Ning_2020.txt

	# get chr and pos
	awk 'NR==1 { print; next } NR==FNR { snp[$2"_"$3]=$2"\t"$3"\t"$4; next } $1"_"$2 in snp { print $3, snp[$1"_"$2] }' OFS='\t' ${targetDir}/Ning_2020.txt data/genetics/chr17/imp_mri/chr17_mri.pvar > ${targetDir}/Ning_2020.tmp.txt

	# any snp missing? - one SNP in same locus (rs145809511). drop it.
	awk 'NR==FNR { snp[$2"_"$3]; next } !($2"_"$3 in snp) { print }' ${targetDir}/Ning_2020.tmp.txt ${targetDir}/Ning_2020.txt 
	\mv ${targetDir}/Ning_2020.tmp.txt ${targetDir}/Ning_2020.txt

	# create munged file
	awk 'NR>1 { print "Ning_2020", "Desikan-Killiany", $1, $4 }' OFS='\t' ${targetDir}/Ning_2020.txt > ${targetDir}/Ning_2020_munged.txt

# change access rights
chmod -R 770 ${targetDir}
