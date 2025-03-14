#!/bin/bash

# ====================================================================
# === prepare 'novelty' (mapping with previous literature results) ===
# ====================================================================

# set working directory
cd /slow/projects/ukb_brainage
conda activate envs/default

# set targetdir
targetDir="data/prevDiscoveries"

# download association files of Smith et al. (2020) and  Ning et al. (2020)
mkdir -p ${targetDir}
wget -O ${targetDir}/Smith_2020.txt https://www.fmrib.ox.ac.uk/ukbiobank/BrainAgingModes/GWAS/GWAS_2019_09_02_1/Table1.txt
wget -O ${targetDir}/Ning_2020.xlsx https://zenodo.org/record/3496206/files/Supplementary_Table_3_SNP_p_value.xlsx?download=1
wget -O ${targetDir}/Ning_2021.xlsx https://zenodo.org/records/3786826/files/Supplementary_Table_1_SNP_RBA_assoc_P_values.xlsx?download=1

# add association results of Jonsson et al. (2019) shown in main article (Table 4)
echo "TRAIT"$'\t'"ID"$'\t'"P"$'\n'\
"jgwm"$'\t'"rs2435204"$'\t'"1.4e-12"$'\n'\
"jgwm"$'\t'"rs1452628"$'\t'"2.3e-09"$'\n'\
"jgwm"$'\t'"rs2790099"$'\t'"8.9e-06"$'\n'\
"jgwm"$'\t'"rs6437412"$'\t'"6.8e-06"$'\n'\
"jgwm"$'\t'"rs2184968"$'\t'"7.5e-05" > ${targetDir}/Jonsson_2019.txt

# add association results of Leonardsen et al. (2023) shown in main article
echo "TRAIT"$'\t'"ID"$'\t'"P"$'\n'\
"t1w"$'\t'"rs73185796"$'\t'"2.53e-08"$'\n'\
"t1w"$'\t'"rs13132853"$'\t'"2.34e-18"$'\n'\
"t1w"$'\t'"rs79107704"$'\t'"1.65e-08"$'\n'\
"t1w"$'\t'"rs2790102"$'\t'"8.92e-09"$'\n'\
"t1w"$'\t'"rs7461069"$'\t'"1.57e-08"$'\n'\
"t1w"$'\t'"rs4880424"$'\t'"3.69e-08"$'\n'\
"t1w"$'\t'"rs17203398"$'\t'"1.42e-10"$'\n'\
"t1w"$'\t'"rs2106786"$'\t'"1.87e-23" > ${targetDir}/Leonardsen_2023.txt

# add association results of Kim et al. (2023) shown in main article (Table 2)
# replace -(deletion) on chr10 with 10:134559486_CGCCTACCCT_C (cat data/genetics/chr10/imp_mri/chr10_mri.pvar | grep 134559486)
# replace -(deletion) on chr7 with 7:120981641_CT_C (cat data/genetics/chr7/imp_mri/chr7_mri.pvar | grep 120981641)
echo "TRAIT"$'\t'"ID"$'\t'"P"$'\n'\
"t1wcnn"$'\t'"rs147431626"$'\t'"5.36e-22"$'\n'\
"t1wcnn"$'\t'"rs201791735"$'\t'"1.14e-15"$'\n'\
"t1wcnn"$'\t'"10:134559486_CGCCTACCCT_C"$'\t'"1.10e-11"$'\n'\
"t1wcnn"$'\t'"rs534102361"$'\t'"3.46e-10"$'\n'\
"t1wcnn"$'\t'"rs536165403"$'\t'"1.64e-08"$'\n'\
"t1wcnn"$'\t'"7:120981641_CT_C"$'\t'"1.90e-08" > ${targetDir}/Kim_2023.txt

# add association results of Wen et al. (2024; https://doi.org/10.1101%2F2023.04.13.536818) shown in supplementary eTable 2
# note that rs534115641 was replaced by rs62064415 (rs534115641 was merged into rs62064415 on March 31, 2015 (Build 144))
# note that rs34051980 was replaced by 8:119948046_AGAGG_A (rs34051980 was merged into rs3081402 on October 12, 2018 (Build 152) -> in UK Biobank: 8:119948046_AGAGG_A)
echo "TRAIT"$'\t'"ID"$'\t'"P"$'\n'\
"wgm"$'\t'"rs61732315"$'\t'"1.63e-08"$'\n'\
"wgm"$'\t'"rs1452628"$'\t'"3.04e-14"$'\n'\
"wgm"$'\t'"rs186399184"$'\t'"1.05e-08"$'\n'\
"wgm"$'\t'"rs10933668"$'\t'"9.41e-09"$'\n'\
"wgm"$'\t'"8:119948046_AGAGG_A"$'\t'"2.02e-08"$'\n'\
"wgm"$'\t'"rs62064415"$'\t'"8.44e-23"$'\n'\
"wwm"$'\t'"rs11118475"$'\t'"3.69e-09"$'\n'\
"wwm"$'\t'"rs61067594"$'\t'"6.01e-19"$'\n'\
"wwm"$'\t'"rs967140"$'\t'"7.26e-12"$'\n'\
"wwm"$'\t'"rs2533872"$'\t'"9.95e-10"$'\n'\
"wwm"$'\t'"rs564819152"$'\t'"9.39e-13"$'\n'\
"wwm"$'\t'"rs12146713"$'\t'"7.67e-14"$'\n'\
"wwm"$'\t'"rs654276"$'\t'"1.96e-08"$'\n'\
"wwm"$'\t'"rs4843550"$'\t'"2.84e-09"$'\n'\
"wwm"$'\t'"rs1894525"$'\t'"2.71e-09"$'\n'\
"wcf"$'\t'"rs5877290"$'\t'"2.31e-08" > ${targetDir}/Wen_2024.txt

# munge sumstats of Smith et al., Jonsson et al., and Leonardsen et al.
awk 'NR==1 { next } $2 != "V0141" { print "Smith_2020", $2, $1, 10^(-$9) }' OFS='\t' ${targetDir}/Smith_2020.txt > ${targetDir}/Smith_2020_munged.txt
awk 'NR==1 { next } { print "Jonsson_2019", $0 }' OFS='\t' ${targetDir}/Jonsson_2019.txt > ${targetDir}/Jonsson_2019_munged.txt
awk 'NR==1 { next } { print "Kim_2023", $0 }' OFS='\t' ${targetDir}/Kim_2023.txt > ${targetDir}/Kim_2023_munged.txt
awk 'NR==1 { next } { print "Leonardsen_2023", $0 }' OFS='\t' ${targetDir}/Leonardsen_2023.txt > ${targetDir}/Leonardsen_2023_munged.txt
awk 'NR==1 { next } { print "Wen_2024", $0 }' OFS='\t' ${targetDir}/Wen_2024.txt > ${targetDir}/Wen_2024_munged.txt

# munge sumstats of Ning et al. 2020
xlsx2csv -s 1 -d "\t" --floatformat %e ${targetDir}/Ning_2020.xlsx | awk 'NR==1 { print } $4 < 5E-8 { print }' OSF='\t' > ${targetDir}/Ning_2020.txt

	# get chr and pos
	awk 'NR==1 { print; next } NR==FNR { snp[$2"_"$3]=$2"\t"$3"\t"$4; next } $1"_"$2 in snp { print $3, snp[$1"_"$2] }' OFS='\t' ${targetDir}/Ning_2020.txt data/genetics/chr17/imp_mri/chr17_mri.pvar > ${targetDir}/Ning_2020.tmp.txt

	# any snp missing? - one SNP in same locus (rs145809511). drop it.
	awk 'NR==FNR { snp[$2"_"$3]; next } !($2"_"$3 in snp) { print }' ${targetDir}/Ning_2020.tmp.txt ${targetDir}/Ning_2020.txt 
	\mv ${targetDir}/Ning_2020.tmp.txt ${targetDir}/Ning_2020.txt

	# create munged file
	awk 'NR>1 { print "Ning_2020", "Desikan-Killiany", $1, $4 }' OFS='\t' ${targetDir}/Ning_2020.txt > ${targetDir}/Ning_2020_munged.txt

# munge sumstats of Ning et al. 2021
xlsx2csv -s 1 -d "\t" --floatformat %e ${targetDir}/Ning_2021.xlsx | awk 'NR==1 { print } $10 < 5E-8 { print }' OSF='\t' > ${targetDir}/Ning_2021.txt

	# get correct rsids (instead of Affx-xxx) based on chr_bp_a1_a2 or chr_bp_a2_a1
	pvar=$(for i in $(cat ${targetDir}/Ning_2021.txt | awk 'NR>1 { print $2 }' | sort -u); do echo data/genetics/chr${i}/imp_mri/chr${i}_mri.pvar; done)
	npvar=$(echo "$pvar" | awk 'END {print NR }'); 
	awk -v npvar=${npvar} 'FNR==1 {filenum++}
		filenum<=npvar { chr_bp_a1_a2[$1"_"$2"_"$4"_"$5]=$3; next }
		filenum>npvar && FNR==1 { print "ID", $0; next }
		$2"_"$3"_"$4"_"$5 in chr_bp_a1_a2 { print chr_bp_a1_a2[$2"_"$3"_"$4"_"$5], $0; next }
		$2"_"$3"_"$5"_"$4 in chr_bp_a1_a2 { print chr_bp_a1_a2[$2"_"$3"_"$5"_"$4], $0 }' OFS='\t' ${pvar} ${targetDir}/Ning_2021.txt  > ${targetDir}/Ning_2021.tmp.txt

	# any snp missing? - same as Ning et al. 2020: one SNP in 17q21.31 locus (rs145809511). drop it.
	awk 'NR==FNR { snp[$2]; next } !($1 in snp) { print }' ${targetDir}/Ning_2021.tmp.txt ${targetDir}/Ning_2021.txt 
	\mv ${targetDir}/Ning_2021.tmp.txt ${targetDir}/Ning_2021.txt

	# create munged file
	awk 'NR>1 { print "Ning_2021", "Desikan-Killiany", $1, $11 }' OFS='\t' ${targetDir}/Ning_2021.txt > ${targetDir}/Ning_2021_munged.txt

	# compare variants in Ning et al. 2020 and 2021 - all variants of Ning et al. 2020 also reported in Ning et al. 2021
	awk 'NR==FNR { snp[$3]; next } !($3 in snp) { print }' ${targetDir}/Ning_2021_munged.txt ${targetDir}/Ning_2020_munged.txt

# change access rights
chmod -R 770 ${targetDir}
