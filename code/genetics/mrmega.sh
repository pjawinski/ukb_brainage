#!/bin/bash

# ============================================
# === run cross-ancestry MR-MEGA and GWAMA ===
# ============================================

# get arguments
inputFiles="$1" # inputFiles="results/gap_gm/replicate/AFR/sumstats.txt.gz,results/gap_gm/replicate/AMR/sumstats.txt.gz,results/gap_gm/replicate/CSA/sumstats.txt.gz,results/gap_gm/replicate/EAS/sumstats.txt.gz,results/gap_gm/replicate/EUR/sumstats.txt.gz,results/gap_gm/replicate/LIFE/sumstats.txt.gz,results/gap_gm/replicate/MID/sumstats.txt.gz"
targetDir="$2" # targetDir="results/gap_gm/replicate/mrmega.all/"
cols="$3" # cols="CHR,BP,ID,A1,A2,A1_FREQ,BETA,SE,N"

# transform variables
IFS=',' read -r -a inputFilesA <<< "${inputFiles}"
IFS=',' read -r -a cols <<< "${cols}"

# echo settings
echo $'\n'"--- MR-MEGA and GWAMA: settings ---"
echo "inputFiles: ${inputFiles[@]}"
echo "targetDir: ${targetDir}"
echo "CHR: ${cols[0]}"
echo "BP: ${cols[1]}"
echo "ID: ${cols[2]}"
echo "A1: ${cols[3]}"
echo "A2: ${cols[4]}"
echo "A1_FREQ: ${cols[5]}"
echo "BETA: ${cols[6]}"
echo "SE: ${cols[7]}"
echo "N: ${cols[8]}"$'\n'
sleep 10

# set targetdir and make folder
mkdir -p "${targetDir}"
targetDir="$(readlink -f "${targetDir}")"

# run MR-MEGA
echo "Running MR-MEGA." 
> ${targetDir}/mrmega.in
for i in "${!inputFilesA[@]}"; do
	echo "${inputFilesA[$i]}" >> ${targetDir}/mrmega.in
done
(MR-MEGA \
--filelist ${targetDir}/mrmega.in \
--qt \
--name_chr "${cols[0]}" \
--name_pos "${cols[1]}" \
--name_marker "${cols[2]}" \
--name_ea "${cols[3]}" \
--name_nea "${cols[4]}" \
--name_eaf "${cols[5]}" \
--name_beta "${cols[6]}" \
--name_se "${cols[7]}" \
--name_n "${cols[8]}" \
--out ${targetDir}/mrmega &

# run random-effects GWAMA 
echo "Running random-effects GWAMA." 
GWAMA \
--filelist ${targetDir}/mrmega.in \
--random \
--quantitative \
--indel_alleles \
--name_marker "${cols[2]}" \
--name_ea "${cols[3]}" \
--name_nea "${cols[4]}" \
--name_eaf "${cols[5]}" \
--name_beta "${cols[6]}" \
--name_se "${cols[7]}" \
--name_n "${cols[8]}" \
--output ${targetDir}/gwama.re &

# run fixed-effects GWAMA 
echo "Running fixed-effects GWAMA" 
GWAMA \
--filelist ${targetDir}/mrmega.in \
--quantitative \
--indel_alleles \
--name_marker "${cols[2]}" \
--name_ea "${cols[3]}" \
--name_nea "${cols[4]}" \
--name_eaf "${cols[5]}" \
--name_beta "${cols[6]}" \
--name_se "${cols[7]}" \
--name_n "${cols[8]}" \
--output ${targetDir}/gwama.fe
wait)

# fix MR-MEGA p
scriptDir=$(dirname "$0")
Rscript "${scriptDir}"/mrmega.fixP.R --args input="${targetDir}"/mrmega.result out="${targetDir}"/mrmega.result

# run qc on MR-MEGA results and remove variants with n < quantile(Nmax, 0.9) / 1.5
echo "Removing variants with low n and low cohort count."
ncriterion=$(awk '$2>0 && $2<24 { print $7 }' "${targetDir}/mrmega.result" | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.90 - 0.5)]/1.5}')
ncriterionY=$(awk '$2==24 { print $7 }' "${targetDir}/mrmega.result" | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.90 - 0.5)]/1.5}')
ncriterionMT=$(awk '$2==0 { print $7 }' "${targetDir}/mrmega.result" | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.90 - 0.5)]/1.5}')
mawk -v ncriterion="${ncriterion}" -v ncriterionY="${ncriterionY}" -v ncriterionMT="${ncriterionMT}" '
	$2==23 { $2="X" } $2==24 { $2="Y" } $2==25 { $2="XY" } $2==0 { $2="MT" } 
	NR==1 || $10!="NA" && (($2!="Y" && $2!="MT" && $7 > ncriterion) || ($2=="Y" && $7 > ncriterionY) || ($2=="MT" && $7 > ncriterionMT)) { print }' "${targetDir}/mrmega.result" OFS='\t' > "${targetDir}/mrmega.qc"

# join results
echo "Adding effect sizes from random- and fixed-effects GWAMA to MR-MEGA output."
mawk '{ gsub("\t"," ",$0) }
	FNR==1 { filenum++; if(filenum==3) { $1="ID";$2="CHR";$3="BP";$4="A1";$5="A2";$6="A1_FREQ";$7="N";gsub("chisq_association","chisq");gsub("ndf_association","ndf");gsub("P-value_association","P"); print $0" gwama_re_beta gwama_re_se gwama_re_pval gwama_fe_beta gwama_fe_se gwama_fe_pval" }; next }
	filenum==1 { re[$1":"$2":"$3]=$5" "$6" "$10; reinverted[$1":"$3":"$2]=-1*$5" "$6" "$10; next }
	filenum==2 { fe[$1":"$2":"$3]=$5" "$6" "$10; feinverted[$1":"$3":"$2]=-1*$5" "$6" "$10; next }
	filenum==3 { output=$0; if($1":"$4":"$5 in re) { output=output" "re[$1":"$4":"$5] }
	 				   else if($1":"$5":"$4 in re) { output=output" "re[$1":"$5":"$4] }
	 					   else { output=output" NA NA NA" }; 
	 					     if($1":"$4":"$5 in fe) { output=output" "fe[$1":"$4":"$5] }
	 					else if($1":"$5":"$4 in fe) { output=output" "fe[$1":"$5":"$4] }
	 					   else { output=output" NA NA NA" }; 
	 			  print output }' \
	 "${targetDir}/gwama.re.out" "${targetDir}/gwama.fe.out" "${targetDir}/mrmega.qc" > "${targetDir}/mrmega.weights"

# clean up folder
echo "Cleaning up folder."
rm -f ${targetDir}/mrmega.qc
pigz -f ${targetDir}/mrmega.result
pigz -f ${targetDir}/mrmega.weights
pigz -f ${targetDir}/gwama.re.out
pigz -f ${targetDir}/gwama.re.err.out
pigz -f ${targetDir}/gwama.re.log.out
pigz -f ${targetDir}/gwama.fe.out
pigz -f ${targetDir}/gwama.fe.err.out
pigz -f ${targetDir}/gwama.fe.log.out
chmod 770 ${targetDir}/mrmega*
chmod 770 ${targetDir}/gwama.re*
chmod 770 ${targetDir}/gwama.fe*
echo "--- Completed: MR-MEGA and GWAMA --- "
