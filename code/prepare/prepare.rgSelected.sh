#!/bin/bash

# ===================================================================================================
# === Selected traits for genetic correlation analyses: Download, harmonize, and munge GWAS files ===
# ===================================================================================================

# set working directory
cd /slow/projects/ukb_brainage
conda activate envs/default
mkdir -p data/sumstats/raw
mkdir -p data/sumstats/harmonized

# ===================
# === Psychiatric ===
# ===================

# PGC_ADHD_2023 (Demontis et al., https://doi.org/10.1038/s41588-022-01285-8), data doi: 10.6084/m9.figshare.22564390
mkdir -p data/sumstats/raw/01_adhd_demontis_2023 && cd data/sumstats/raw/01_adhd_demontis_2023
wget --content-disposition https://figshare.com/ndownloader/files/40036699 # https://figshare.com/articles/dataset/adhd2022/22564390?file=40036699
wget --content-disposition https://figshare.com/ndownloader/files/40036684 # https://figshare.com/articles/dataset/adhd2022/22564390?file=40036684
cd -
sumstats="data/sumstats/raw/01_adhd_demontis_2023/ADHD2022_iPSYCH_deCODE_PGC.meta.gz" # GWAS file with summary statistics
outFile="data/sumstats/harmonized/01_adhd_demontis_2023" # set outputFile (suffix .gz and .log will be attached)
skipLines=0 # needs to be set of header is not first line
sep="auto" # set delimiter manually if necessary
colsIn="CHR,BP,SNP,A1,A2,FRQ_U_186843,INFO,OR,SE,P,Nca,Nco" # input column names
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,INFO,OR,SE,P,Nca,Nco" # output column names (use column names CHR,BP,ID,A1,A2,A1_FREQ,INFO,OR,BETA,SE,P,Nca,Nco,Neff_half, and N; otherwise functions like "addRSID" may not work)
halveNeff=FALSE # divide Neff by 2 if it was calculated by 4/(1/Nca + 1/Nco)
setN=-1 # set N manually if no N column is provided
setNca=-1 # set number of cases
setNco=-1 # set number of controls
compNeff_half=TRUE # compute Neff_half as 2/(1/Nca + 1/Nco)
compN=TRUE # compute total N as Nca+Nco
compMAF=TRUE # compute MAF from A1_FREQ
invertCol=FALSE # invert effect size if necessary (e.g. Z column for Baselmans et al. neuroticism GWAS, which was coded in the direction of 'Well-being')
rmDuplicatesCol=FALSE # define output column for excluding duplicated variant ids (first occurence will be kept)
addRSID=FALSE # add RSID if not provided (requires preprocessed 1kgp files)
kgpFiles=FALSE # 1kgp file handler (requires preprocessed 1kgp files)
keepKGPfreq=FALSE # keep 1kgp columns or drop them (1KGP REF allele, ALT allele, CEU allele frequency)
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# PGC_ASD_2019 (Grove et al., https://doi.org/10.1038/s41588-019-0344-8), data doi: 10.6084/m9.figshare.14671989 
mkdir -p data/sumstats/raw/01_asd_grove_2019 && cd data/sumstats/raw/01_asd_grove_2019
wget --content-disposition https://figshare.com/ndownloader/files/28169289 # https://figshare.com/articles/dataset/asd2019/14671989?file=28169289
wget --content-disposition https://figshare.com/ndownloader/files/28169292 # https://figshare.com/articles/dataset/asd2019/14671989?file=28169292
cd -
sumstats="data/sumstats/raw/01_asd_grove_2019/iPSYCH-PGC_ASD_Nov2017.gz"
outFile="data/sumstats/harmonized/01_asd_grove_2019"
skipLines=0
sep="auto"
colsIn="CHR,BP,SNP,A1,A2,INFO,OR,SE,P"
colsOut="CHR,BP,ID,A1,A2,INFO,OR,SE,P"
halveNeff=FALSE
setN=-1
setNca=8382
setNco=27969
compNeff_half=TRUE
compN=TRUE
compMAF=FALSE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# PGC_BIP_2021 (Mullins et al., https://doi.org/10.1038/s41588-021-00857-4) data doi:  10.6084/m9.figshare.14102594
mkdir -p data/sumstats/raw/01_bip_mullins_2021 && cd data/sumstats/raw/01_bip_mullins_2021
wget --content-disposition https://figshare.com/ndownloader/files/26603681 # https://figshare.com/articles/dataset/PGC3_bipolar_disorder_GWAS_summary_statistics/14102594?file=26603681
wget --content-disposition https://figshare.com/ndownloader/files/40036705 # no ukbb: https://figshare.com/articles/dataset/bip2021_noUKBB/22564402?file=40036705
cd -

sumstats="data/sumstats/raw/01_bip_mullins_2021/pgc-bip2021-all.vcf.tsv.gz"
outFile="data/sumstats/harmonized/01_bip_mullins_2021"
skipLines=72
sep="auto"
colsIn="#CHROM,POS,ID,A1,A2,FCON,IMPINFO,BETA,SE,PVAL,NCAS,NCON"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,INFO,BETA,SE,P,Nca,Nco"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=TRUE
compN=TRUE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

sumstats="data/sumstats/raw/01_bip_mullins_2021/daner_bip_pgc3_nm_noukbiobank.gz"
outFile="data/sumstats/harmonized/01_bip_mullins_2021_noUKB"
skipLines=0
sep="auto"
colsIn="CHR,BP,SNP,A1,A2,FRQ_U_313436,INFO,OR,SE,P,Nca,Nco,Neff_half"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,INFO,OR,SE,P,Nca,Nco,Neff_half"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=TRUE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# PGC_AN_2019 (Watson et al., 10.1038/s41588-019-0439-2) data doi: 10.6084/m9.figshare.14671980 
mkdir -p data/sumstats/raw/01_an_watson_2019 && cd data/sumstats/raw/01_an_watson_2019
wget --content-disposition https://figshare.com/ndownloader/files/28169268 # https://figshare.com/articles/dataset/an2019/14671980?file=28169268
wget --content-disposition https://figshare.com/ndownloader/files/28169271 # https://figshare.com/articles/dataset/an2019/14671980?file=28169271
cd -

sumstats="data/sumstats/raw/01_an_watson_2019/pgcAN2.2019-07.vcf.tsv.gz"
outFile="data/sumstats/harmonized/01_an_watson_2019"
skipLines=70
sep="auto"
colsIn="CHROM,POS,ID,ALT,REF,IMPINFO,BETA,SE,PVAL,NCAS,NCON"
colsOut="CHR,BP,ID,A1,A2,INFO,BETA,SE,P,Nca,Nco"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=TRUE
compN=TRUE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# PGC_MDD_2018 (Wray et al., https://pubmed.ncbi.nlm.nih.gov/29700475) data doi: 10.6084/m9.figshare.21655784
mkdir -p data/sumstats/raw/01_mdd_wray_2018 && cd data/sumstats/raw/01_mdd_wray_2018
wget --content-disposition https://figshare.com/ndownloader/files/28169502 # https://figshare.com/articles/dataset/mdd2018/14672085?file=28169502
wget --content-disposition https://figshare.com/ndownloader/files/28169508 # https://figshare.com/articles/dataset/mdd2018/14672085?file=28169508
wget --content-disposition https://figshare.com/ndownloader/files/34427408 # no ukbb: https://figshare.com/articles/dataset/mdd2018/14672085?file=34427408

sumstats="data/sumstats/raw/01_mdd_wray_2018/MDD2018_ex23andMe.gz"
outFile="data/sumstats/harmonized/01_mdd_wray_2018"
skipLines=0
sep="auto"
colsIn="CHR,BP,SNP,A1,A2,FRQ_U_113154,INFO,OR,SE,P,Nca,Nco,Neff"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,INFO,OR,SE,P,Nca,Nco,Neff_half"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=TRUE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

sumstats="data/sumstats/raw/01_mdd_wray_2018/daner_pgc_mdd_meta_w2_no23andMe_rmUKBB.gz"
outFile="data/sumstats/harmonized/01_mdd_wray_2018_noUKB"
skipLines=0
sep="auto"
colsIn="CHR,BP,SNP,A1,A2,FRQ_U_97250,INFO,OR,SE,P,Nca,Nco,Neff_half"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,INFO,OR,SE,P,Nca,Nco,Neff"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=TRUE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# PGC_MDD_2019 (Howard et al., https://pubmed.ncbi.nlm.nih.gov/30718901) data doi: 10.7488/ds/2458
mkdir -p data/sumstats/raw/01_mdd_howard_2019 && cd data/sumstats/raw/01_mdd_howard_2019 
wget --content-disposition https://datashare.ed.ac.uk/download/DS_10283_3203.zip
unzip DS_10283_3203.zip
cd -

sumstats="data/sumstats/raw/01_mdd_howard_2019/PGC_UKB_depression_genome-wide.txt"
outFile="data/sumstats/harmonized/01_mdd_howard_2019"
skipLines=0
sep="auto"
colsIn="MarkerName,A1,A2,LogOR,StdErrLogOR,P"
colsOut="ID,A1,A2,BETA,SE,P"
halveNeff=FALSE
setN=-1
setNca=170756
setNco=329443
compNeff_half=TRUE
compN=TRUE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# PGC_PTSD_2019 (Nievergelt et al., 10.1038/s41467-019-12576-w) data doi:
mkdir -p data/sumstats/raw/01_ptsd_nievergelt_2019 && cd data/sumstats/raw/01_ptsd_nievergelt_2019
wget --content-disposition https://figshare.com/ndownloader/files/28169610 # https://figshare.com/articles/dataset/ptsd2019/14672133?file=28169610
wget --content-disposition https://figshare.com/ndownloader/files/28169727 # https://figshare.com/articles/dataset/ptsd2019/14672133?file=28169727
cd -
sumstats="data/sumstats/raw/01_ptsd_nievergelt_2019/pts_eur_freeze2_overall.results.gz"
outFile="data/sumstats/harmonized/01_ptsd_nievergelt_2019"
skipLines=0
sep="auto"
colsIn="CHR,BP,SNP,A1,A2,FRQ_U_151447,INFO,OR,SE,P,Nca,Nco,Neff"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,INFO,OR,SE,P,Nca,Nco,Neff_half"
halveNeff=TRUE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=TRUE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# PGC_SCZ_2022 (Trubetskoy et al., 10.1038/s41586-022-04434-5) data doi: 10.6084/m9.figshare.19426775
mkdir -p data/sumstats/raw/01_scz_trubetskoy_2022 && cd data/sumstats/raw/01_scz_trubetskoy_2022
wget --content-disposition https://figshare.com/ndownloader/files/34865091 # https://figshare.com/articles/dataset/scz2022/19426775?file=34865091
wget --content-disposition https://figshare.com/ndownloader/files/34517828 # https://figshare.com/articles/dataset/scz2022/19426775?file=34517828
cd -
sumstats="data/sumstats/raw/01_scz_trubetskoy_2022/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"
outFile="data/sumstats/harmonized/01_scz_trubetskoy_2022"
skipLines=73
sep="auto"
colsIn="CHROM,POS,ID,A1,A2,FCON,IMPINFO,BETA,SE,PVAL,NCAS,NCON,NEFF"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,INFO,BETA,SE,P,Nca,Nco,Neff"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=TRUE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# =====================
# === substance use ===
# =====================

# drinks per week (Saunders et al., https://doi.org/10.1038/s41586-022-05477-4)
mkdir -p data/sumstats/raw/02_dpw_saunders_2022 && cd data/sumstats/raw/02_dwp_saunders_2022
wget --content-disposition https://conservancy.umn.edu/bitstream/handle/11299/241912/README.txt
wget https://conservancy.umn.edu/bitstream/handle/11299/241912/EUR_stratified.zip
unzip EUR_stratified.zip
rm -f EUR_stratified.zip
rm -rf __MACOSX
mkdir other
mv GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR.txt.gz other/
mv GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR.txt.gz other/
cd -
sumstats="data/sumstats/raw/02_dpw_saunders_2022/GSCAN_DrnkWk_2022_GWAS_SUMMARY_STATS_EUR.txt.gz"
outFile="data/sumstats/harmonized/02_dpw_saunders_2022"
skipLines=0
sep="auto"
colsIn="CHR,POS,RSID,EFFECT_ALLELE,OTHER_ALLELE,AF_1000G,BETA,SE,P,N"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# cigarettes per day (Saunders et al., https://doi.org/10.1038/s41586-022-05477-4)
mkdir -p data/sumstats/raw/02_cpd_saunders_2022 && cd data/sumstats/raw/02_cpd_saunders_2022
mv ../02_dpw_saunders_2022/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt.gz .
cd -
sumstats="data/sumstats/raw/02_cpd_saunders_2022/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt.gz"
outFile="data/sumstats/harmonized/02_cpd_saunders_2022"
skipLines=0
sep="auto"
colsIn="CHR,POS,RSID,EFFECT_ALLELE,OTHER_ALLELE,AF_1000G,BETA,SE,P,N"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# smoking initiation (Saunders et al., https://doi.org/10.1038/s41586-022-05477-4)
mkdir -p data/sumstats/raw/02_smkInit_saunders_2022 && cd data/sumstats/raw/02_smkInit_saunders_2022
mv ../02_alc_saunders_2022/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.txt.gz .
sumstats="data/sumstats/raw/02_smkInit_saunders_2022/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.txt.gz"
outFile="data/sumstats/harmonized/02_smkInit_saunders_2022"
skipLines=0
sep="auto"
colsIn="CHR,POS,RSID,EFFECT_ALLELE,OTHER_ALLELE,AF_1000G,BETA,SE,P,N"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# cannabis use disorder (Johnson et al., https://doi.org/10.1016/S2215-0366(20)30339-4) data doi: 10.6084/m9.figshare.14842692
mkdir -p data/sumstats/raw/02_cud_johnson_2020 && cd data/sumstats/raw/02_cud_johnson_2020
wget --content-disposition https://figshare.com/ndownloader/files/28570821 # https://figshare.com/articles/dataset/sud2020-cud/14842692?file=28570821 
wget --content-disposition https://figshare.com/ndownloader/files/28570824 # https://figshare.com/articles/dataset/sud2020-cud/14842692?file=28570824
cd -
sumstats="data/sumstats/raw/02_cud_johnson_2020/CUD_EUR_full_public_11.14.2020.gz"
outFile="data/sumstats/raw/02_cud_johnson_2020"
skipLines=0
sep="auto"
colsIn="CHR,BP,SNP,A1,A2,Z,P,N_CAS,N_CON,N"
colsOut="CHR,BP,ID,A1,A2,Z,P,Nca,Nco,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=TRUE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# coffee consumption (Zhong et al., https://doi.org/10.1093/hmg/ddz061)
mkdir -p data/sumstats/raw/02_cof_zhong_2019 && cd data/sumstats/raw/02_cof_zhong_2019
wget --content-disposition "https://prism.northwestern.edu/records/ze65d-hcg18/files/coffee.assess.nobmi.gz?download=1"
cd -
sumstats="data/sumstats/raw/02_cof_zhong_2019/coffee.assess.nobmi.gz"
outFile="data/sumstats/harmonized/02_cof_zhong_2019"
skipLines=0
sep="auto"
colsIn="SNP,ID,ALT,REF,ALT_FREQ,MACH_R2,BETA,SE,P,OBS_CT"
colsOut="CHR_BP,ID,A1,A2,A1_FREQ,INFO,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# ====================
# === neurological ===
# ====================

# PGC_Alzheimer's_2021 (Wightman et al., https://doi.org/10.1038/s41588-021-00921-z)
mkdir -p data/sumstats/raw/03_ad_wightman_2021 && cd data/sumstats/raw/03_ad_wightman_2021
wget --content-disposition https://ctg.cncr.nl/documents/p1651/PGCALZ2sumstatsExcluding23andMe.txt.gz
cd -
sumstats="data/sumstats/raw/03_ad_wightman_2021/PGCALZ2sumstatsExcluding23andMe.txt.gz"
outFile="data/sumstats/harmonized/03_ad_wightman_2021"
skipLines=0
sep="auto"
colsIn="chr,PosGRCh37,testedAllele,otherAllele,z,p,N"
colsOut="CHR,BP,A1,A2,Z,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=FALSE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=TRUE
kgpFiles='data/1kgp/v5a/ALL.chr\$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.meta.vcf.gz'
keepKGPfreq=TRUE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# generalized epilepsy 2022 (ILAE consortium, https://doi.org/10.1101/2022.06.08.22276120)
mkdir -p data/sumstats/raw/03_epiGen_ilae_2022 && cd data/sumstats/raw/03_epiGen_ilae_2022
wget --content-disposition http://www.epigad.org/download/final_sumstats.zip # https://www.epigad.org/
unzip -p final_sumstats.zip ILAE3_Caucasian_GGE_final.tbl > ILAE3_Caucasian_GGE_final.tbl
unzip -p final_sumstats.zip Readme.txt > Readme.txt
cd -
sumstats="data/sumstats/raw/03_epiGen_ilae_2022/ILAE3_Caucasian_GGE_final.tbl"
outFile="data/sumstats/harmonized/03_epiGen_ilae_2022"
skipLines=0
sep="tab"
colsIn="CHR,BP,MarkerName,Allele1,Allele2,Freq1,Z-score,P-value,Effective_N"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,Z,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# focal epilepsy 2022 (ILAE consortium, https://doi.org/10.1101/2022.06.08.22276120)
mkdir -p data/sumstats/raw/03_epiFoc_ilae_2022 && cd data/sumstats/raw/03_epiFoc_ilae_2022
unzip -p ../03_epiGen_ilae_2022/final_sumstats.zip ILAE3_Caucasian_focal_epilepsy_final.tbl > ILAE3_Caucasian_focal_epilepsy_final.tbl
unzip -p ../03_epiGen_ilae_2022/final_sumstats.zip Readme.txt > Readme.txt
cd -
sumstats="data/sumstats/raw/03_epiFoc_ilae_2022/ILAE3_Caucasian_focal_epilepsy_final.tbl"
outFile="data/sumstats/harmonized/03_epiFoc_ilae_2022"
skipLines=0
sep="tab"
colsIn="CHR,BP,MarkerName,Allele1,Allele2,Freq1,Z-score,P-value,Effective_N"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,Z,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# amyotrophic lateral sclerosis (van Rheenen et al., 10.1016/j.neuron.2018.02.027)
mkdir -p data/sumstats/raw/03_als_vanRheenen_2021 && cd data/sumstats/raw/03_als_vanRheenen_2021
wget --content-disposition http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027164/GCST90027164_buildGRCh37.tsv.gz
cd -
sumstats="data/sumstats/raw/03_als_vanRheenen_2021/GCST90027164_buildGRCh37.tsv.gz"
outFile="data/sumstats/harmonized/03_als_vanRheenen_2021"
skipLines=0
sep="auto"
colsIn="chromosome,base_pair_location,rsid,effect_allele,other_allele,effect_allele_frequency,beta,standard_error,p_value,N_effective"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=FALSE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# stroke (Mishra et al, 10.1038/s41586-022-05165-3)
mkdir -p data/sumstats/raw/03_stroke_mishra_2022 && cd data/sumstats/raw/03_stroke_mishra_2022
wget --content-disposition http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/GCST90104539/README.txt
wget --content-disposition http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/GCST90104539/GCST90104539_buildGRCh37.tsv.gz
cd -
sumstats="data/sumstats/raw/03_stroke_mishra_2022/GCST90104539_buildGRCh37.tsv.gz"
outFile="data/sumstats/harmonized/03_stroke_mishra_2022"
skipLines=0
sep="auto"
colsIn="chromosome,base_pair_location,effect_allele,other_allele,effect_allele_frequency,beta,standard_error,p_value"
colsOut="CHR,BP,A1,A2,A1_FREQ,BETA,SE,P"
halveNeff=FALSE
setN=-1
setNca=73652
setNco=1234808
compNeff_half=TRUE
compN=TRUE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=TRUE
kgpFiles='data/1kgp/v5a/ALL.chr\$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.meta.vcf.gz'
keepKGPfreq=TRUE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# ===================
# === personality ===
# ===================

# well-being spectrum (Baselmans et al., 10.1038/s41588-018-0320-8): files https://surfdrive.surf.nl/files/index.php/s/Ow1qCDpFT421ZOO
mkdir -p data/sumstats/raw/04_well_baselmans_2019 && cd data/sumstats/raw/04_well_baselmans_2019
wget --content-disposition "https://surfdrive.surf.nl/files/index.php/s/Ow1qCDpFT421ZOO/download?path=%2FMultivariate_GWAMA_sumstats%2FN_GWAMA&files=read.me"
wget --content-disposition "https://surfdrive.surf.nl/files/index.php/s/Ow1qCDpFT421ZOO/download?path=%2FMultivariate_GWAMA_sumstats%2FN_GWAMA&files=Read.me.UPDATE.docx"
wget --content-disposition "https://surfdrive.surf.nl/files/index.php/s/Ow1qCDpFT421ZOO/download?path=%2FMultivariate_GWAMA_sumstats%2FN_GWAMA&files=N_GWAMA_WBspectrum_no23andME.txt_04072019.txt.gz"
cd -
sumstats="data/sumstats/raw/04_well_baselmans_2019/N_GWAMA_WBspectrum_no23andME.txt_04072019.txt.gz"
outFile="data/sumstats/harmonized/04_well_baselmans_2019"
skipLines=0
sep="auto"
colsIn="CHR,BP,MarkerName,A1,A2,EAF,Z,PVAL,N"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,Z,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# neuroticism (Baselmans et al., 10.1038/s41588-018-0320-8): files https://surfdrive.surf.nl/files/index.php/s/Ow1qCDpFT421ZOO
mkdir -p data/sumstats/raw/04_neur_baselmans_2019 && cd data/sumstats/raw/04_neur_baselmans_2019
wget --content-disposition "https://surfdrive.surf.nl/files/index.php/apps/richdocuments/public?fileId=6754782050&shareToken=Ow1qCDpFT421ZOO"
wget --content-disposition "https://surfdrive.surf.nl/files/index.php/s/Ow1qCDpFT421ZOO/download?path=%2FMultivariate_GWAMA_sumstats%2FUnivariate&files=NEUno23and_18022020.txt.gz"
cd -
sumstats="data/sumstats/raw/04_neur_baselmans_2019/NEUno23and_18022020.txt.gz"
outFile="data/sumstats/harmonized/04_neur_baselmans_2019"
skipLines=0
sep="auto"
colsIn="CHR,BP,MarkerName,A1,A2,EAF,Z,PVAL,N"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,Z,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol="Z"
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# loneliness (Day et al., https://doi.org/10.1038/s41467-018-04930-1): files https://surfdrive.surf.nl/files/index.php/s/Ow1qCDpFT421ZOO
mkdir -p data/sumstats/raw/04_lone_day_2018 && cd data/sumstats/raw/04_lone_day_2018
wget --content-disposition "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006923/README.txt"
wget --content-disposition "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006924/MTAG_results.txt.gz"
cd -
sumstats="data/sumstats/raw/04_lone_day_2018/MTAG_results.txt.gz"
outFile="data/sumstats/harmonized/04_lone_day_2018"
skipLines=0
sep="auto"
colsIn="chr,bpos,snpid,a1,a2,freq,mtag_beta,mtag_se,mtag_pval,n"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# general risk-tolerance (Karlsson Linner et al., https://doi.org/10.1038/s41588-018-0309-3)
mkdir -p data/sumstats/raw/04_risk_karlssonlinner_2019 && cd data/sumstats/raw/04_risk_karlssonlinner_2019
wget -O README_RISK.txt "https://ssgac.s3.amazonaws.com/README_RISK.txt?response-content-disposition=attachment&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIA2FXDXDVZD43Y6ZGT%2F20230424%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20230424T144503Z&X-Amz-Expires=30&X-Amz-SignedHeaders=host&X-Amz-Signature=40bb8986bb5b81a6161e2f46f86507c58f6f823351f97cc097f371e02059949f"
wget -O RISK_GWAS_MA_UKB+replication.txt "https://ssgac.s3.amazonaws.com/RISK_GWAS_MA_UKB%2Breplication.txt?response-content-disposition=attachment&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIA2FXDXDVZD43Y6ZGT%2F20230421%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20230421T191443Z&X-Amz-Expires=30&X-Amz-SignedHeaders=host&X-Amz-Signature=5bed2399746c6086c591aee6a488b5bd283cc8ebde3c63f538fbf7242753a86a"
wget -O RISK_GWAS_MA_UKB+23andMe+replication.txt "https://ssgac.s3.amazonaws.com/RISK_GWAS_MA_UKB%2B23andMe%2Breplication.txt?response-content-disposition=attachment&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIA2FXDXDVZD43Y6ZGT%2F20230421%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20230421T191807Z&X-Amz-Expires=30&X-Amz-SignedHeaders=host&X-Amz-Signature=fb378ed3de54839929489abb6fbad9dbf1a0a426853696198a0367b145fec8a1"
cd -
sumstats="data/sumstats/raw/04_risk_karlssonlinner_2019/RISK_GWAS_MA_UKB+replication.txt"
outFile="data/sumstats/harmonized/04_risk_karlssonlinner_2019"
skipLines=0
sep="auto"
colsIn="CHR,POS,MarkerName,A1,A2,EAF_A1,Beta,SE,Pval"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,BETA,SE,P"
halveNeff=FALSE
setN=939908
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# =============
# === sleep ===
# =============

# insomnia (Jansen et al., 10.1038/s41588-018-0333-3)
mkdir -p data/sumstats/raw/05_ins_jansen_2019 && cd data/sumstats/raw/05_ins_jansen_2019
wget --content-disposition https://ctg.cncr.nl/documents/p1651/Insomnia_sumstats_Jansenetal.readme.txt
wget --content-disposition https://ctg.cncr.nl/documents/p1651/Insomnia_sumstats_Jansenetal.txt.gz
cd -
sumstats="data/sumstats/raw/05_ins_jansen_2019/Insomnia_sumstats_Jansenetal.txt.gz"
outFile="data/sumstats/harmonized/05_ins_jansen_2019"
skipLines=0
sep="auto"
colsIn="CHR,BP,SNP,A1,A2,INFO,OR,SE,P,N,MAF"
colsOut="CHR,BP,ID,A1,A2,INFO,OR,SE,P,N,MAF"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=FALSE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# chronotype (Jones, 2019, 10.1038/s41467-018-08259-7)
mkdir -p data/sumstats/raw/05_chron_jones_2019 && cd data/sumstats/raw/05_chron_jones_2019
wget --content-disposition https://s3.amazonaws.com/broad-portal-resources/sleep/chronotype_raw_README.txt
wget --content-disposition https://personal.broadinstitute.org/ryank/chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt.gz
sumstats="data/sumstats/raw/05_chron_jones_2019/chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt.gz"
outFile="data/sumstats/harmonized/05_chron_jones_2019"
skipLines=0
sep="auto"
colsIn="CHR,BP,SNP,ALLELE1,ALLELE0,A1FREQ,INFO,BETA,SE,P_BOLT_LMM"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,INFO,BETA,SE,P"
halveNeff=FALSE
setN=449734
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# daytime sleepiness (Wang et al., 10.1038/s41467-019-11456-7)
mkdir -p data/sumstats/raw/05_ds_wang_2019 && cd data/sumstats/raw/05_ds_wang_2019
wget --content-disposition https://broad-portal-resources.s3.amazonaws.com/sleep/Saxena.fullUKBB.DaytimeSleepiness.README.txt
wget --content-disposition https://personal.broadinstitute.org/ryank/Saxena.fullUKBB.DaytimeSleepiness.sumstats.zip
unzip Saxena.fullUKBB.DaytimeSleepiness.sumstats.zip
cd -
sumstats="data/sumstats/raw/05_ds_wang_2019/Saxena.fullUKBB.DaytimeSleepiness.sumstats.txt"
outFile="data/sumstats/harmonized/05_ds_wang_2019"
skipLines=0
sep="auto"
colsIn="CHR,BP,SNP,ALLELE1,ALLELE0,A1FREQ,INFO,BETA,SE,P"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,INFO,BETA,SE,P"
halveNeff=FALSE
setN=452071
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# sleep duration (Dashti et al., https://doi.org/10.1038/s41467-019-08917-4)
mkdir -p data/sumstats/raw/05_sd_dashti_2019 && cd data/sumstats/raw/05_sd_dashti_2019
wget --content-disposition https://s3.amazonaws.com/broad-portal-resources/sleep/Saxena_fullUKBB_Sleepduration_summary_stats_README
wget --content-disposition https://personal.broadinstitute.org/ryank/sleepdurationsumstats.txt.zip
unzip sleepdurationsumstats.txt.zip
cd -
sumstats="data/sumstats/raw/05_sd_dashti_2019/sleepdurationsumstats.txt"
outFile="data/sumstats/harmonized/05_sd_dashti_2019"
skipLines=0
sep="auto"
colsIn="CHR,BP,SNP,ALLELE1,ALLELE0,A1FREQ,INFO,BETA_SLEEPDURATION,SE_SLEEPDURATION,P_SLEEPDURATION"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,INFO,BETA,SE,P"
halveNeff=FALSE
setN=446118
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# =================
# === cognition ===
# =================

# educational attainment (Okbay et al., https://doi.org/10.1038/s41588-022-01016-z) - data resource: https://thessgac.com/papers/14
mkdir -p data/sumstats/raw/06_ea_okbay_2022 && cd data/sumstats/raw/06_ea_okbay_2022
wget -O ReadMe_EA4.txt "https://ssgac.s3.amazonaws.com/ReadMe_EA4.txt?response-content-disposition=attachment&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIA2FXDXDVZD43Y6ZGT%2F20230421%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20230421T175932Z&X-Amz-Expires=30&X-Amz-SignedHeaders=host&X-Amz-Signature=ae2de949782b3aed03e0d3490198075a478d21093f33a23df2167eaf8811d1d6"
wget -O EA4_additive_excl_23andMe.txt.gz "https://ssgac.s3.amazonaws.com/EA4_additive_excl_23andMe.txt.gz?response-content-disposition=attachment&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIA2FXDXDVZD43Y6ZGT%2F20230421%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20230421T180040Z&X-Amz-Expires=30&X-Amz-SignedHeaders=host&X-Amz-Signature=297d8c3cc3a0dea3800dfee5c0b3aa0df397b911f1f15b98c722392ed8e715d5"
cd -
sumstats="data/sumstats/raw/06_ea_okbay_2022/EA4_additive_excl_23andMe.txt.gz"
outFile="data/sumstats/harmonized/06_ea_okbay_2022"
skipLines=0
sep="auto"
colsIn="Chr,BP,rsID,Effect_allele,Other_allele,EAF_HRC,Beta,SE,P"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,BETA,SE,P"
halveNeff=FALSE
setN=765283
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# intelligence (Savage et al., https://doi.org/10.1038/s41588-018-0152-6)
mkdir -p data/sumstats/raw/06_int_savage_2018 && cd data/sumstats/raw/06_int_savage_2018
wget --content-disposition http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006250/SavageJansen_IntMeta_sumstats.zip
unzip SavageJansen_IntMeta_sumstats.zip
cd -
sumstats="data/sumstats/raw/06_int_savage_2018/sumstats/SavageJansen_2018_intelligence_metaanalysis.txt"
outFile="data/sumstats/harmonized/06_int_savage_2018"
skipLines=0
sep="auto"
colsIn="CHR,POS,SNP,A1,A2,EAF_HRC,minINFO,stdBeta,SE,P,N_analyzed"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,INFO,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# cognitive performance (Lee et al. https://doi.org/10.1038/s41588-018-0147-3)
mkdir -p data/sumstats/raw/06_cp_lee_2018 && cd data/sumstats/raw/06_cp_lee_2018
wget -O ReadMe_EA3.txt "https://ssgac.s3.amazonaws.com/README_EA3.txt?response-content-disposition=attachment&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIA2FXDXDVZD43Y6ZGT%2F20230421%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20230421T181947Z&X-Amz-Expires=30&X-Amz-SignedHeaders=host&X-Amz-Signature=cbef67064fd3b6089ab6414f4e2e0d210d62bf05ba4ae91153872fb6fe9ac6d0"
wget -O GWAS_CP_all.txt "https://ssgac.s3.amazonaws.com/GWAS_CP_all.txt?response-content-disposition=attachment&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIA2FXDXDVZD43Y6ZGT%2F20230421%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20230421T182013Z&X-Amz-Expires=30&X-Amz-SignedHeaders=host&X-Amz-Signature=9d5b0d6b6515b62d86103fe8f6bc3f8bda256e7d966968060c0d42d266e738de"
cd -
sumstats="data/sumstats/raw/06_cp_lee_2018/GWAS_CP_all.txt"
outFile="data/sumstats/harmonized/06_cp_lee_2018"
skipLines=0
sep="auto"
colsIn="CHR,POS,MarkerName,A1,A2,EAF,Beta,SE,Pval"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,BETA,SE,P"
halveNeff=FALSE
setN=257828
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# reaction time (Davies, https://doi.org/10.1038/s41467-018-04362-x)
mkdir -p data/sumstats/raw/06_rt_davies_2018 && cd data/sumstats/raw/06_rt_davies_2018
wget --content-disposition http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006268/Davies2018_UKB_RT_summary_results_29052018.txt
cd -
sumstats="data/sumstats/raw/06_rt_davies_2018/Davies2018_UKB_RT_summary_results_29052018.txt"
outFile="data/sumstats/harmonized/06_rt_davies_2018"
skipLines=0
sep="auto"
colsIn="CHR,BP,MarkerName,Effect_allele,Other_allele,Beta,P"
colsOut="CHR,BP,ID,A1,A2,BETA,P"
halveNeff=FALSE
setN=330069
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=FALSE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# memory performance (Davies, https://doi.org/10.1038/mp.2016.45)
mkdir -p data/sumstats/raw/06_mem_davies_2016 && cd data/sumstats/raw/06_mem_davies_2016
wget --content-disposition http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003497/Davies2016_README
wget --content-disposition http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003497/Davies2016_UKB_Memory_summary_results_22072016.txt
cd -
sumstats="data/sumstats/raw/06_mem_davies_2016/Davies2016_UKB_Memory_summary_results_22072016.txt"
outFile="data/sumstats/harmonized/06_mem_davies_2016"
skipLines=0
sep="auto"
colsIn="Chromosome,Position,Markername,Effect_allele,Other_allele,Beta,P-value"
colsOut="CHR,BP,ID,A1,A2,BETA,P"
halveNeff=FALSE
setN=112067
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=FALSE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# ======================
# === anthropometric ===
# ======================

# height (Yengo et al. 2022, )
mkdir -p data/sumstats/raw/07_height_yengo_2022 && cd data/sumstats/raw/07_height_yengo_2022
wget --content-disposition "https://www.joelhirschhornlab.org/_files/ugd/1d0101_1cefe701e51246d29a2ec85f7f76db29.txt"
wget --content-disposition "https://504394d8-624a-4827-9f25-95a83cd9675a.filesusr.com/archives/1d0101_394c2d3120ba4d0a9f6326ed56ff8854.gz?dn=GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz"
cd -
sumstats="data/sumstats/raw/07_height_yengo_2022/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz"
outFile="data/sumstats/harmonized/07_height_yengo_2022"
skipLines=0
sep="auto"
colsIn="CHR,POS,RSID,EFFECT_ALLELE,OTHER_ALLELE,EFFECT_ALLELE_FREQ,BETA,SE,P,N"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# bmi (Yengo et al. 2018, 10.1093/hmg/ddy271)
mkdir -p data/sumstats/raw/07_bmi_yengo_2018 && cd data/sumstats/raw/07_bmi_yengo_2018
wget --content-disposition "https://portals.broadinstitute.org/collaboration/giant/images/0/01/README_summary_statistics_Yengo_et_al_2018.txt"
wget --content-disposition "https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz"
cd -
sumstats="data/sumstats/raw/07_bmi_yengo_2018/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz"
outFile="data/sumstats/harmonized/07_bmi_yengo_2018"
skipLines=0
sep="auto"
colsIn="CHR,POS,SNP,Tested_Allele,Other_Allele,Freq_Tested_Allele_in_HRS,BETA,SE,P,N"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# waist-hip-ratio (Pulit et al., 2019, https://doi.org/10.1093/hmg/ddy327) all files - including those without BMI adjustment - at https://doi.org/10.5281/zenodo.1251813
mkdir -p data/sumstats/raw/07_whr_pulit_2019 && cd data/sumstats/raw/07_whr_pulit_2019
wget --content-disposition "https://portals.broadinstitute.org/collaboration/giant/images/a/a4/README_summary_statistics_Pulit_et_al_2018.txt"
wget --content-disposition "https://portals.broadinstitute.org/collaboration/giant/images/6/6e/Whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt.gz" # wget --content-disposition "https://zenodo.org/record/1251813/files/whr.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1"
cd -
sumstats="data/sumstats/raw/07_whr_pulit_2019/Whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt.gz"
outFile="data/sumstats/harmonized/07_whr_pulit_2019"
skipLines=0
sep="auto"
colsIn="CHR,POS,SNP,Tested_Allele,Other_Allele,Freq_Tested_Allele,INFO,BETA,SE,P,N"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,INFO,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# =======================
# === cardio-vascular ===
# =======================

# coronary artery disease (Aragam et al. 2022, 10.1038/s41588-022-01233-6)
mkdir -p data/sumstats/raw/08_cad_aragam_2022 && cd data/sumstats/raw/08_cad_aragam_2022
wget --content-disposition "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132314/README.txt"
wget --content-disposition "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132314/GCST90132314_buildGRCh37.tsv"
cd -
sumstats="data/sumstats/raw/08_cad_aragam_2022/GCST90132314_buildGRCh37.tsv"
outFile="data/sumstats/harmonized/08_cad_aragam_2022"
skipLines=0
sep="auto"
colsIn="chromosome,base_pair_location,effect_allele,other_allele,effect_allele_frequency,beta,standard_error,p_value,n"
colsOut="CHR,BP,A1,A2,A1_FREQ,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=TRUE
kgpFiles='data/1kgp/v5a/ALL.chr\$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.meta.vcf.gz'
keepKGPfreq=TRUE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}
	
# systolic blood pressure (Evangelou et al. 2018, 10.1038/s41588-018-0205-x)
mkdir -p data/sumstats/raw/08_sbp_evangelou_2018 && cd data/sumstats/raw/08_sbp_evangelou_2018
wget --content-disposition "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006624/README_Evangelou_30224653.txt"
wget --content-disposition "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006624/Evangelou_30224653_SBP.txt.gz"
cd -
sumstats="data/sumstats/raw/08_sbp_evangelou_2018/Evangelou_30224653_SBP.txt.gz"
outFile="data/sumstats/harmonized/08_sbp_evangelou_2018"
skipLines=0
sep="auto"
colsIn="MarkerName,Allele1,Allele2,Freq1,Effect,StdErr,P,TotalSampleSize,N_effective"
colsOut="CHR_BP,A1,A2,A1_FREQ,BETA,SE,P,N,Neff_half"
halveNeff=TRUE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=TRUE
kgpFiles='data/1kgp/v5a/ALL.chr\$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.meta.vcf.gz'
keepKGPfreq=TRUE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# diastolic blood pressure (Evangelou et al. 2018, 10.1038/s41588-018-0205-x)
mkdir -p data/sumstats/raw/08_dbp_evangelou_2018 && cd data/sumstats/raw/08_dbp_evangelou_2018
wget --content-disposition "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006630/README_Evangelou_30224653.txt"
wget --content-disposition "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006630/Evangelou_30224653_DBP.txt.gz"
cd -
sumstats="data/sumstats/raw/08_dbp_evangelou_2018/Evangelou_30224653_DBP.txt.gz"
outFile="data/sumstats/harmonized/08_dbp_evangelou_2018"
skipLines=0
sep="auto"
colsIn="MarkerName,Allele1,Allele2,Freq1,Effect,StdErr,P,TotalSampleSize,N_effective"
colsOut="CHR_BP,A1,A2,A1_FREQ,BETA,SE,P,N,Neff_half"
halveNeff=TRUE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=TRUE
kgpFiles='data/1kgp/v5a/ALL.chr\$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.meta.vcf.gz'
keepKGPfreq=TRUE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# myocardial infarction (Hartiala et al. 2021, https://doi.org/10.1093/eurheartj/ehaa1040) N according to https://www.ebi.ac.uk/gwas/studies/GCST011365
mkdir -p data/sumstats/raw/08_myocardial_hartiala_2021 && cd data/sumstats/raw/08_myocardial_hartiala_2021
wget --content-disposition http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST011001-GCST012000/GCST011365/GCST011365_buildGRCh37.tsv
cd -
sumstats="data/sumstats/raw/08_myocardial_hartiala_2021/GCST011365_buildGRCh37.tsv"
outFile="data/sumstats/harmonized/08_myocardial_hartiala_2021"
skipLines=0
sep="auto"
colsIn="chromosome,base_pair_location,variant_id,effect_allele,other_allele,effect_allele_frequency,beta,standard_error,p_value"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,BETA,SE,P"
halveNeff=FALSE
setN=639221
setNca=61505
setNco=577716
compNeff_half=TRUE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# Type-2-Diabetes (Xue et al. 2018 https://doi.org/10.1038/s41467-018-04951-w)
mkdir -p data/sumstats/raw/08_diabetes_xue_2018 && cd data/sumstats/raw/08_diabetes_xue_2018
wget --content-disposition "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006867/Xue_et_al_T2D_META_Nat_Commun_2018.pdf"
wget --content-disposition "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006867/Xue_et_al_T2D_META_Nat_Commun_2018.gz"
cd -
sumstats="data/sumstats/raw/08_diabetes_xue_2018/Xue_et_al_T2D_META_Nat_Commun_2018.gz"
outFile="data/sumstats/harmonized/08_diabetes_xue_2018"
skipLines=0
sep="auto"
colsIn="CHR,BP,SNP,A1,A2,frq_A1,b,se,P,N"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# =======================================
# === Munge harmonized files for LDSC ===
# =======================================

# settings
conda activate envs/ldsc
targetDir=data/sumstats/munged
mkdir -p ${targetDir}
LDsnplist=/fast/software/ldsc/resources/w_hm3.noMHC.snplist
files=$(ls data/sumstats/harmonized/*.gz)

# loop over files
for i in ${files}; do
    outFile=$(echo ${i} | sed 's%.*/%%g' | sed 's/.gz//g')
    munge_sumstats.py \
        --sumstats ${i} \
        --merge-alleles "${LDsnplist}" \
        --snp ID \
        --out "${targetDir}/${outFile}"
    chmod -R 750 ${targetDir}/${outFile}*
done



