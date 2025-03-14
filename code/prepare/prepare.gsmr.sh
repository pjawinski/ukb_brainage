#!/bin/bash

# ==============================================================================================
# === Selected traits for mendelian randomization: Download, harmonize, and munge GWAS files ===
# ==============================================================================================
# take same exposure files as Zhu et al. 2018 (10.1038/s41467-017-02317-2; supplementary Table 3)

# set working directory
cd /slow/projects/ukb_brainage
conda activate envs/default
mkdir -p data/sumstats/raw
mkdir -p data/sumstats/harmonized
mkdir -p data/sumstats/munged

# EduYears (Okbay et al. 2016, https://doi.org/10.1038/nature17671), data: http://www.thessgac.org/data
mkdir -p data/sumstats/raw/mr_edu_okbay_2016 && cd data/sumstats/raw/mr_edu_okbay_2016
wget -O ReadMe.txt "https://ssgac.s3.amazonaws.com/ReadMe.txt?response-content-disposition=attachment&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIA2FXDXDVZD43Y6ZGT%2F20230523%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20230523T101517Z&X-Amz-Expires=30&X-Amz-SignedHeaders=host&X-Amz-Signature=abed11783ac7d269ccff93ead3c2fcafbe7be2b53b449a4b2081f60cae9bcc95"
wget -O EduYears_Main.txt "https://ssgac.s3.amazonaws.com/EduYears_Main.txt?response-content-disposition=attachment&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIA2FXDXDVZD43Y6ZGT%2F20230523%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20230523T101537Z&X-Amz-Expires=30&X-Amz-SignedHeaders=host&X-Amz-Signature=b201f81711e16fafdb95c1500b7b4c9d3c6a3a6d49e54de2133ec3844c16313b"
wget -O EduYears_Discovery_5000.txt "https://ssgac.s3.amazonaws.com/EduYears_Discovery_5000.txt?response-content-disposition=attachment&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIA2FXDXDVZD43Y6ZGT%2F20230523%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20230523T102145Z&X-Amz-Expires=30&X-Amz-SignedHeaders=host&X-Amz-Signature=f137c69899e959f4d8949cf5f04db068e9ce1f8dc87c9eac28af14c2e1968f9d"
cd -
sumstats="data/sumstats/raw/mr_edu_okbay_2016/EduYears_Discovery_5000.txt"
outFile="data/sumstats/harmonized/mr_eduPruned_okbay_2016"
skipLines=0
sep="auto"
colsIn="CHR,POS,A1,A2,EAF,Beta,SE,Pval"
colsOut="CHR,BP,A1,A2,A1_FREQ,BETA,SE,P"
halveNeff=FALSE
setN=293723
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
or2beta=FALSE
compBetaSe=FALSE
stdize=FALSE
addRSID=TRUE
kgpFiles='data/1kgp/v5a/ALL.chr\$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.meta.vcf.gz'
keepKGPfreq=TRUE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${or2beta} ${compBetaSe} ${stdize} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# BMI (Locke et al. 2015, https://doi.org/10.1038/nature14177), data: http://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files
mkdir -p data/sumstats/raw/mr_bmi_locke_2015 && cd data/sumstats/raw/mr_bmi_locke_2015
wget "http://portals.broadinstitute.org/collaboration/giant/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz"
cd -
sumstats="data/sumstats/raw/mr_bmi_locke_2015/SNP_gwas_mc_merge_nogc.tbl.uniq.gz"
outFile="data/sumstats/harmonized/mr_bmi_locke_2015"
skipLines=0
sep="auto"
colsIn="SNP,A1,A2,Freq1.Hapmap,b,se,p,N"
colsOut="ID,A1,A2,A1_FREQ,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
or2beta=FALSE
compBetaSe=FALSE
stdize=TRUE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${or2beta} ${compBetaSe} ${stdize} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# WHRadjBMI (Shungin et al. 2015, https://doi.org/10.1038/nature14132), data: http://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files
mkdir -p data/sumstats/raw/mr_whr_shungin_2015 && cd data/sumstats/raw/mr_whr_shungin_2015
wget "http://portals.broadinstitute.org/collaboration/giant/images/e/eb/GIANT_2015_WHRadjBMI_COMBINED_EUR.txt.gz"
zcat GIANT_2015_WHRadjBMI_COMBINED_EUR.txt.gz | awk 'BEGIN { print "MarkerName\tAllele1\tAllele2\tFreqAllele1HapMapCEU\tb\tse\tp\tN" } NR > 8 { print }' | gzip > GIANT_2015_WHRadjBMI_COMBINED_EUR.edited.txt.gz
cd -
sumstats="data/sumstats/raw/mr_whr_shungin_2015/GIANT_2015_WHRadjBMI_COMBINED_EUR.edited.txt.gz"
outFile="data/sumstats/harmonized/mr_whr_shungin_2015"
skipLines=0
sep="auto"
colsIn="MarkerName,Allele1,Allele2,FreqAllele1HapMapCEU,b,se,p,N"
colsOut="ID,A1,A2,A1_FREQ,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
or2beta=FALSE
compBetaSe=FALSE
stdize=TRUE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${or2beta} ${compBetaSe} ${stdize} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# LDL cholesterol (Willer et al. 2013, 10.1038/ng.2797)
mkdir -p data/sumstats/raw/mr_ldlc_willer_2013 && cd data/sumstats/raw/mr_ldlc_willer_2013
wget "http://csg.sph.umich.edu/abecasis/public/lipids2013/jointGwasMc_LDL.txt.gz"
cd -
sumstats="data/sumstats/raw/mr_ldlc_willer_2013/jointGwasMc_LDL.txt.gz"
outFile="data/sumstats/harmonized/mr_ldlc_willer_2013"
skipLines=0
sep="auto"
colsIn="SNP_hg19,A1,A2,Freq.A1.1000G.EUR,beta,se,P-value,N"
colsOut="CHR_BP,A1,A2,A1_FREQ,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
or2beta=FALSE
compBetaSe=FALSE
stdize=TRUE
addRSID=TRUE
kgpFiles='data/1kgp/v5a/ALL.chr\$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.meta.vcf.gz'
keepKGPfreq=TRUE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${or2beta} ${compBetaSe} ${stdize} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# HDL cholesterol (Willer et al. 2013, 10.1038/ng.2797)
mkdir -p data/sumstats/raw/mr_hdlc_willer_2013 && cd data/sumstats/raw/mr_hdlc_willer_2013
wget "http://csg.sph.umich.edu/abecasis/public/lipids2013/jointGwasMc_HDL.txt.gz"
cd -
sumstats="data/sumstats/raw/mr_hdlc_willer_2013/jointGwasMc_HDL.txt.gz"
outFile="data/sumstats/harmonized/mr_hdlc_willer_2013"
skipLines=0
sep="auto"
colsIn="SNP_hg19,A1,A2,Freq.A1.1000G.EUR,beta,se,P-value,N"
colsOut="CHR_BP,A1,A2,A1_FREQ,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
or2beta=FALSE
compBetaSe=FALSE
stdize=TRUE
addRSID=TRUE
kgpFiles='data/1kgp/v5a/ALL.chr\$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.meta.vcf.gz'
keepKGPfreq=TRUE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${or2beta} ${compBetaSe} ${stdize} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# Triglycerides (Willer et al. 2013, 10.1038/ng.2797)
mkdir -p data/sumstats/raw/mr_trigl_willer_2013 && cd data/sumstats/raw/mr_trigl_willer_2013
wget "http://csg.sph.umich.edu/abecasis/public/lipids2013/jointGwasMc_TG.txt.gz"
cd -
sumstats="data/sumstats/raw/mr_trigl_willer_2013/jointGwasMc_TG.txt.gz"
outFile="data/sumstats/harmonized/mr_trigl_willer_2013"
skipLines=0
sep="auto"
colsIn="SNP_hg19,A1,A2,Freq.A1.1000G.EUR,beta,se,P-value,N"
colsOut="CHR_BP,A1,A2,A1_FREQ,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
or2beta=FALSE
compBetaSe=FALSE
stdize=TRUE
addRSID=TRUE
kgpFiles='data/1kgp/v5a/ALL.chr\$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.meta.vcf.gz'
keepKGPfreq=TRUE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${or2beta} ${compBetaSe} ${stdize} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# Height (Wood et al. 2014, https://doi.org/10.1038/ng.3097)
mkdir -p data/sumstats/raw/mr_height_wood_2014 && cd data/sumstats/raw/mr_height_wood_2014
wget "http://portals.broadinstitute.org/collaboration/giant/images/0/01/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz"
cd -
sumstats="data/sumstats/raw/mr_height_wood_2014/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz"
outFile="data/sumstats/harmonized/mr_height_wood_2014"
skipLines=0
sep="auto"
colsIn="MarkerName,Allele1,Allele2,Freq.Allele1.HapMapCEU,b,SE,p,N"
colsOut="ID,A1,A2,A1_FREQ,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
or2beta=FALSE
compBetaSe=FALSE
stdize=TRUE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${or2beta} ${compBetaSe} ${stdize} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# PGC_SCZ_2022 (Trubetskoy et al., 10.1038/s41586-022-04434-5) data doi: 10.6084/m9.figshare.19426775
mkdir -p data/sumstats/raw/mr_scz_trubetskoy_2022 && cd data/sumstats/raw/mr_scz_trubetskoy_2022
wget --content-disposition https://figshare.com/ndownloader/files/34865091 # https://figshare.com/articles/dataset/scz2022/19426775?file=34865091
wget --content-disposition https://figshare.com/ndownloader/files/34517828 # https://figshare.com/articles/dataset/scz2022/19426775?file=34517828
cd -
sumstats="data/sumstats/raw/mr_scz_trubetskoy_2022/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"
outFile="data/sumstats/harmonized/mr_scz_trubetskoy_2022"
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
or2beta=FALSE
compBetaSe=FALSE
stdize=FALSE
or2beta=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${or2beta} ${compBetaSe} ${stdize} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# coronary artery disease (Nikpay et al. 2015, http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003116/cad.add.160614.website.txt)
mkdir -p data/sumstats/raw/mr_cad_nikpay_2015 && cd data/sumstats/raw/mr_cad_nikpay_2015
wget --content-disposition http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003116/cad.add.160614.website.txt
wget --content-disposition http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003116/cad.add.readme
cd -
sumstats="data/sumstats/raw/mr_cad_nikpay_2015/cad.add.160614.website.txt"
outFile="data/sumstats/harmonized/mr_cad_nikpay_2015"
skipLines=0
sep="auto"
colsIn="chr,bp_hg19,markername,effect_allele,noneffect_allele,effect_allele_freq,beta,se_dgc,p_dgc"
colsOut="CHR,BP,ID,A1,A2,A1_FREQ,BETA,SE,P"
halveNeff=FALSE
setN=184305 # 60,801 cases and 123,504 controls from 48 studies
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=FALSE
invertCol=FALSE
rmDuplicatesCol=FALSE
or2beta=FALSE
compBetaSe=FALSE
stdize=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${or2beta} ${compBetaSe} ${stdize} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# type-2-diabetes (Scott et al. 2017, https://doi.org/10.2337/db16-1253)
mkdir -p data/sumstats/raw/mr_t2d_scott_2017 && cd data/sumstats/raw/mr_t2d_scott_2017
wget --content-disposition http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004773/DIAGRAM_1000G_GWAS.pdf
wget --content-disposition http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004773/METAANALYSIS_DIAGRAM_SE1.txt
cd -
sumstats="data/sumstats/raw/mr_t2d_scott_2017/METAANALYSIS_DIAGRAM_SE1.txt"
outFile="data/sumstats/harmonized/mr_t2d_scott_2017"
skipLines=0
sep="auto"
colsIn="Chr:Position,Allele1,Allele2,Effect,StdErr,P-value,TotalSampleSize"
colsOut="CHR_BP,A1,A2,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=FALSE
invertCol=FALSE
rmDuplicatesCol=FALSE
or2beta=FALSE
compBetaSe=FALSE
stdize=FALSE
addRSID=TRUE
kgpFiles='data/1kgp/v5a/ALL.chr\$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.meta.vcf.gz'
keepKGPfreq=TRUE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${or2beta} ${compBetaSe} ${stdize} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# systolic blood pressure (Evangelou et al. 2018, 10.1038/s41588-018-0205-x)
mkdir -p data/sumstats/raw/mr_sbp_evangelou_2018 && cd data/sumstats/raw/mr_sbp_evangelou_2018
wget --content-disposition "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006624/README_Evangelou_30224653.txt"
wget --content-disposition "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006624/Evangelou_30224653_SBP.txt.gz"
cd -
sumstats="data/sumstats/raw/08_sbp_evangelou_2018/Evangelou_30224653_SBP.txt.gz"
outFile="data/sumstats/harmonized/mr_sbp_evangelou_2018"
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
or2beta=FALSE
compBetaSe=FALSE
stdize=TRUE
addRSID=TRUE
kgpFiles='data/1kgp/v5a/ALL.chr\$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.meta.vcf.gz'
keepKGPfreq=TRUE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${or2beta} ${compBetaSe} ${stdize} ${addRSID} ${kgpFiles} ${keepKGPfreq}

mkdir -p data/sumstats/raw/mr_sbp_noUKB_evangelou_2018 && cd data/sumstats/raw/mr_sbp_noUKB_evangelou_2018
wget --content-disposition "https://collect.qmul.ac.uk/down?t=6T6JFAI24NKREUCN/6H7JBELPF3V98MDMLJ4QM3G"
cd -
sumstats="data/sumstats/raw/mr_sbp_noUKB_evangelou_2018/ICBP_SBPsummaryResults.txt.gz"
outFile="data/sumstats/harmonized/mr_sbp_noUKB_evangelou_2018"
skipLines=0
sep="auto"
colsIn="markername,allele1,allele2,freq1,effect,stderr,pvalue,totalsamplesize"
colsOut="CHR_BP,A1,A2,A1_FREQ,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
or2beta=FALSE
compBetaSe=FALSE
stdize=TRUE
addRSID=TRUE
kgpFiles='data/1kgp/v5a/ALL.chr\$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.meta.vcf.gz'
keepKGPfreq=TRUE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${or2beta} ${compBetaSe} ${stdize} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# diastolic blood pressure (Evangelou et al. 2018, 10.1038/s41588-018-0205-x)
mkdir -p data/sumstats/raw/mr_dbp_evangelou_2018 && cd data/sumstats/raw/mr_dbp_evangelou_2018
wget --content-disposition "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006630/README_Evangelou_30224653.txt"
wget --content-disposition "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006630/Evangelou_30224653_DBP.txt.gz"
cd -
sumstats="data/sumstats/raw/mr_dbp_evangelou_2018/Evangelou_30224653_DBP.txt.gz"
outFile="data/sumstats/harmonized/mr_dbp_evangelou_2018"
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
or2beta=FALSE
compBetaSe=FALSE
stdize=TRUE
addRSID=TRUE
kgpFiles='data/1kgp/v5a/ALL.chr\$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.meta.vcf.gz'
keepKGPfreq=TRUE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${or2beta} ${compBetaSe} ${stdize} ${addRSID} ${kgpFiles} ${keepKGPfreq}

mkdir -p data/sumstats/raw/mr_dbp_noUKB_evangelou_2018 && cd data/sumstats/raw/mr_dbp_noUKB_evangelou_2018
wget --content-disposition "https://collect.qmul.ac.uk/down?t=R8MS8AU1KPK3DUGK/651D7MIP30FF9EEM9ILKEOG"
cd -
sumstats="data/sumstats/raw/mr_dbp_noUKB_evangelou_2018/ICBP_DBPsummaryResults.txt.gz"
outFile="data/sumstats/harmonized/mr_dbp_noUKB_evangelou_2018"
skipLines=0
sep="auto"
colsIn="markername,allele1,allele2,freq1,effect,stderr,pvalue,totalsamplesize"
colsOut="CHR_BP,A1,A2,A1_FREQ,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
or2beta=FALSE
compBetaSe=FALSE
stdize=TRUE
addRSID=TRUE
kgpFiles='data/1kgp/v5a/ALL.chr\$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.meta.vcf.gz'
keepKGPfreq=TRUE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${or2beta} ${compBetaSe} ${stdize} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# pulse pressure (Evangelou et al. 2018, 10.1038/s41588-018-0205-x)
mkdir -p data/sumstats/raw/mr_pp_noUKB_evangelou_2018 && cd data/sumstats/raw/mr_pp_noUKB_evangelou_2018
wget --content-disposition "https://collect.qmul.ac.uk/down?t=4HTHH1D7J1F3PV0Q/45PHT3SK8N4QUPK3S40P6C8"
cd -
sumstats="data/sumstats/raw/mr_pp_noUKB_evangelou_2018/ICBP_PPsummaryResults.txt.gz"
outFile="data/sumstats/harmonized/mr_pp_noUKB_evangelou_2018"
skipLines=0
sep="auto"
colsIn="markername,allele1,allele2,freq1,effect,stderr,pvalue,totalsamplesize"
colsOut="CHR_BP,A1,A2,A1_FREQ,BETA,SE,P,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=TRUE
invertCol=FALSE
rmDuplicatesCol=FALSE
or2beta=FALSE
compBetaSe=FALSE
stdize=TRUE
addRSID=TRUE
kgpFiles='data/1kgp/v5a/ALL.chr\$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.meta.vcf.gz'
keepKGPfreq=TRUE
Rscript code/genetics/sumstats.harmonize.R ${sumstats} ${outFile} ${skipLines} ${sep} ${colsIn} ${colsOut} ${halveNeff} ${setNca} ${setNco} ${setN} ${compNeff_half} ${compN} ${compMAF} ${invertCol} ${rmDuplicatesCol} ${or2beta} ${compBetaSe} ${stdize} ${addRSID} ${kgpFiles} ${keepKGPfreq}

# ===============================
# === prepare files for GSMR  ===
# ===============================

# specify target directory and get list of files
targetDir=data/sumstats/mr
mkdir -p ${targetDir}
files=$(ls data/sumstats/harmonized/mr_*.gz)

# settings
skipLines=0
sep="auto"
colsIn="ID,A1,A2,A1_FREQ,BETA,SE,P,N"
colsOut="SNP,A1,A2,freq,b,se,p,N"
halveNeff=FALSE
setN=-1
setNca=-1
setNco=-1
compNeff_half=FALSE
compN=FALSE
compMAF=FALSE
invertCol=FALSE
rmDuplicatesCol=SNP
or2beta=FALSE
compBetaSe=FALSE
stdize=FALSE
addRSID=FALSE
kgpFiles=FALSE
keepKGPfreq=FALSE

# loop over files
for i in ${files}; do
    outFile=$(echo "${targetDir}"/$(echo "${i}" | sed 's%.*/%%g' | sed 's/.gz//g'))
	Rscript code/genetics/sumstats.harmonize.R "${i}" "${outFile}" "${skipLines}" "${sep}" "${colsIn}" "${colsOut}" "${halveNeff}" "${setNca}" "${setNco}" "${setN}" "${compNeff_half}" "${compN}" "${compMAF}" "${invertCol}" "${rmDuplicatesCol}" "${or2beta}" "${compBetaSe}" "${stdize}" "${addRSID}" "${kgpFiles}" "${keepKGPfreq}"
    chmod -R 750 "${outFile}"*
done

# use 1kgp frequency column
i=data/sumstats/harmonized/mr_t2d_scott_2017.gz
colsIn="ID,A1,A2,KGP_A1_FREQ,BETA,SE,P,N"
outFile=$(echo $targetDir/$(echo ${i} | sed 's%.*/%%g' | sed 's/.gz//g'))
Rscript code/genetics/sumstats.harmonize.R "${i}" "${outFile}" "${skipLines}" "${sep}" "${colsIn}" "${colsOut}" "${halveNeff}" "${setNca}" "${setNco}" "${setN}" "${compNeff_half}" "${compN}" "${compMAF}" "${invertCol}" "${rmDuplicatesCol}" "${or2beta}" "${compBetaSe}" "${stdize}" "${addRSID}" "${kgpFiles}" "${keepKGPfreq}"
chmod -R 750 "${outFile}"*

