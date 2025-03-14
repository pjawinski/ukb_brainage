#!/bin/bash

# ==================================
# === install conda environments ===
# ==================================

for env in default eqtl finemap genesis gofuncr ldsc locuszoom mendelian phesant pointdensity pops power xgb; do
	if [ ! -d "envs/${env}" ]; then
		mamba env create --file envs/${env}.yml -p envs/${env}
	fi
done

# ====================================
# ==== create trait and covs files ===
# ====================================

# set working directory and activate conda environment
cd /slow/projects/ukb_brainage
conda activate envs/default

# create discovery trait files
mkdir -p data/traits
for trait in gap_gm gap_wm gap_gwm; do
	subsFile='data/genetics/chr1/imp_mri_qc/chr1_mri_qc.psam'
	varsFile="results/mri/brainage.discovery.txt"
	idCol="IID"
	vars="FID,IID,brainage_${trait}_stk"
	varsRename="FID,IID,${trait}"
	outFile="data/traits/${trait}.txt"
	Rscript code/genetics/getvars.R ${subsFile} ${varsFile} ${idCol} ${vars} ${varsRename} ${outFile} 
done

# create discovery covs file
mv "data/traits/covs.txt" "data/traits/covs.backup.txt"
subsFile='data/genetics/chr1/imp_mri_qc/chr1_mri_qc.psam'
varsFile='results/genetics/r2024.vars.pca.txt'
idCol='IID'
vars='FID,IID,sex,t1.age,t1.age2,t1.ac1,t1.ac2,array,t1.TIV,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20'
varsRename='FID,IID,sex,age,age2,ac1,ac2,array,TIV,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20'
outFile="data/traits/covs.txt"
Rscript code/genetics/getvars.R ${subsFile} ${varsFile} ${idCol} ${vars} ${varsRename} ${outFile} 

# create all-in-one replication file
subsFile='results/mri/r2024.vars.replication.txt'
varsFile='results/mri/brainage.replication.singlemodel.txt'
idCol='IID'
vars='FID,IID,pan,sex,t1_age,t1_age2,t1_ac1,t1_ac2,t1_ac3,t1_TIV,array,PanC1,PanC2,PanC3,PanC4,PanC5,PanC6,PanC7,PanC8,PanC9,PanC10,PanC11,PanC12,PanC13,PanC14,PanC15,PanC16,PanC17,PanC18,PanC19,PanC20,brainage_gap_gm_stack,brainage_gap_wm_stack,brainage_gap_gwm_stack'
varsRename='FID,IID,pan,sex,age,age2,ac1,ac2,ac3,TIV,array,PanC1,PanC2,PanC3,PanC4,PanC5,PanC6,PanC7,PanC8,PanC9,PanC10,PanC11,PanC12,PanC13,PanC14,PanC15,PanC16,PanC17,PanC18,PanC19,PanC20,gap_gm,gap_wm,gap_gwm'
outFile="data/traits/replicate.txt"
Rscript code/genetics/getvars.R ${subsFile} ${varsFile} ${idCol} ${vars} ${varsRename} ${outFile} 

# ====================================================
# === run PHESANT cross-trait association analysis ===
# ====================================================

# set working directory and activate conda environment
cd /slow/projects/ukb_brainage
conda activate envs/phesant

# settings
covs="sex,age,age2,ac1,ac2,TIV"
phenoFile="data/basket/20210205_2007685/data/ukb45233_imaging_phesant.csv" # phenoFile="data/basket/20200409_2007685/data/ukb41573_imaging_wba_phesant.csv"
phesantDir="/fast/software/PHESANT" # phesantDir="/home/groups/markett/software/PHESANT"
nparts=1
standardise="FALSE"

# run analysis
trait="gap_gm"; traitFile="data/${trait}/${trait}.txt"; covsFile="data/${trait}/covs.txt"; targetDir="results/${trait}/discovery/phesant"
	taskset -c 0-24 ./code/genetics/phesant.sh "${trait}" "${traitFile}" "${covsFile}" "${covs}" "${targetDir}" "${phenoFile}" "${phesantDir}" "${nparts}" "${standardise}"
trait="gap_wm"; traitFile="data/${trait}/${trait}.txt"; covsFile="data/${trait}/covs.txt"; targetDir="results/${trait}/discovery/phesant"
	taskset -c 25-49 ./code/genetics/phesant.sh "${trait}" "${traitFile}" "${covsFile}" "${covs}" "${targetDir}" "${phenoFile}" "${phesantDir}" "${nparts}" "${standardise}"
trait="gap_gwm"; traitFile="data/${trait}/${trait}.txt"; covsFile="data/${trait}/covs.txt"; targetDir="results/${trait}/discovery/phesant"
	taskset -c 50-74 ./code/genetics/phesant.sh "${trait}" "${traitFile}" "${covsFile}" "${covs}" "${targetDir}" "${phenoFile}" "${phesantDir}" "${nparts}" "${standardise}"

# create plots and get result summary
imaging=FALSE
ylim=60
ysteps=10
multipleTesting=both
width=8.7
height=5.0
repel_nudge_y=10
(for trait in gap_gm gap_wm gap_gwm; do (
	phesantResults="results/${trait}/discovery/phesant/phesant.output/results-combined.txt"
	targetDir="results/${trait}/discovery/phesant/"
	code/genetics/phesant.output.R "${trait}" "${phesantResults}" "${targetDir}" "${imaging}" "${multipleTesting}" "${ylim}" "${ysteps}" "${width}" "${height}" "${repel_nudge_y}" 
	) &
done)

# combine phesant plots and result summaries across traits
conda activate envs/default
traits="gap_gm,gap_wm,gap_gwm"
phesantSummary="results/gap_gm/discovery/phesant/phesant.summary.txt,results/gap_wm/discovery/phesant/phesant.summary.txt,results/gap_gwm/discovery/phesant/phesant.summary.txt"
phesantPlot="results/gap_gm/discovery/phesant/phesant.png,results/gap_wm/discovery/phesant/phesant.png,results/gap_gwm/discovery/phesant/phesant.png"
outputFile="results/combined/discovery.phewas"
code/genetics/phesant.combine.R "${traits}" "${phesantSummary}" "${phesantPlot}" "${outputFile}"

# =======================================================================
# === correlations with FreeSurfer variables (Desikan-Killiany Atlas) ===
# =======================================================================

# set working directory and activate conda environment
cd /slow/projects/ukb_brainage
conda activate envs/default

# calculate correlations
covNames="sex,age,age2,ac1,ac2,TIV"
freesurferFile="results/mri/freesurfer.tsv.gz"
(for trait in gap_gm gap_wm gap_gwm; do (
	traitFile="data/traits/${trait}.txt"
	covsFile="data/traits/covs.txt"
	targetDir="results/${trait}/discovery/surfcorr"
	code/mri/surfcorr.R "${trait}" "${traitFile}" "${covsFile}" "${covNames}" "${freesurferFile}" "${targetDir}"
	) &
done)

# combine surfcorr results
traits="gap_gm,gap_wm,gap_gwm"
surfcorrFiles="results/gap_gm/discovery/surfcorr/surfcorr.txt,results/gap_wm/discovery/surfcorr/surfcorr.txt,results/gap_gwm/discovery/surfcorr/surfcorr.txt"
outFile="results/combined/discovery.surfcorr.txt"
code/mri/surfcorr.combine.R "${traits}" "${surfcorrFiles}" "${outFile}"

# make plots
(for trait in gap_gm gap_wm gap_gwm; do (
	matlab -r "\
		workingDir = '$(pwd)';\
		ENIGMAtoolboxPath = '/fast/software/ENIGMA/matlab';\
		surfCorr = 'results/combined/discovery.surfcorr.txt';\
		mappingFileSurfArea = 'code/mri/surfplot.mapping.surfarea.txt';\
		mappingFileThickAvg = 'code/mri/surfplot.mapping.thickavg.txt';\
		mappingFileGrayVol = 'code/mri/surfplot.mapping.grayvol.txt';\
		mappingFileSubcortical = 'code/mri/surfplot.mapping.subcortical.txt';\
		rCol = '${trait}_rho';\
		pCol = '${trait}_fdr';\
		barTitle = 'Correlation (r)';\
		outFile = 'results/${trait}/discovery/surfcorr/surfplot.png';\
		colorBar = 'horizontal';\
		cblim = '0.35';\
		run code/mri/surfplot.m;\
		exit"
		) &
done
wait)

# combine surfplots for supplementum
traits="gap_gm,gap_wm,gap_gwm"
surfPlots="results/gap_gm/discovery/surfcorr/surfplot.png,results/gap_wm/discovery/surfcorr/surfplot.png,results/gap_gwm/discovery/surfcorr/surfplot.png"
outFile="results/combined/discovery.surfcorr.png"
width=10
height=8.85
Rscript code/mri/surfplot.combine.R "${traits}" "${surfPlots}" "${outFile}" "${width}" "${height}"

# combine gwm prediction accuracy, phesant, and surface plots for main article
accuracyplot="results/mri/accuracy.gwm.png"
phesantplot="results/gap_gwm/discovery/phesant/phesant.png"
surfplot="results/gap_gwm/discovery/surfcorr/surfplot.png"
outFile="results/combined/discovery.phenotypic.png"
Rscript code/mri/phenotypic.combine.R "${accuracyplot}" "${phesantplot}" "${surfplot}" "${outFile}"

# reviewer request: how predictive are individual brain measures for age?
# a) create trait file for chronological age
subsFile='data/genetics/chr1/imp_mri_qc/chr1_mri_qc.psam'
varsFile="results/mri/brainage.discovery.txt"
idCol="IID"
vars="FID,IID,t1_age"
varsRename="FID,IID,chronAge"
outFile="data/traits/chronAge.txt"
Rscript code/genetics/getvars.R ${subsFile} ${varsFile} ${idCol} ${vars} ${varsRename} ${outFile} 

	# calculate correlations
	trait="chronAge"
	covNames="sex,ac1,ac2,TIV"
	freesurferFile="results/mri/freesurfer.tsv.gz"
	traitFile="data/traits/chronAge.txt"
	covsFile="data/traits/covs.txt"
	targetDir="results/combined/chronAge"
	code/mri/surfcorr.R "${trait}" "${traitFile}" "${covsFile}" "${covNames}" "${freesurferFile}" "${targetDir}"

	# add to merged output table
	awk 'NR==FNR { if($3 < 1E-307) { $3=1E-307 }; id[$1]=$2"\t"$3; next }
		FNR==1 { print $0, "NA", "chronAge_rho", "chronAge_p"; next }
		FNR> 1 { print $0, "", id[$1] }' OFS='\t' results/combined/chronAge/surfcorr.txt results/combined/discovery.surfcorr.txt  \
		> results/combined/discovery.surfcorr.chronAge.txt


# ==========================================================
# === discovery sample: genome-wide association analysis ===
# ==========================================================

# set working directory
cd /slow/projects/ukb_brainage

# settings
covsFile="data/traits/covs.txt"
chrFileHandle="data/genetics/chr\${i}/imp_mri_qc/chr\${i}_mri_qc"
snpQC="data/genetics/meta/ukb_snp_qc.txt"
impMFI="data/genetics/ukb_imp_mfi"
maf=0.01
threads=100
for trait in gap_gm gap_wm gap_gwm; do
	traitFile="data/traits/${trait}.txt"
	targetDir="results/${trait}/discovery/gwas"
	code/genetics/gwas.sh "${trait}" "${traitFile}" "${covsFile}" "${targetDir}" "${chrFileHandle}" "${snpQC}" "${impMFI}" "${maf}" "${threads}" 
done

# ==================================================
# === discovery sample: run conditional analysis ===
# ==================================================

# set working directory and load conda environment
cd /slow/projectss/ukb_brainage
conda activate envs/default

# run conditional analysis and do positional annotation (!! do not run in parallel across traits)
subsFile='data/genetics/chr1/imp_mri_qc/chr1_mri_qc.psam'
chrFileHandle="data/genetics/chr\${i}/imp_mri_qc/bed/chr\${i}_mri_qc"
cols="ID,A1,A2,A1_FREQ,BETA,SE,P,N,CHR,BP,Z,ETA2,INFO,TYPED" # first 10 columns have the following order ID,A1,A2,A1_FREQ,BETA,SE,P,N,CHR,BP
pthresh=1E-6
threads=100
humandb="data/annovar/humandb"
refseq="data/refseq/GRCh37_latest_genomic.edit.gff.gz"
for trait in gap_gm gap_wm gap_gwm; do
	targetDir="results/${trait}/discovery/conditional"
	sumstats="results/${trait}/discovery/gwas/sumstats.txt.gz"
	sumstatsPLINK="results/${trait}/discovery/gwas/sumstats_plink.txt.gz"
	./code/genetics/conditional.sh "${subsFile}" "${targetDir}" "${chrFileHandle}" "${sumstats}" "${cols}" "${sumstatsPLINK}" "${pthresh}" "${threads}" "${humandb}" "${refseq}"
done

# =================================================
# === discovery sample: run ld score regression ===
# =================================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage
conda activate envs/ldsc

# ld score regression
sumstatsCols="CHR,BP,ID,A1,A2,BETA,SE,P,N"
sumstatsColNames="CHR,POS,SNP,A1,A2,BETA,SE,P,N"
LDsnplist="data/ldsc/resources/w_hm3.noMHC.snplist"
LDchr="data/ldsc/resources/eur_w_ld_chr/"
LDbaselineH2="data/ldsc/resources/1000G_Phase3_baselineLD_v2.2/baselineLD."
LDbaselineTau="data/ldsc/resources/1000G_Phase3_baseline_v1.2/baseline."
LDweights="data/ldsc/resources/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
LDfrqfiles="data/ldsc/resources/1000G_Phase3_frq/1000G.EUR.QC."
partitioned=""
celltypes=""
celltypegroups=""
ctglabels=""
(for trait in gap_gm gap_wm gap_gwm; do (
    targetDir="results/${trait}/discovery/ldsc"
    sumstats="results/${trait}/discovery/gwas/sumstats.txt.gz"
	./code/genetics/ldsc.sh "${trait}" "${targetDir}" "${sumstats}" "${sumstatsCols}" "${sumstatsColNames}" "${LDsnplist}" "${LDchr}" "${LDbaselineH2}" "${LDbaselineTau}" "${LDweights}" "${LDfrqfiles}" "${partitioned}" "${celltypes}" "${celltypegroups}" "${ctglabels}"
	) &
done
wait)

# combine ldsc estimates in one file
conda activate envs/default
ldscFiles="results/gap_gm/discovery/ldsc/ldsc.h2.results results/gap_wm/discovery/ldsc/ldsc.h2.results results/gap_gwm/discovery/ldsc/ldsc.h2.results"
outFile="results/combined/discovery.ldsc.h2.txt"
./code/genetics/ldsc.combine.sh "${ldscFiles}" "${outFile}"

# ===============================================
# === discovery sample: create regional plots ===
# ===============================================

# set working directory and load conda environment
cd /slow/projectss/ukb_brainage
conda activate envs/locuszoom

# general settings
chrFilehandler='data/genetics/chr${i}/imp_mri_qc/vcf/chr${i}_mri_qc.vcf.gz'
pthresh=5E-8
flank="500kb"

# gap_gm
trait="gap_gm"
traitDescription="grey matter"
targetDir="results/${trait}/discovery/locuszoom/"
sumstats="results/${trait}/discovery/gwas/sumstats.txt.gz"
id="ID"
p="P"
conditionalFile="results/${trait}/discovery/conditional/conditional.cleaned.tophits.annovar.txt"

	# optional settings
	snpsManual="rs796226228 rs376899426 rs377158217 rs76122535 rs56072903"
	flankManual=""
	chrManual="1 2 2 2 17"
	startManual="214800000 197500000 200200000 203000000 43000000"
	endManual="216200000 199500000 201600000 205000000 45500000"
	./code/genetics/locuszoom.sh "${trait}" "${traitDescription}" "${targetDir}" "${chrFilehandler}" "${sumstats}" "${id}" "${p}" "${conditionalFile}" "${pthresh}" "${flank}" "${snpsManual}" "${flankManual}" "${chrManual}" "${startManual}" "${endManual}"

# gap_wm
trait="gap_wm"
traitDescription="white matter"
targetDir="results/${trait}/discovery/locuszoom/"
sumstats="results/${trait}/discovery/gwas/sumstats.txt.gz"
id="ID"
p="P"
conditionalFile="results/${trait}/discovery/conditional/conditional.cleaned.tophits.annovar.txt"

	# optional settings
	snpsManual="rs10932008 rs12287076 rs111854640"
	flankManual=""
	chrManual="2 11 17"
	startManual="203000000 46800000 43000000"
	endManual="205000000 48600000 45500000"
	./code/genetics/locuszoom.sh "${trait}" "${traitDescription}" "${targetDir}" "${chrFilehandler}" "${sumstats}" "${id}" "${p}" "${conditionalFile}" "${pthresh}" "${flank}" "${snpsManual}" "${flankManual}" "${chrManual}" "${startManual}" "${endManual}" &

# gap_gwm
trait="gap_gwm"
traitDescription="grey and white matter"
targetDir="results/${trait}/discovery/locuszoom/"
sumstats="results/${trait}/discovery/gwas/sumstats.txt.gz"
id="ID"
p="P"
conditionalFile="results/${trait}/discovery/conditional/conditional.cleaned.tophits.annovar.txt"

	# optional settings
	snpsManual="s550778111 rs377158217 rs76122535 rs199737620"
	flankManual=""
	chrManual="1 2 2 17"
	startManual="43300000 200200000 203000000 43000000"
	endManual="445000000 201600000 205000000 45500000"
	./code/genetics/locuszoom.sh "${trait}" "${traitDescription}" "${targetDir}" "${chrFilehandler}" "${sumstats}" "${id}" "${p}" "${conditionalFile}" "${pthresh}" "${flank}" "${snpsManual}" "${flankManual}" "${chrManual}" "${startManual}" "${endManual}"

# =======================================================================================
# === discovery sample: run cross-trait clumping to identify independent associations ===
# =======================================================================================

# set working directory
cd /slow/projects/ukb_brainage/
conda activate envs/default

# settings
traitlist="gap_gm gap_wm gap_gwm"
conditionalFiles="results/gap_gm/discovery/conditional/conditional.cleaned.tophits.annovar.txt results/gap_wm/discovery/conditional/conditional.cleaned.tophits.annovar.txt results/gap_gwm/discovery/conditional/conditional.cleaned.tophits.annovar.txt"
nonsynFiles="" # leave empty - use for GWAS meta-analysis results across discovery and replication
catalogFiles="" # leave empty - use for GWAS meta-analysis results across discovery and replication
finemapFiles="" # leave empty - use for GWAS meta-analysis results across discovery and replication
smrFiles="" # leave empty - use for GWAS meta-analysis results across discovery and replication
eqtlFiles="" # leave empty - use for GWAS meta-analysis results across discovery and replication
popsFiles="" # leave empty - use for GWAS meta-analysis results across discovery and replication
targetDir="results/combined/discovery.snplevel"
chrFilehandler='data/genetics/chr${i}/imp_mri_qc_EURjoined/bed/chr${i}_mri_qc'
LDsample="data/gap_gm/gap_gm.txt"
pthresh=1E-6
./code/genetics/snplevel.sh "${traitlist}" "${conditionalFiles}" "${nonsynFiles}" "${catalogFiles}" "${finemapFiles}" "${smrFiles}" "${eqtlFiles}" "${popsFiles}" "${targetDir}" "${chrFilehandler}" "${LDsample}" "${pthresh}"

# =======================================================
# === discovery sample: create manhattan and qq plots ===
# =======================================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage/
conda activate envs/default

# create manhattan plots
nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5
annotationThresh=5E-100 # pthresh=5E-8
sig=5E-8 # sig=5E-8
yend=35 # upper y axis limit
ysteps=5 # y axis breaks
width=11 # plot width (10.67 inch)
height=4 # plot height (4 inch)
preview=FALSE # only 1% of all variations is plotted
(for trait in gap_gm gap_wm gap_gwm; do (
	targetDir="results/${trait}/discovery/manhattan/"
	sumstats="results/${trait}/discovery/gwas/sumstats.txt.gz"
	annotationFile="results/${trait}/discovery/conditional/conditional.cleaned.tophits.annovar.txt"
	Rscript code/genetics/manhattan.R "${trait}" "${targetDir}" "${sumstats}" "${nCriterion}" "${annotationFile}" "${annotationThresh}" "${sig}" "${yend}" "${ysteps}" "${width}" "${height}" "${preview}"
	) &
done
wait)

# create qq-plots
pCol="P" # column containing p-values
nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5
prune=FALSE # ommit overplotting of variants with pval > 0.01
drawCI=FALSE # draw confidence interval?
drawLambda=FALSE # draw LambdaGC?
xend=8 # upper x axis limit
xsteps=2 # x axis breaks
yend=35 # upper y axis limit
ysteps=5 # y axis breaks
width=3 # plot width (4 inch)
height=4 # plot height (4 inch)
(for trait in gap_gm gap_wm gap_gwm; do (
	targetDir="results/${trait}/discovery/qqplot/"
	sumstats="results/${trait}/discovery/gwas/sumstats.txt.gz"
	Rscript code/genetics/qqplot.R "${trait}" "${targetDir}" "${sumstats}" "${pCol}" "${nCriterion}" "${prune}" "${drawCI}" "${drawLambda}" "${xend}" "${xsteps}" "${yend}" "${ysteps}" "${width}" "${height}"
	) &
done
wait)

# combine Manhattan and qq-plots
traits="gap_gm,gap_wm,gap_gwm"
plotTitles="" # annotation will be a-f if no plot titles are provided
manhattanPlots="results/gap_gm/discovery/manhattan/manhattan.png,results/gap_wm/discovery/manhattan/manhattan.png,results/gap_gwm/discovery/manhattan/manhattan.png"
qqPlots="results/gap_gm/discovery/qqplot/qqplot.png,results/gap_wm/discovery/qqplot/qqplot.png,results/gap_gwm/discovery/qqplot/qqplot.png"
outputFile="results/combined/discovery.manhattan.qq.png"
width=14
height=12
Rscript code/genetics/manhattan.qq.combine.R "${traits}" "${plotTitles}" "${manhattanPlots}" "${qqPlots}" "${outputFile}" "${width}" "${height}"

# ==============================================
# === replication sample: run GWAS and GWAMA ===
# ==============================================

# set working directory
cd "/slow/projects/ukb_brainage"

# calculate expected number of successfull replications (accounting for winner's curse bias)
conda activate envs/power
discoveryFile="results/combined/discovery.snplevel/snplevel.gws.txt"
traits="gap_gm,gap_wm,gap_gwm"
traitNames="grey matter,white matter,grey and white matter"
colTrait="TRAIT"
colBETA="BETA"
colSE="SE"
colN="N"
colFREQ="A1_FREQ"
cols2keep="LOCUS_COUNT,TRAIT,CHR,BP,ID,NEAREST_GENE,A1,A2,A1_FREQ,BETA,SE,Z,P,N"
filterCol="DISCOV_COUNT"
filterVal=1
replicationN=22256
replicationP=0.1 # set 0.1 for one-tailed p = 0.05
outFile="results/combined/replication.power"
Rscript ./code/genetics/replicate.power.R "${discoveryFile}" "${traits}" "${traitNames}" "${colTrait}" "${colBETA}" "${colSE}" "${colN}" "${colFREQ}" "${cols2keep}" "${filterCol}" "${filterVal}" "${replicationN}" "${replicationP}" "${outFile}"

# run GWAS in UKB replication samples
conda activate envs/default
chrFileHandle="data/genetics/chr\${i}/imp_mri_qc_\${ancestry}/chr\${i}_mri_qc"
data="data/traits/replicate.txt"
ancestriesCol="pan"
snpQC="data/genetics/meta/ukb_snp_qc.txt"
impMFI="data/genetics/ukb_imp_mfi"
maf=0.01
threads=100
ancestries="AFR,AMR,CSA,EAS,MID"
covs="sex,age,age2,ac1,ac2,ac3,TIV,array,PanC{1..4}"
for trait in gap_gm gap_wm gap_gwm; do
	targetDir="results/${trait}/replicate"
	taskset -c 0-99 code/genetics/replicate.sh "${trait}" "${targetDir}" "${chrFileHandle}" "${data}" "${covs}" "${ancestries}" "${ancestriesCol}" "${snpQC}" "${impMFI}" "${maf}" "${threads}" 
done
ancestries="EUR"
covs="sex,age,age2,ac1,ac2,ac3,TIV,array,PanC{1..20}"
for trait in gap_gm gap_wm gap_gwm; do
	targetDir="results/${trait}/replicate"
	taskset -c 0-99 code/genetics/replicate.sh "${trait}" "${targetDir}" "${chrFileHandle}" "${data}" "${covs}" "${ancestries}" "${ancestriesCol}" "${snpQC}" "${impMFI}" "${maf}" "${threads}" 
done

# add LIFE sumstats (harmonize)
ukbSNPs="data/genetics/qc_snplist/ukb_info_08_nodups.txt"
(for trait in gap_gm gap_wm gap_gwm; do (
	sumstatsLIFE="data/life/gwas/s303_1_GWAS_brain_noCovStats.brainage_${trait}_stack.glm.linear"
	targetDir="results/${trait}/replicate/LIFE/"
	code/genetics/harmonizeLIFE.sh "${sumstatsLIFE}" "${ukbSNPs}" "${targetDir}"
	) & 
done)
wait

# create manhattan plots
nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5 
conditionalFile="" # no conditional analysis
annotationThresh=5E-8 # pthresh=5E-8
sig=5E-8 # sig=5E-8
yend=12 # upper y axis limit
ysteps=4 # y axis breaks
width=11 # plot width (10.67 inch)
height=2.5 # plot height (4 inch)
preview=FALSE # only 1% of all variations is plotted
(for ancestry in AFR AMR CSA EAS EUR MID LIFE; do (
	for trait in gap_gm gap_wm gap_gwm; do (
		targetDir="results/${trait}/replicate/${ancestry}/"
		sumstats="results/${trait}/replicate/${ancestry}/sumstats.txt.gz"
		Rscript code/genetics/manhattan.R "${trait}" "${targetDir}" "${sumstats}" "${nCriterion}" "${conditionalFile}" "${annotationThresh}" "${sig}" "${yend}" "${ysteps}" "${width}" "${height}" "${preview}"
		) &
	done ) &
done
wait)

# create qq-plots
pCol="P" # column containing p-values
nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5 
prune=FALSE # ommit overplotting of variants with pval > 0.01
drawCI=FALSE # draw confidence interval?
drawLambda=FALSE # draw LambdaGC?
xend=8 # upper x axis limit
xsteps=2 # x axis breaks
yend=12 # upper y axis limit
ysteps=4 # y axis breaks
width=3 # plot width (3 inch)
height=2.5 # plot height (4 inch)
(for ancestry in AFR AMR CSA EAS EUR MID LIFE; do (
	for trait in gap_gm gap_wm gap_gwm; do (
		targetDir="results/${trait}/replicate/${ancestry}/"
		sumstats="results/${trait}/replicate/${ancestry}/sumstats.txt.gz"
		Rscript code/genetics/qqplot.R "${trait}" "${targetDir}" "${sumstats}" "${pCol}" "${nCriterion}" "${prune}" "${drawCI}" "${drawLambda}" "${xend}" "${xsteps}" "${yend}" "${ysteps}" "${width}" "${height}"
		) &
	done ) &
done
wait)

# combine manhattan and qq plots
sampleList="AFR,AMR,CSA,EAS,EUR,MID,LIFE"
plotTitles="UKB_African,UKB_American,UKB_Central/South_Asian,UKB_East_Asian,UKB_European,UKB_Middle_Eastern,LIFE_European" # underscore "_" will be replaced by white space
width=14
height=17.5
(for trait in gap_gm gap_wm gap_gwm; do (
	manhattanPlots=$(echo $(for sample in $(echo $sampleList | sed 's/,/ /g'); do echo results/${trait}/replicate/${sample}/manhattan.png; done) | sed 's/ /,/g' )
	qqPlots=$(echo $(for sample in $(echo $sampleList | sed 's/,/ /g'); do echo results/${trait}/replicate/${sample}/qqplot.png; done) | sed 's/ /,/g' )
	outputFile="results/${trait}/replicate/replicate.qqplot.manhattan.png"
	Rscript code/genetics/manhattan.qq.combine.R "${sampleList}" "${plotTitles}" "${manhattanPlots}" "${qqPlots}" "${outputFile}" "${width}" "${height}"
	) &
done
wait)

# run METAL across EUR ancestry replication cohorts
cohorts="EUR,LIFE" # set output file handler & make sure that the order matches the list of sumstats (below)
type="ivweight" # type="nweight"
idCol="ID"
carryOver="CHR,BP"
leaveOneOut=0
(for trait in gap_gm gap_wm gap_gwm; do (
	targetDir="results/${trait}/replicate/metal.eur/"; mkdir -p ${targetDir}
	sumstats=$(echo $(for i in $(echo $cohorts | sed 's/,/ /g'); do echo "results/${trait}/replicate/${i}/sumstats.txt.gz"; done) | sed 's/ /,/g')
	code/genetics/metal.sh ${trait} ${targetDir} ${cohorts} ${sumstats} ${type} ${idCol} ${carryOver} ${leaveOneOut}
	) &
done
wait)

	# quality-check METAL results
	nCriterion=TRUE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5
	hetP=1E-6 # heterogeneity p threshold for excluding SNPs
	(for trait in gap_gm gap_wm gap_gwm; do (
		sumstats="results/${trait}/replicate/metal.eur/metal.ivweight.gz"
		outFile="results/${trait}/replicate/metal.eur/metal.ivweight.qc"
		Rscript code/genetics/metal.qc.R ${sumstats} ${outFile} ${nCriterion} ${hetP}
		) &
	done
	wait)

	# run conditional analysis and do positional annotation (!! do not run in parallel across traits)
	subsFile="data/genetics/chr1/imp_mri_qc_EUR/chr1_mri_qc.psam"
	chrFileHandle="data/genetics/chr\$i/imp_mri_qc_EUR/bed/chr\${i}_mri_qc"
	cols="ID,A1,A2,A1_FREQ,BETA,SE,P,N,CHR,BP,A1_FREQ_SE,A1_FREQ_MIN,A1_FREQ_MAX,DIRECTION,HET_ISQ,HetChiSq,HET_DF,HET_P" # first 10 columns have the following order ID,A1,A2,A1_FREQ,BETA,SE,P,N,CHR,BP
	pthresh=1E-6
	threads=100
	humandb="data/annovar/humandb"
	refseq="data/refseq/GRCh37_latest_genomic.edit.gff.gz"
	for trait in gap_gm gap_wm gap_gwm; do
		targetDir="results/${trait}/replicate/metal.eur/conditional"
		sumstats="results/${trait}/replicate/metal.eur/metal.ivweight.qc.gz"
		sumstatsPLINK="results/${trait}/replicate/EUR/sumstats_plink.txt.gz"
		./code/genetics/conditional.sh "${subsFile}" "${targetDir}" "${chrFileHandle}" "${sumstats}" "${cols}" "${sumstatsPLINK}" "${pthresh}" "${threads}" "${humandb}" "${refseq}"
	done

	# create GWAMA manhattan plots
	nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5 
	annotationThresh=5E-100 # pthresh=5E-8
	sig=5E-8 # sig=5E-8
	yend=20 # upper y axis limit
	ysteps=5 # y axis breaks
	width=11 # plot width (10.67 inch)
	height=2.5 # plot height (4 inch)
	preview=FALSE # only 1% of all variations is plotted
	(for trait in gap_gm gap_wm gap_gwm; do (
		targetDir="results/${trait}/replicate/metal.eur/"
		sumstats="results/${trait}/replicate/metal.eur/metal.ivweight.qc.gz"
		annotationFile="results/${trait}/replicate/metal.eur/conditional/conditional.cleaned.tophits.annovar.txt"
		Rscript code/genetics/manhattan.R "${trait}" "${targetDir}" "${sumstats}" "${nCriterion}" "${annotationFile}" "${annotationThresh}" "${sig}" "${yend}" "${ysteps}" "${width}" "${height}" "${preview}"
		) &
	done
	wait)

	# create GWAMA qq-plots
	pCol="P" # column containing p-values
	nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5 
	prune=FALSE # ommit overplotting of variants with pval > 0.01
	drawCI=FALSE # draw confidence interval?
	drawLambda=FALSE # draw LambdaGC?
	xend=8 # upper x axis limit
	xsteps=2 # x axis breaks
	yend=20 # upper y axis limit
	ysteps=5 # y axis breaks
	width=3 # plot width (3 inch)
	height=2.5 # plot height (4 inch)
	(for trait in gap_gm gap_wm gap_gwm; do (
		targetDir="results/${trait}/replicate/metal.eur/"
		sumstats="results/${trait}/replicate/metal.eur/metal.ivweight.qc.gz"
		Rscript code/genetics/qqplot.R "${trait}" "${targetDir}" "${sumstats}" "${pCol}" "${nCriterion}" "${prune}" "${drawCI}" "${drawLambda}" "${xend}" "${xsteps}" "${yend}" "${ysteps}" "${width}" "${height}"
		) &
	done
	wait)

	# combine GWAMA manhattan and qq plots
	traitList="gap_gm,gap_wm,gap_gwm"
	plotTitles="Grey_matter,White_matter,Grey_and_white_matter" # underscore "_" will be replaced by white space
	width=14
	height=8
	manhattanPlots=$(echo $(for trait in $(echo $traitList | sed 's/,/ /g'); do echo results/${trait}/replicate/metal.eur/manhattan.png; done) | sed 's/ /,/g' )
	qqPlots=$(echo $(for trait in $(echo $traitList | sed 's/,/ /g'); do echo results/${trait}/replicate/metal.eur/qqplot.png; done) | sed 's/ /,/g' )
	outputFile="results/combined/replicate.eur.qqplot.manhattan.png"
	Rscript code/genetics/manhattan.qq.combine.R "${traitList}" "${plotTitles}" "${manhattanPlots}" "${qqPlots}" "${outputFile}" "${width}" "${height}"

	# get replication results of discovery index variants
	replicationNcovs=28 # sex,age,age2,TIV,PC1-4 (28 in case of ac1-3,PC1-20,array)
	(for trait in gap_gm gap_wm gap_gwm; do (
		discoveryFile="results/${trait}/discovery/conditional/conditional.cleaned.tophits.annovar.txt"
		replicationFile="results/${trait}/replicate/metal.eur/metal.ivweight.qc.gz"
		outFile="results/${trait}/replicate/metal.eur/replicateDiscoveries"
		Rscript code/genetics/replicateDiscoveries.R "${discoveryFile}" "${replicationFile}" "${replicationNcovs}" "${outFile}"
		) &
	done
	wait)
	
	# combine replication results of discovery index variants (merge loci across traits)
	snplevelFile="results/combined/discovery.snplevel/snplevel.txt"
	replicateFiles="results/gap_gm/replicate/metal.eur/replicateDiscoveries.txt,results/gap_wm/replicate/metal.eur/replicateDiscoveries.txt,results/gap_gwm/replicate/metal.eur/replicateDiscoveries.txt"
	outFile="results/combined/replicateDiscoveries.eur"
	traits="gap_gm,gap_wm,gap_gwm"
	traitNames="grey matter,white matter,grey and white matter"
	Rscript code/genetics/replicateDiscoveries.combine.R "${snplevelFile}" "${replicateFiles}" "${outFile}" "${traits}" "${traitNames}"

		# draw qq plots
		inputFile="results/combined/replicateDiscoveries.eur.txt"
		pCol="REPLIC_P"
		criterionCols="DISCOV_TRAITNAME"
		duplicatedCol=""
		prune=FALSE
		drawCI=TRUE
		xend=2
		xsteps=1
		yend=35
		ysteps=5
		width=2.5
		height=4
		outFile=("results/gap_gm/replicate/metal.eur/replicateDiscoveries.eur.gm.png" "results/gap_wm/replicate/metal.eur/replicateDiscoveries.eur.wm.png" "results/gap_gwm/replicate/metal.eur/replicateDiscoveries.eur.gwm.png")
		criterionVals=("grey matter" "white matter" "grey and white matter")
		for (( i=0; i<${#outFile[@]}; i++ )); do
			Rscript code/genetics/replicateDiscoveries.qq.R "${inputFile}" "${outFile[i]}" "${pCol}" "${criterionCols}" "${criterionVals[i]}" "${duplicatedCol}" "${prune}" "${drawCI}" "${xend}" "${xsteps}" "${yend}" "${ysteps}" "${width}" "${height}"
		done
		xend=2.2
		outFile="results/combined/replicateDiscoveries.eur.crosstrait.png"
		criterionCols=""
		criterionVals=""
		duplicatedCol="DISCOV_LOCUS_COUNT"
		Rscript code/genetics/replicateDiscoveries.qq.R "${inputFile}" "${outFile}" "${pCol}" "${criterionCols}" "${criterionVals}" "${duplicatedCol}" "${prune}" "${drawCI}" "${xend}" "${xsteps}" "${yend}" "${ysteps}" "${width}" "${height}"

		# combine plots
		plotTitles="Grey_matter,White_matter,Grey_and_white_matter"
		plots="results/gap_gm/replicate/metal.eur/replicateDiscoveries.eur.gm.png,results/gap_wm/replicate/metal.eur/replicateDiscoveries.eur.wm.png,results/gap_gwm/replicate/metal.eur/replicateDiscoveries.eur.gwm.png"
		outFile="results/combined/replicateDiscoveries.eur.qq.png"
		width=7.5
		height=4
		Rscript code/genetics/replicateDiscoveries.qq.combine.R "${plotTitles}" "${plots}" "${outFile}" "${width}" "${height}"

	# run ld score regression across eur ancestry replication results
	conda activate envs/ldsc
	LDsnplist="data/ldsc/resources/w_hm3.noMHC.snplist"
	LDchr="data/ldsc/resources/eur_w_ld_chr/"
	LDbaseline="data/ldsc/resources/baselineLD/baselineLD."
	LDweights="data/ldsc/resources/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
	LDfrqfiles="data/ldsc/resources/1000G_Phase3_frq/1000G.EUR.QC."
	partitioned=0
	sumstatsCols="CHR,BP,ID,A1,A2,BETA,SE,P,N"
	sumstatsColNames="CHR,POS,SNP,A1,A2,BETA,SE,P,N"
	for trait in gap_gm gap_wm gap_gwm; do
	    targetDir="results/${trait}/replicate/metal.eur/ldsc"
	    sumstats="results/${trait}/replicate/metal.eur/metal.ivweight.qc.gz"
		./code/genetics/ldsc.sh "${trait}" "${targetDir}" "${sumstats}" "${sumstatsCols}" "${sumstatsColNames}" "${LDsnplist}" "${LDchr}" "${LDbaseline}" "${LDweights}" "${LDfrqfiles}" "${partitioned}"
	done

		# combine estimates in one file
		conda activate envs/default
		traitList="gap_gm gap_wm gap_gwm"
		ldscFiles="results/gap_gm/replicate/metal.eur/ldsc/ldsc.h2.results results/gap_wm/replicate/metal.eur/ldsc/ldsc.h2.results results/gap_gwm/replicate/metal.eur/ldsc/ldsc.h2.results"
		outFile="results/combined/replicate.eur.ldsc.h2.txt"
		./code/genetics/ldsc.combine.sh "${ldscFiles}" "${outFile}"

		# calculate rg between discovery and replication results
		conda activate envs/ldsc
		LDchr="data/ldsc/resources/eur_w_ld_chr/"
		xprefix="discovery."
		yprefix="replication."
		threads=100 # set number of parallel workers
		for trait in gap_gm gap_wm gap_gwm; do
		    sumstatsX="results/${trait}/discovery/ldsc/ldsc.${trait}.sumstats.gz"
		    sumstatsY="results/${trait}/replicate/metal.eur/ldsc/ldsc.${trait}.sumstats.gz"
		    targetDir="results/${trait}/replicate/metal.eur/rgdiscov"
			./code/genetics/rg.sh "${targetDir}" "${sumstatsX}" "${sumstatsY}" "${xprefix}" "${yprefix}" "${LDchr}" "${threads}"
		done

		# combine estimates of rg between discovery and replication
		awk 'NR==1 { print } FNR==1 { next } { print }' "results/gap_gm/replicate/metal.eur/rgdiscov/rg.results" "results/gap_wm/replicate/metal.eur/rgdiscov/rg.results" "results/gap_gwm/replicate/metal.eur/rgdiscov/rg.results" \
			> results/combined/replicate.eur.rgdiscov

		# create supplementum table with h2 and rg
		header=trait$'\t'$'\t'discovery.snp_count$'\t'discovery.h2$'\t'discovery.h2_se$'\t'discovery.lambda_GC$'\t'discovery.chi2$'\t'discovery.intercept$'\t'discovery.intercept_se$'\t'discovery.ratio$'\t'discovery.ratio_se$'\t'$'\t'replication.snp_count$'\t'replication.h2$'\t'replication.h2_se$'\t'replication.lambda_GC$'\t'replication.chi2$'\t'replication.intercept$'\t'replication.intercept_se$'\t'replication.ratio$'\t'replication.ratio_se$'\t'$'\t'rg$'\t'se$'\t'z$'\t'p$'\t'gcov_int$'\t'gcov_int_se
		awk -v header="${header}" 'BEGIN { print header } 
			FNR==1 { filenum++; next}
			filenum==1 { trait[FNR-1]=$1"\t\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 }
			filenum==2 { trait[FNR-1]=trait[FNR-1]"\t\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 }
			filenum==3 { trait[FNR-1]=trait[FNR-1]"\t\t"$3"\t"$4"\t"$5"\t"$6"\t"$11"\t"$12; print trait[FNR-1] }' results/combined/discovery.ldsc.h2.txt results/combined/replicate.eur.ldsc.h2.txt results/combined/replicate.eur.rgdiscov \
			> results/combined/replicate.eur.rgdiscov.suppl.txt

# run MR-MEGA across multi-ancestry replication cohorts (also run random and fixed-effects GWAMA and add weights to MR-MEGA results)
(for trait in gap_gm gap_wm gap_gwm; do (
	inputFiles=$(echo results/${trait}/replicate/{AFR,AMR,CSA,EAS,EUR,LIFE,MID}/sumstats.txt.gz | sed 's/ /,/g')
	targetDir="results/${trait}/replicate/mrmega.all/"
	cols="CHR,BP,ID,A1,A2,A1_FREQ,BETA,SE,N"
	./code/genetics/mrmega.sh "${inputFiles}" "${targetDir}" "${cols}"
	) &
done
wait)

	# run positional clumping
	idCol="ID"
	chrCol="CHR"
	bpCol="BP"
	pCol="P"
	maxP=5e-8
	decrease=FALSE
	windowSize=3000000
	printLeadOnly=TRUE
	(for trait in gap_gm gap_wm gap_gwm; do (
		inFile="results/${trait}/replicate/mrmega.all/mrmega.weights.gz"
		outFile="results/${trait}/replicate/mrmega.all/mrmega.weights.clumped.txt"
		code/genetics/positional.clumping.R "${inFile}" "${outFile}" "${idCol}" "${chrCol}" "${bpCol}" "${pCol}" "${maxP}" "${decrease}" "${windowSize}" "${printLeadOnly}"
		) &
	done)

	# create MR-MEGA manhattan plots	
	nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5 
	annotationThresh=5E-100 # pthresh=5E-8
	sig=5E-8 # sig=5E-8
	yend=20 # upper y axis limit
	ysteps=5 # y axis breaks
	width=11 # plot width (10.67 inch)
	height=2.5 # plot height (4 inch)
	preview=FALSE # only 1% of all variations is plotted
	(for trait in gap_gm gap_wm gap_gwm; do (
		targetDir="results/${trait}/replicate/mrmega.all/"
		sumstats="results/${trait}/replicate/mrmega.all/mrmega.weights.gz"
		awk -F'\t' 'NR==1 { print "ID", "SNP_CNT", "P", "LEAD_SNP_pJ", "NEAREST_GENE" }
			NR>1 { print $1, $(NF-1), $20, $20, "BLANK" }' OFS='\t' "results/${trait}/replicate/mrmega.all/mrmega.weights.clumped.txt" > "${targetDir}/manhattan.annotation.txt"
		annotationFile="${targetDir}/manhattan.annotation.txt"
		Rscript code/genetics/manhattan.R "${trait}" "${targetDir}" "${sumstats}" "${nCriterion}" "${annotationFile}" "${annotationThresh}" "${sig}" "${yend}" "${ysteps}" "${width}" "${height}" "${preview}"
		rm -f "${annotationFile}"
		) &
	done
	wait)

	# create MR-MEGA qq-plots
	pCol="P" # column containing p-values
	nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5 
	prune=FALSE # ommit overplotting of variants with pval > 0.01
	drawCI=FALSE # draw confidence interval?
	drawLambda=FALSE # draw LambdaGC?
	xend=8 # upper x axis limit
	xsteps=2 # x axis breaks
	yend=20 # upper y axis limit
	ysteps=5 # y axis breaks
	width=3 # plot width (3 inch)
	height=2.5 # plot height (4 inch)
	(for trait in gap_gm gap_wm gap_gwm; do (
		targetDir="results/${trait}/replicate/mrmega.all/"
		sumstats="results/${trait}/replicate/mrmega.all/mrmega.weights.gz"
		Rscript code/genetics/qqplot.R "${trait}" "${targetDir}" "${sumstats}" "${pCol}" "${nCriterion}" "${prune}" "${drawCI}" "${drawLambda}" "${xend}" "${xsteps}" "${yend}" "${ysteps}" "${width}" "${height}"
		) &
	done
	wait)

	# combine MR-MEGA manhattan and qq plots
	traitList="gap_gm,gap_wm,gap_gwm"
	plotTitles="Grey_matter,White_matter,Grey_and_white_matter" # underscore "_" will be replaced by white space
	width=14
	height=8
	manhattanPlots=$(echo $(for trait in $(echo $traitList | sed 's/,/ /g'); do echo results/${trait}/replicate/mrmega.all/manhattan.png; done) | sed 's/ /,/g' )
	qqPlots=$(echo $(for trait in $(echo $traitList | sed 's/,/ /g'); do echo results/${trait}/replicate/mrmega.all/qqplot.png; done) | sed 's/ /,/g' )
	outputFile="results/combined/replicate.all.qqplot.manhattan.png"
	Rscript code/genetics/manhattan.qq.combine.R "${traitList}" "${plotTitles}" "${manhattanPlots}" "${qqPlots}" "${outputFile}" "${width}" "${height}"

	# get MR-MEGA replication results of discovery index variants
	replicationNcovs=28 # sex,age,age2,TIV,PC1-4 (28 in case of ac1-3,PC1-20,array)
	(for trait in gap_gm gap_wm gap_gwm; do (
		discoveryFile="results/${trait}/discovery/conditional/conditional.cleaned.tophits.annovar.txt"
		replicationFile="results/${trait}/replicate/mrmega.all/mrmega.weights.gz"
		outFile="results/${trait}/replicate/mrmega.all/replicateDiscoveries"
		Rscript code/genetics/replicateDiscoveries.mrmega.R "${discoveryFile}" "${replicationFile}" "${replicationNcovs}" "${outFile}"
		) &
	done
	wait)
	
	# combine MR-MEGA replication results of discovery index variants (merge loci across traits)
	snplevelFile="results/combined/discovery.snplevel/snplevel.txt"
	replicateFiles="results/gap_gm/replicate/mrmega.all/replicateDiscoveries.txt,results/gap_wm/replicate/mrmega.all/replicateDiscoveries.txt,results/gap_gwm/replicate/mrmega.all/replicateDiscoveries.txt"
	outFile="results/combined/replicateDiscoveries.all"
	traits="gap_gm,gap_wm,gap_gwm"
	traitNames="grey matter,white matter,grey and white matter"
	Rscript code/genetics/replicateDiscoveries.mrmega.combine.R "${snplevelFile}" "${replicateFiles}" "${outFile}" "${traits}" "${traitNames}"

		# draw qq plots
		inputFile="results/combined/replicateDiscoveries.all.txt"
		pCol="REPLIC_P"
		criterionCols="DISCOV_TRAITNAME"
		duplicatedCol=""
		prune=FALSE
		drawCI=TRUE
		xend=2
		xsteps=1
		yend=35
		ysteps=5
		width=2.5
		height=4
		outFile=("results/gap_gm/replicate/mrmega.all/replicateDiscoveries.png" "results/gap_wm/replicate/mrmega.all/replicateDiscoveries.png" "results/gap_gwm/replicate/mrmega.all/replicateDiscoveries.png")
		criterionVals=("grey matter" "white matter" "grey and white matter")
		for (( i=0; i<${#outFile[@]}; i++ )); do
			Rscript code/genetics/replicateDiscoveries.qq.R "${inputFile}" "${outFile[i]}" "${pCol}" "${criterionCols}" "${criterionVals[i]}" "${duplicatedCol}" "${prune}" "${drawCI}" "${xend}" "${xsteps}" "${yend}" "${ysteps}" "${width}" "${height}"
		done
		xend=2.2
		outFile="results/combined/replicateDiscoveries.all.crosstrait.png"
		criterionCols=""
		criterionVals=""
		duplicatedCol="DISCOV_LOCUS_COUNT"
		Rscript code/genetics/replicateDiscoveries.qq.R "${inputFile}" "${outFile}" "${pCol}" "${criterionCols}" "${criterionVals}" "${duplicatedCol}" "${prune}" "${drawCI}" "${xend}" "${xsteps}" "${yend}" "${ysteps}" "${width}" "${height}"

		# combine plots
		plotTitles="Grey_matter,White_matter,Grey_and_white_matter"
		plots="results/gap_gm/replicate/mrmega.all/replicateDiscoveries.png,results/gap_wm/replicate/mrmega.all/replicateDiscoveries.png,results/gap_gwm/replicate/mrmega.all/replicateDiscoveries.png"
		outFile="results/combined/replicateDiscoveries.all.qq.png"
		width=7.5
		height=4
		Rscript code/genetics/replicateDiscoveries.qq.combine.R "${plotTitles}" "${plots}" "${outFile}" "${width}" "${height}"

# show number of variants
for ancestry in AFR AMR CSA EAS EUR MID LIFE; do
	echo ${ancestry} $(zcat results/gap_gm/replicate/${ancestry}/sumstats.txt.gz | wc -l)
done

# show sample size
for ancestry in AFR AMR CSA EAS EUR MID LIFE; do
	echo ${ancestry} $(zcat results/gap_gm/replicate/${ancestry}/sumstats.txt.gz | head -2 | awk -F'\t' 'NR==2 { print $11 }')
done

# ==================================================
# === Discovery + Replication: run meta-analysis ===
# ==================================================

# set working directory
cd "/slow/projects/ukb_brainage"
conda activate envs/default

# run METAL across EUR ancestry discovery and replication
cohorts="DISCOV,EUR,LIFE" # set output file handler & make sure that the order matches the list of sumstats (below)
type="ivweight" # type="nweight"
idCol="ID"
carryOver="CHR,BP"
leaveOneOut=0
(for trait in gap_gm gap_wm gap_gwm; do (
	targetDir="results/${trait}/gwama/eur/"; mkdir -p ${targetDir}
	sumstats="results/${trait}/gwas/sumstats.txt.gz,results/${trait}/replicate/EUR/sumstats.txt.gz,results/${trait}/replicate/LIFE/sumstats.txt.gz"
	code/genetics/metal.sh ${trait} ${targetDir} ${cohorts} ${sumstats} ${type} ${idCol} ${carryOver} ${leaveOneOut}
	) &
done
wait)

	# quality-check METAL results
	nCriterion=TRUE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5
	hetP=1E-6 # heterogeneity p threshold for excluding SNPs
	(for trait in gap_gm gap_wm gap_gwm; do (
		sumstats="results/${trait}/gwama/eur/metal.ivweight.gz"
		outFile="results/${trait}/gwama/eur/metal.ivweight.qc"
		Rscript code/genetics/metal.qc.R ${sumstats} ${outFile} ${nCriterion} ${hetP}
		) &
	done
	wait)

	# run conditional analysis and do positional annotation (!! do not run in parallel across traits)
	conda activate envs/default
	subsFile="data/genetics/chr1/imp_mri_qc_EURjoined/chr1_mri_qc.psam"
	chrFileHandle="data/genetics/chr\${i}/imp_mri_qc_EURjoined/bed/chr\${i}_mri_qc"
	cols="ID,A1,A2,A1_FREQ,BETA,SE,P,N,CHR,BP,A1_FREQ_SE,A1_FREQ_MIN,A1_FREQ_MAX,DIRECTION,HET_ISQ,HetChiSq,HET_DF,HET_P" # first 10 columns have the following order ID,A1,A2,A1_FREQ,BETA,SE,P,N,CHR,BP
	pthresh=1E-6
	threads=100
	humandb="data/annovar/humandb"
	refseq="data/refseq/GRCh37_latest_genomic.edit.gff.gz"
	for trait in gap_gm gap_wm gap_gwm; do
		targetDir="results/${trait}/gwama/eur/conditional"
		sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz"
		sumstatsPLINK="results/${trait}/gwas/sumstats_plink.txt.gz"
		./code/genetics/conditional.sh "${subsFile}" "${targetDir}" "${chrFileHandle}" "${sumstats}" "${cols}" "${sumstatsPLINK}" "${pthresh}" "${threads}" "${humandb}" "${refseq}"
	done

	# create GWAMA manhattan plots
	nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5 
	annotationThresh=1E-1-0 # pthresh=5E-8
	sig=5E-8 # sig=5E-8
	yend=40 # upper y axis limit
	ysteps=10 # y axis breaks
	width=11 # plot width (10.67 inch)
	height=4 # plot height (4 inch)
	preview=FALSE # only 1% of all variations is plotted
	(for trait in gap_gm gap_wm gap_gwm; do (
		targetDir="results/${trait}/gwama/eur/"
		sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz"
		annotationFile="results/${trait}/gwama/eur/conditional/conditional.cleaned.tophits.annovar.txt"
		Rscript code/genetics/manhattan.R "${trait}" "${targetDir}" "${sumstats}" "${nCriterion}" "${annotationFile}" "${annotationThresh}" "${sig}" "${yend}" "${ysteps}" "${width}" "${height}" "${preview}"
		) &
	done
	wait)

	# create GWAMA qq-plots
	pCol="P" # column containing p-values
	nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5 
	prune=FALSE # ommit overplotting of variants with pval > 0.01
	drawCI=FALSE # draw confidence interval?
	drawLambda=FALSE # draw LambdaGC?
	xend=8 # upper x axis limit
	xsteps=2 # x axis breaks
	yend=40 # upper y axis limit
	ysteps=10 # y axis breaks
	width=3 # plot width (3 inch)
	height=4 # plot height (4 inch)
	(for trait in gap_gm gap_wm gap_gwm; do (
		targetDir="results/${trait}/gwama/eur/"
		sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz"
		Rscript code/genetics/qqplot.R "${trait}" "${targetDir}" "${sumstats}" "${pCol}" "${nCriterion}" "${prune}" "${drawCI}" "${drawLambda}" "${xend}" "${xsteps}" "${yend}" "${ysteps}" "${width}" "${height}"
		) &
	done
	wait)

	# combine GWAMA manhattan and qq plots
	traitList="gap_gm,gap_wm,gap_gwm"
	plotTitles="" # underscore "_" will be replaced by white space; leave empty for a-f annotation
	width=14
	height=12
	manhattanPlots=$(echo $(for trait in $(echo $traitList | sed 's/,/ /g'); do echo results/${trait}/gwama/eur/manhattan.png; done) | sed 's/ /,/g' )
	qqPlots=$(echo $(for trait in $(echo $traitList | sed 's/,/ /g'); do echo results/${trait}/gwama/eur/qqplot.png; done) | sed 's/ /,/g' )
	outputFile="results/combined/gwama.eur.qqplot.manhattan.png"
	Rscript code/genetics/manhattan.qq.combine.R "${traitList}" "${plotTitles}" "${manhattanPlots}" "${qqPlots}" "${outputFile}" "${width}" "${height}"

# run MR-MEGA across multi-ancestry discovery and replication cohorts (also run random and fixed-effects GWAMA and add weights to MR-MEGA results)
(for trait in gap_gm gap_wm gap_gwm; do (
	mv "results/${trait}/gwama/all" "results/${trait}/gwama/all.backup"
	inputFiles=$(echo "results/${trait}/discovery/gwas/sumstats.txt.gz,"$(echo results/${trait}/replicate/{AFR,AMR,CSA,EAS,EUR,LIFE,MID}/sumstats.txt.gz) | sed 's/ /,/g')
	targetDir="results/${trait}/gwama/all/"
	cols="CHR,BP,ID,A1,A2,A1_FREQ,BETA,SE,N"
	./code/genetics/mrmega.sh "${inputFiles}" "${targetDir}" "${cols}"
	) &
done
wait)

	# run positional clumping
	idCol="ID"
	chrCol="CHR"
	bpCol="BP"
	pCol="P"
	maxP=5e-8
	decrease=FALSE
	windowSize=3000000
	printLeadOnly=TRUE
	(for trait in gap_gm gap_wm gap_gwm; do (
		inFile="results/${trait}/gwama/all/mrmega.weights.gz"
		outFile="results/${trait}/gwama/all/mrmega.weights.clumped.txt"
		code/genetics/positional.clumping.R "${inFile}" "${outFile}" "${idCol}" "${chrCol}" "${bpCol}" "${pCol}" "${maxP}" "${decrease}" "${windowSize}" "${printLeadOnly}"
		) &
	done
	wait)

	# create MR-MEGA manhattan plots	
	nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5 
	annotationThresh=5E-100 # pthresh=5E-8
	sig=5E-8 # sig=5E-8
	yend=40 # upper y axis limit
	ysteps=10 # y axis breaks
	width=11 # plot width (10.67 inch)
	height=4 # plot height (4 inch)
	preview=FALSE # only 1% of all variations is plotted
	(for trait in gap_gm gap_wm gap_gwm; do (
		targetDir="results/${trait}/gwama/all/"
		sumstats="results/${trait}/gwama/all/mrmega.weights.gz"
		awk -F'\t' 'NR==1 { print "ID", "SNP_CNT", "P", "LEAD_SNP_pJ", "NEAREST_GENE" }
			NR>1 { print $1, $(NF-1), $20, $20, "BLANK" }' OFS='\t' "results/${trait}/gwama/all/mrmega.weights.clumped.txt" > "${targetDir}/manhattan.annotation.txt"
		annotationFile="${targetDir}/manhattan.annotation.txt"
		Rscript code/genetics/manhattan.R "${trait}" "${targetDir}" "${sumstats}" "${nCriterion}" "${annotationFile}" "${annotationThresh}" "${sig}" "${yend}" "${ysteps}" "${width}" "${height}" "${preview}"
		rm -f "${annotationFile}"
		) &
	done) &

	# create MR-MEGA qq-plots
	pCol="P" # column containing p-values
	nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5 
	prune=FALSE # ommit overplotting of variants with pval > 0.01
	drawCI=FALSE # draw confidence interval?
	drawLambda=FALSE # draw LambdaGC?
	xend=8 # upper x axis limit
	xsteps=2 # x axis breaks
	yend=40 # upper y axis limit
	ysteps=10 # y axis breaks
	width=3 # plot width (3 inch)
	height=4 # plot height (4 inch)
	(for trait in gap_gm gap_wm gap_gwm; do (
		targetDir="results/${trait}/gwama/all/"
		sumstats="results/${trait}/gwama/all/mrmega.weights.gz"
		Rscript code/genetics/qqplot.R "${trait}" "${targetDir}" "${sumstats}" "${pCol}" "${nCriterion}" "${prune}" "${drawCI}" "${drawLambda}" "${xend}" "${xsteps}" "${yend}" "${ysteps}" "${width}" "${height}"
		) &
	done
	wait)

	# combine MR-MEGA manhattan and qq plots
	traitList="gap_gm,gap_wm,gap_gwm"
	plotTitles="Grey_matter,White_matter,Grey_and_white_matter" # underscore "_" will be replaced by white space
	width=14
	height=12
	manhattanPlots=$(echo $(for trait in $(echo $traitList | sed 's/,/ /g'); do echo results/${trait}/gwama/all/manhattan.png; done) | sed 's/ /,/g' )
	qqPlots=$(echo $(for trait in $(echo $traitList | sed 's/,/ /g'); do echo results/${trait}/gwama/all/qqplot.png; done) | sed 's/ /,/g' )
	outputFile="results/combined/gwama.all.qqplot.manhattan.png"
	Rscript code/genetics/manhattan.qq.combine.R "${traitList}" "${plotTitles}" "${manhattanPlots}" "${qqPlots}" "${outputFile}" "${width}" "${height}"

# show number of variants
for trait in gap_gm gap_wm gap_gwm; do
	echo ${trait} $(zcat results/${trait}/gwama/eur/metal.ivweight.qc.gz | wc -l)
done

for trait in gap_gm gap_wm gap_gwm; do
	echo ${trait} $(zcat results/${trait}/gwama/all/mrmega.weights.gz | wc -l)
done


# ======================================================
# === Discovery + replication: create regional plots ===
# ======================================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage
conda activate envs/locuszoom

# general settings
chrFilehandler='data/genetics/chr${i}/imp_mri_qc_EURjoined/vcf/chr${i}_mri_qc.vcf.gz'
pthresh=5E-8
flank="500kb"

# gap_gm
trait="gap_gm"
traitDescription="grey matter"
targetDir="results/${trait}/gwama/eur/locuszoom/"
sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz"
id="ID"
p="P"
conditionalFile="results/${trait}/gwama/eur/conditional/conditional.cleaned.tophits.annovar.txt"

	# optional settings
	snpsManual="rs1452628 rs12619333 rs1367858 rs530464314 rs2253491 6:126690257_AT_A rs12547781 rs76839250 17:27962571_GCCC_G rs56072903 X:107888149_CAA_C X:133781440_CTG_C"
	flankManual=""
	chrManual="1 2 2 2 6 6 8 16 17 17 X X"
	startManual="214800000 197500000 200200000 203000000 30500000 126200000 130400000 89200000 27400000 43000000 107000000 133000000"
	endManual="216200000 199500000 201600000 205000000 32000000 127400000 131600000 90400000 28800000 45500000 108800000 134200000"
	./code/genetics/locuszoom.sh "${trait}" "${traitDescription}" "${targetDir}" "${chrFilehandler}" "${sumstats}" "${id}" "${p}" "${conditionalFile}" "${pthresh}" "${flank}" "${snpsManual}" "${flankManual}" "${chrManual}" "${startManual}" "${endManual}" &

# gap_wm
trait="gap_wm"
traitDescription="white matter"
targetDir="results/${trait}/gwama/eur/locuszoom/"
sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz"
conditionalFile="results/${trait}/gwama/eur/conditional/conditional.cleaned.tophits.annovar.txt"

	# optional settings
	snpsManual="rs78300504 2:204127320_AAC_A rs2253491 rs10096381 rs755594165 rs111513543 rs2260227"
	flankManual=""
	chrManual="1 2 6 8 9 17 17"
	startManual="214600000 203000000 30500000 9500000 127400000 19400000 43000000"
	endManual="216000000 205000000 32000000 12000000 128800000 20800000 45500000"
	./code/genetics/locuszoom.sh "${trait}" "${traitDescription}" "${targetDir}" "${chrFilehandler}" "${sumstats}" "${id}" "${p}" "${conditionalFile}" "${pthresh}" "${flank}" "${snpsManual}" "${flankManual}" "${chrManual}" "${startManual}" "${endManual}" &

# gap_gwm
trait="gap_gwm"
traitDescription="grey and white matter"
targetDir="results/${trait}/gwama/eur/locuszoom/"
sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz"
conditionalFile="results/${trait}/gwama/eur/conditional/conditional.cleaned.tophits.annovar.txt"

	# optional settings
	snpsManual="rs55775495 rs1367858 rs530464314 rs2548331 rs2253491 rs9398805 rs76839250 17:19824272_CT_C rs2289629 rs4327091"
	flankManual=""
	chrManual="2 2 2 5 6 6 16 17 17 17"
	startManual="197800000 200200000 203000000 71600000 30500000 126200000 89400000 19400000 27400000 43000000"
	endManual="199200000 201600000 205000000 72800000 32000000 127600000 90400000 20600000 28800000 45500000"
	./code/genetics/locuszoom.sh "${trait}" "${traitDescription}" "${targetDir}" "${chrFilehandler}" "${sumstats}" "${id}" "${p}" "${conditionalFile}" "${pthresh}" "${flank}" "${snpsManual}" "${flankManual}" "${chrManual}" "${startManual}" "${endManual}"

# =================================================================
# === Discovery + replication: Get 95% credible set of variants ===
# =================================================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage/
conda activate envs/finemap

# run SBayesRC
cols="ID,A1,A2,A1_FREQ,BETA,SE,P,N" # set input column names that correspond to header columns SNP A1 A2 freq b se p N
ldmFolder="/fast/software/gctb/resources/ukbEUR_32k/ldm"
annotFile="/fast/software/gctb/resources/ukbEUR_32k/ldm/annot_baselineLF_v2.2.UKB_32k.txt"
threads=30
sbayes="sbayesrc"
imputation=1
mcmc=1
trait="gap_gm"; targetDir="results/${trait}/gwama/eur/sbayes"; sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz";
	taskset -c 0-29 ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}" "${mcmc}"
trait="gap_wm"; targetDir="results/${trait}/gwama/eur/sbayes"; sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz";
	taskset -c 40-79 ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}" "${mcmc}"
trait="gap_gwm"; targetDir="results/${trait}/gwama/eur/sbayes"; sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz";
	taskset -c 80-109 ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}" "${mcmc}"

	# get credible sets
	ldFile="/fast/software/gctb/resources/ukbEUR_32k/ldm/pairwise0.5.ld.txt" # file with pairwise ld > 0.5
	pip=0.95 # posterior inclusion probability
	pep=0.50 # posterior exclusion probability (proportion of random sets of variants with similar h2snp)
	threads=30 # number of threads to use
	pthresh=5e-8 # pvalue threshold for assigning discovered loci
	(for trait in gap_gm gap_wm gap_gwm; do (
		out="results/${trait}/gwama/eur/sbayes/sbayesrc" # output path
		sbayesRes="results/${trait}/gwama/eur/sbayes/sbayesrc" # path to sbayesrc results
		conditionalFile="results/${trait}/gwama/eur/conditional/conditional.cleaned.tophits.txt" # set conditional file to assign discovered loci
		code/genetics/sbayes.cs.sh "${out}" "${sbayesRes}" "${ldFile}" "${pip}" "${pep}" "${threads}" "${conditionalFile}" "${pthresh}"
		) &
	done
	wait)

	# map genes by position
	chr_bp_ref_alt="Chrom,Position,A1,A2"
	humandb="data/annovar/humandb"
	refseq="data/refseq/GRCh37_latest_genomic.edit.gff.gz"
	(for trait in gap_gm gap_wm gap_gwm; do (
		inputFile="results/${trait}/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.snpRes.test"
		outputFile="results/${trait}/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.snpRes.test.annot.txt"
		code/genetics/annotation.sh "${inputFile}" "${outputFile}" "${chr_bp_ref_alt}" "${humandb}" "${refseq}"
	) &
	done
	wait)

	# add relevant genes (sorted by cumulative PP of positionally mapped variants) to summary file
	(for trait in gap_gm gap_wm gap_gwm; do (
		annotFile="results/${trait}/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.snpRes.annot.txt"
		summaryFile="results/${trait}/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.joined"
		outputFile="results/${trait}/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.joined.genes.txt"
		code/genetics/sbayes.genes.R "${annotFile}" "${summaryFile}" "${outputFile}"
	) &
	done
	wait)

	# create summary file for exonic variants
	locusId="LEAD_SNP"
	id="Name"
	beta_se="A1Effect,SE"
	pthresh=1
	(for trait in gap_gm gap_wm gap_gwm; do (
		nsynFile="results/${trait}/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.snpRes.annot.nonsynonymous.txt"
		outputFile="results/${trait}/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.snpRes.annot.nonsynonymous.summary.txt"
		code/genetics/annotation.exonic.R "${nsynFile}" "${locusId}" "${id}" "${beta_se}" "${pthresh}" "${outputFile}" 
	) &
	done
	wait)

# run SusieR
bedFiles='data/genetics/chr${i}/imp_mri_qc_EURjoined/bed/chr${i}_mri_qc'
bgenFiles='data/genetics/chr${i}/imp_mri_qc_EURjoined/bgen/chr${i}_mri_qc'
LDsample="data/genetics/chr1/imp_mri_qc_EURjoined/chr1_mri_qc.psam"
sumstatsCols="CHR,BP,ID,A1,A2,A1_FREQ,BETA,SE,P,N" # requires the following columns: "CHR,POS,SNPID,A1,A2,A1_FREQ,BETA,SE,P,N" 
pthresh=5E-8
nthreads=100
(for trait in gap_gm gap_wm gap_gwm; do (
	rm -rf "results/${trait}/gwama/eur/susieR"
	targetDir="results/${trait}/gwama/eur/susieR"
	conditionalFile="results/${trait}/gwama/eur/conditional/conditional.cleaned.tophits.annovar.txt"
	sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz"
	taskset -c 0-99 ./code/genetics/susieR.sh "${trait}" "${targetDir}" "${bedFiles}" "${bgenFiles}" "${LDsample}" "${conditionalFile}" "${sumstats}" "${sumstatsCols}" "${pthresh}" "${nthreads}"
	) &
done
wait)

	# map genes by position
	chr_bp_ref_alt="chromosome,position,allele1,allele2"
	humandb="data/annovar/humandb"
	refseq="data/refseq/GRCh37_latest_genomic.edit.gff.gz"
	(for trait in gap_gm gap_wm gap_gwm; do (
		inputFile="results/${trait}/gwama/eur/susieR/susieR.df.purity50.txt"
		outputFile="results/${trait}/gwama/eur/susieR/susieR.df.purity50.annot.txt"
		code/genetics/annotation.sh "${inputFile}" "${outputFile}" "${chr_bp_ref_alt}" "${humandb}" "${refseq}"
	) &
	done
	wait)

	# add relevant genes (sorted by cumulative PP of variants from finemapping) to credible.ls.txt file
	(for trait in gap_gm gap_wm gap_gwm; do (
		dfFile="results/${trait}/gwama/eur/susieR/susieR.df.purity50.annot.txt"
		lsFile="results/${trait}/gwama/eur/susieR/susieR.ls.purity50.txt"
		outputFile="results/${trait}/gwama/eur/susieR/susieR.ls.purity50.genes.txt"
		code/genetics/susieR.genes.R "${dfFile}" "${lsFile}" "${outputFile}"
	) &
	done
	wait)

	# create summary file for exonic variants
	locusId="index_rsid"
	id="rsid"
	beta_se="beta,se"
	pthresh=5e-8
	(for trait in gap_gm gap_wm gap_gwm; do (
		nsynFile="results/${trait}/gwama/eur/susieR/susieR.df.purity50.annot.nonsynonymous.txt"
		outputFile="results/${trait}/gwama/eur/susieR/susieR.df.purity50.annot.nonsynonymous.summary.txt"
		code/genetics/annotation.exonic.R "${nsynFile}" "${locusId}" "${id}" "${beta_se}" "${pthresh}" "${outputFile}" 
	) &
	done
	wait)

# run FINEMAP
bedFiles='data/genetics/chr${i}/imp_mri_qc_EURjoined/bed/chr${i}_mri_qc'
bgenFiles='data/genetics/chr${i}/imp_mri_qc_EURjoined/bgen/chr${i}_mri_qc'
LDsample="data/genetics/chr1/imp_mri_qc_EURjoined/chr1_mri_qc.psam"
sumstatsCols="CHR,BP,ID,A1,A2,A1_FREQ,BETA,SE,P,N" # requires the following columns: "CHR,POS,SNPID,A1,A2,A1_FREQ,BETA,SE,P,N" 
pthresh=5E-8
nthreads=100
(for trait in gap_gm gap_wm gap_gwm; do (
	rm -rf "results/${trait}/gwama/eur/finemap_backup"
	\mv "results/${trait}/gwama/eur/finemap" "results/${trait}/gwama/eur/finemap_backup"
	targetDir="results/${trait}/gwama/eur/finemap"
	conditionalFile="results/${trait}/gwama/eur/conditional/conditional.cleaned.tophits.annovar.txt"
	sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz"
	taskset -c 0-99 ./code/genetics/finemap.sh "${trait}" "${targetDir}" "${bedFiles}" "${bgenFiles}" "${LDsample}" "${conditionalFile}" "${sumstats}" "${sumstatsCols}" "${pthresh}" "${nthreads}"
	) &
done
wait)

	# map genes by position
	chr_bp_ref_alt="chromosome,bp,allele1,allele2"
	humandb="data/annovar/humandb"
	refseq="data/refseq/GRCh37_latest_genomic.edit.gff.gz"
	(for trait in gap_gm gap_wm gap_gwm; do (
		inputFile="results/${trait}/gwama/eur/finemap/finemap.df.txt"
		outputFile="results/${trait}/gwama/eur/finemap/finemap.df.annot.txt"
		code/genetics/annotation.sh "${inputFile}" "${outputFile}" "${chr_bp_ref_alt}" "${humandb}" "${refseq}"
	) &
	done
	wait)

	# add relevant genes (sorted by cumulative PP of variants from finemapping) to credible.ls.txt file
	(for trait in gap_gm gap_wm gap_gwm; do (
		dfFile="results/${trait}/gwama/eur/finemap/finemap.df.annot.txt"
		lsFile="results/${trait}/gwama/eur/finemap/finemap.ls.txt"
		outputFile="results/${trait}/gwama/eur/finemap/finemap.ls.genes.txt"
		code/genetics/finemap.genes.R "${dfFile}" "${lsFile}" "${outputFile}"
	) &
	done
	wait)

	# create summary file for exonic variants
	locusId="index_rsid"
	id="rsid"
	beta_se="beta,se"
	pthresh=5e-8
	(for trait in gap_gm gap_wm gap_gwm; do (
		nsynFile="results/${trait}/gwama/eur/finemap/finemap.df.annot.nonsynonymous.txt"
		outputFile="results/${trait}/gwama/eur/finemap/finemap.df.annot.nonsynonymous.summary.txt"
		code/genetics/annotation.exonic.R "${nsynFile}" "${locusId}" "${id}" "${beta_se}" "${pthresh}" "${outputFile}" 
	) &
	done
	wait)

# ========================================================
# === Discovery + replication: run ld score regression ===
# ========================================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage
conda activate envs/ldsc

# ld score regression
sumstatsCols="CHR,BP,ID,A1,A2,BETA,SE,P,N"
sumstatsColNames="CHR,POS,SNP,A1,A2,BETA,SE,P,N"
LDsnplist="data/ldsc/resources/w_hm3.noMHC.snplist"
LDchr="data/ldsc/resources/eur_w_ld_chr/"
LDbaselineH2="data/ldsc/resources/1000G_Phase3_baselineLD_v2.2/baselineLD."
LDbaselineTau="data/ldsc/resources/1000G_Phase3_baseline_v1.2/baseline."
LDweights="data/ldsc/resources/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
LDfrqfiles="data/ldsc/resources/1000G_Phase3_frq/1000G.EUR.QC."
partitioned=1
celltypes=""
celltypegroups=$(echo /fast/software/ldsc/resources/1000G_Phase3_cell_type_groups/cell_type_group.{1..10}. | sed 's/ /,/g')
ctglabels="adrenal.pancreas,cardiovascular,cns,connective.bone,gi,hematopoietic,kidney,liver,other,skeletalmuscle"
(for trait in gap_gm gap_wm gap_gwm; do (
    targetDir="results/${trait}/gwama/eur/ldsc"
    rm -rf $targetDir
    sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz"
	./code/genetics/ldsc.sh "${trait}" "${targetDir}" "${sumstats}" "${sumstatsCols}" "${sumstatsColNames}" "${LDsnplist}" "${LDchr}" "${LDbaselineH2}" "${LDbaselineTau}" "${LDweights}" "${LDfrqfiles}" "${partitioned}" "${celltypes}" "${celltypegroups}" "${ctglabels}"
	) &
done
wait)

	# combine ldsc estimates in one file
	conda activate envs/default
	traitList="gap_gm gap_wm gap_gwm"
	ldscFiles="results/gap_gm/gwama/eur/ldsc/ldsc.h2.results results/gap_wm/gwama/eur/ldsc/ldsc.h2.results results/gap_gwm/gwama/eur/ldsc/ldsc.h2.results"
	outFile="results/combined/gwama.eur.ldsc.h2.txt"
	./code/genetics/ldsc.combine.sh "${ldscFiles}" "${outFile}"

	# make partitioned heritability plot and get summary table
	ymin=-10
	ymax=25
	ysteps=10
	width=8.1
	height=3.7
	oneTailed=TRUE
	(for trait in gap_gm gap_wm gap_gwm; do (
	    ldscFile="results/${trait}/gwama/eur/ldsc/ldsc.partitioned.results"
	    Rscript code/genetics/ldsc.partitioned.R "${ldscFile}" "${oneTailed}" "${ymin}" "${ymax}" "${ysteps}" "${width}" "${height}"
	    ) &
	done
	wait)

		# combine plots and summary tables
		traits="gap_gm,gap_wm,gap_gwm"
		plotTitles="Grey_matter,White_matter,Grey_and_white_matter"
		ldscPlots="results/gap_gm/gwama/eur/ldsc/ldsc.partitioned.results.png,results/gap_wm/gwama/eur/ldsc/ldsc.partitioned.results.png,results/gap_gwm/gwama/eur/ldsc/ldsc.partitioned.results.png"
		ldscSummary="results/gap_gm/gwama/eur/ldsc/ldsc.partitioned.results.summary.txt,results/gap_wm/gwama/eur/ldsc/ldsc.partitioned.results.summary.txt,results/gap_gwm/gwama/eur/ldsc/ldsc.partitioned.results.summary.txt"
		outputFile="results/combined/gwama.eur.ldsc.partitioned"
		width=8.1
		height=11.1
		Rscript code/genetics/ldsc.partitioned.combined.R "${traits}" "${plotTitles}" "${ldscPlots}" "${ldscSummary}" "${outputFile}" "${width}" "${height}"

	# make cell-type-group plot and get summary table
	ymin=0
	ymax=6
	ysteps=2
	width=5.1
	height=3.0
	oneTailed=TRUE
	(for trait in gap_gm gap_wm gap_gwm; do (
	    ldscFile="results/${trait}/gwama/eur/ldsc/ldsc.ctg.results"
	    Rscript code/genetics/ldsc.ctg.R "${ldscFile}" "${oneTailed}" "${ymin}" "${ymax}" "${ysteps}" "${width}" "${height}"
	    ) &
	done
	wait)

		# combine plots and summary tables
		traits="gap_gm,gap_wm,gap_gwm"
		plotTitles="Grey_matter,White_matter,Grey_and_white_matter"
		ldscPlots="results/gap_gm/gwama/eur/ldsc/ldsc.ctg.results.png,results/gap_wm/gwama/eur/ldsc/ldsc.ctg.results.png,results/gap_gwm/gwama/eur/ldsc/ldsc.ctg.results.png"
		ldscSummary="results/gap_gm/gwama/eur/ldsc/ldsc.ctg.results.summary.txt,results/gap_wm/gwama/eur/ldsc/ldsc.ctg.results.summary.txt,results/gap_gwm/gwama/eur/ldsc/ldsc.ctg.results.summary.txt"
		outputFile="results/combined/gwama.eur.ldsc.ctg"
		width=5.1
		height=9
		Rscript code/genetics/ldsc.partitioned.combined.R "${traits}" "${plotTitles}" "${ldscPlots}" "${ldscSummary}" "${outputFile}" "${width}" "${height}"

	# get genetic correlation beween traits
	conda activate envs/ldsc
	LDchr="data/ldsc/resources/eur_w_ld_chr/"
	xprefix=""
	yprefix=""
	threads=100 # set number of parallel workers
	sumstatsX="results/gap_gm/gwama/eur/ldsc/ldsc.gap_gm.sumstats.gz results/gap_wm/gwama/eur/ldsc/ldsc.gap_wm.sumstats.gz results/gap_gwm/gwama/eur/ldsc/ldsc.gap_gwm.sumstats.gz"
	sumstatsY="${sumstatsX}"
	targetDir="results/combined/gwama.eur.rgtissue/"
	./code/genetics/rg.sh "${targetDir}" "${sumstatsX}" "${sumstatsY}" "${xprefix}" "${yprefix}" "${LDchr}" "${threads}"

	# for comparison, get phenotypic correlation
	R
	gm = read.delim('data/traits/gap_gm.txt', header = T)
	wm = read.delim('data/traits/gap_wm.txt', header = T)
	gwm = read.delim('data/traits/gap_gwm.txt', header = T)
	discov = dplyr::inner_join(gm,wm)
	discov = dplyr::inner_join(discov,gwm)
	replic = read.delim('data/traits/replicate.txt', header = T)
	replic = replic[replic$pan=='EUR' & !is.na(replic$pan),c('FID','IID','gap_gm','gap_wm','gap_gwm')]
	df = rbind(discov,replic)
	covs = read.delim('results/mri/r2024.vars.txt', header = T)
	covs = covs[,c('FID','IID','sex','t1.age','t1.age2','t1.ac1','t1.ac2','t1.ac3','t1.TIV')]
	df = dplyr::inner_join(df,covs)
	pcor_result = ppcor::pcor.test(df$gap_gm,df$gap_wm,df[,6:10])
		r = pcor_result$estimate
		k = pcor_result$gp
		n = nrow(df)
		se = sqrt((1-r^2)/(n-k-1))
		message(sprintf('r_gm_wm: %.3f, se_gm_wm: %.3f',r,se))
	pcor_result = ppcor::pcor.test(df$gap_gm,df$gap_gwm,df[,6:10])
		r = pcor_result$estimate
		k = pcor_result$gp
		n = nrow(df)
		se = sqrt((1-r^2)/(n-k-1))
		message(sprintf('r_wm_gwm: %.3f, se_wm_gwm: %.3f',r,se))
	pcor_result = ppcor::pcor.test(df$gap_wm,df$gap_gwm,df[,6:10])
		r = pcor_result$estimate
		k = pcor_result$gp
		n = nrow(df)
		se = sqrt((1-r^2)/(n-k-1))
		message(sprintf('r_wm_gwm: %.3f, se_wm_gwm: %.3f',r,se))
	q()

# ====================================================================
# === Discovery + replication: Run ANNOVAR enrichment test in FUMA ===
# ====================================================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage
conda activate envs/default

# prepare summary statistics for FUMA upload
sumstatsCols="CHR,BP,ID,A1,A2,BETA,SE,P,N" # first column must be chromosome column due to filtering for Y, XY, and MT
sumstatsColNames="CHR,BP,SNP,A1,A2,BETA,SE,P,N" # column names must be in accordance with fuma input https://fuma.ctglab.nl/tutorial#prepare-input-files
(for trait in gap_gm gap_wm gap_gwm; do (
	targetDir="results/${trait}/gwama/eur/fuma"
	sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz"
    ./code/genetics/fuma.sh "${trait}" "${targetDir}" "${sumstats}" "${sumstatsCols}" "${sumstatsColNames}"
    ) &
done
wait)

# combine results of ANNOVAR enrichment test
traits="gap_gm,gap_wm,gap_gwm"
inFiles="results/gap_gm/gwama/eur/fuma/FUMA_job525095/annov.stats.txt,results/gap_wm/gwama/eur/fuma/FUMA_job525096/annov.stats.txt,results/gap_gwm/gwama/eur/fuma/FUMA_job525098/annov.stats.txt"
outFile="results/combined/gwama.eur.annovar.txt"
./code/genetics/fuma.annov.combine.R "${traits}" "${inFiles}" "${outFile}"

# make annovar enrichment plot for all traits
annovarStats="results/gap_gm/gwama/eur/fuma/FUMA_job525095.zip,results/gap_wm/gwama/eur/fuma/FUMA_job525096.zip,results/gap_gwm/gwama/eur/fuma/FUMA_job525098.zip"
outputFile="results/combined/gwama.eur.annovar.png"
width=8.04
height=7.50
Rscript code/genetics/fuma.annov.R "${annovarStats}" "${outputFile}" "${width}" "${height}"

# =====================================================
# === Discovery + replication: run ApoE-e4 analysis ===
# =====================================================

# set working directory and load conda environment
cd /slow/projectss/ukb_brainage
conda activate envs/default

# run apoe-e4 association
geneticsDir="data/genetics"
maf=0.01
threads=50
(for trait in gap_gm gap_wm gap_gwm; do (
	targetDir="results/${trait}/apoe"
	./code/genetics/apoe.sh "${trait}" "${targetDir}" "${geneticsDir}" "${maf}" "${threads}" 
	) &
done)
wait

# ===========================================================
# === Discovery + replication: run eqtl/sqtl SMR analysis ===
# ===========================================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage/
conda activate envs/default

# run analysis
chrFilehandler='data/genetics/chr${i}/imp_mri_qc_EURjoined/bed/chr${i}_mri_qc'
brainmetaFolder="data/smr/BrainMeta"
sumstatsCols="ID,A1,A2,A1_FREQ,BETA,SE,P,N" # "SNP,A1,A2,freq,b,se,p,N"
pthresh=1E-6
LDsample="data/genetics/chr1/imp_mri_qc_EURjoined/chr1_mri_qc.psam"
threads=100
clumpR2=0.8
smrCorrect="fdr"

for trait in gap_gm gap_wm gap_gwm; do
	targetDir="results/${trait}/gwama/eur/smr"
	sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz"
	conditionalFile="results/${trait}/gwama/eur/conditional/conditional.cleaned.tophits.annovar.txt"
	./code/genetics/smr.sh "${trait}" "${targetDir}" "${chrFilehandler}" "${brainmetaFolder}" "${sumstats}" "${sumstatsCols}" "${conditionalFile}" "${pthresh}" "${LDsample}" "${threads}" "${clumpR2}" "${smrCorrect}" 
done

# ======================================================================
# === Discovery + replication: run GWAS catalog screening on tophits ===
# ======================================================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage/
conda activate envs/default

# run catalog screening
cols2keep="LOCUS_CNT,SNP_CNT,CHR,BP,ID,CYTOBAND,LEAD_SNP,LEAD_SNP_P,LEAD_SNP_KB,LEAD_SNP_RSQ,LEAD_SNP_ALLELES,REGION,NEAREST_GENE,NEAREST_GENE_DESCRIPTION,NEAREST_GENE_BIOTYPE,DISTANCE,A1,A2,A1_FREQ,BETA,SE,P,N"
snpCol="ID"
catalogFile="data/gwas_catalog/gwas_catalog_v1.0-associations_e112_r2024-09-13.tsv"
catalogThresh="5E-8"
(for trait in gap_gm gap_wm gap_gwm; do (
	inFile="results/${trait}/gwama/eur/conditional/conditional.cleaned.tophits.annovar.txt"
	outFile="results/${trait}/gwama/eur/catalog/catalog.txt"
	mkdir -p "results/${trait}/gwama/eur/catalog/"
	awk -F'\t' 'NR==1 { print } $16 < 5E-8 && $35 < 5E-8 { print }' OFS='\t' ${inFile} > ${inFile}.tmp
	./code/genetics/catalog.sh "${inFile}.tmp" "${outFile}" "${cols2keep}" "${snpCol}" "${catalogFile}" "${catalogThresh}" # only perform direct variant lookups
	rm -f ${inFile}.tmp
	) &
done
wait)

# create per-locus-summary file
locusCol="LOCUS_CNT"
leadsnpCol="LEAD_SNP"
catalogTraitCol="catalog_trait"
(for trait in gap_gm gap_wm gap_gwm; do (
	inFile="results/${trait}/gwama/eur/catalog/catalog.txt"
	outFile="results/${trait}/gwama/eur/catalog/catalog"
	./code/genetics/catalog.output.R "${inFile}" "${outFile}" "${locusCol}" "${leadsnpCol}" "${catalogTraitCol}"
	) &
done
wait)

# rearrange columns
compress="none"
skipLines=0
sep="auto"
colsIn="LOCUS_CNT,SNP_CNT,CHR,BP,ID,CYTOBAND,LEAD_SNP,LEAD_SNP_P,LEAD_SNP_KB,LEAD_SNP_RSQ,LEAD_SNP_ALLELES,REGION,NEAREST_GENE,NEAREST_GENE_DESCRIPTION,NEAREST_GENE_BIOTYPE,DISTANCE,A1,A2,A1_FREQ,BETA,SE,P,N,catalog_trait,catalog_pubmedid,catalog_author,catalog_date,catalog_a1,catalog_pval"
colsRename="${colsIn}"
colsOut="${colsIn}"
(for trait in gap_gm gap_wm gap_gwm; do (
	inFile="results/${trait}/gwama/eur/catalog/catalog.txt"
	outFile="results/${trait}/gwama/eur/catalog/catalog.suppl.txt"
	./code/genetics/arrange.R "${inFile}" "${outFile}" "${compress}" "${skipLines}" "${sep}" "${colsIn}" "${colsRename}" "${colsOut}"
	) &
done
wait)

# =====================================================
# === Discovery + replication: run GTEx eQTL lookup ===
# =====================================================

# set working directory
cd /slow/projects/ukb_brainage/
conda activate envs/eqtl

# run gtex eqtl lookup
rsid="ID"
chr="CHR"
bp="BP"
ref="A2"
alt="A1"
locus="LOCUS_CNT"
assoc="SNP_CNT"
leadsnp="LEAD_SNP"
gtexDir="data/gtex"
(for trait in gap_gm gap_wm gap_gwm; do (
	inFile="results/${trait}/gwama/eur/conditional/conditional.cleaned.tophits.annovar.txt"
	outFile="results/${trait}/gwama/eur/eqtl/eqtl"
	mkdir -p "results/${trait}/gwama/eur/eqtl/"
	awk -F'\t' 'NR==1 { print } $16 < 5E-8 && $35 < 5E-8 { print }' OFS='\t' ${inFile} > ${inFile}.tmp
	./code/genetics/eqtl.sh "${inFile}.tmp" "${outFile}" "${rsid}" "${chr}" "${bp}" "${ref}" "${alt}" "${locus}" "${assoc}" "${leadsnp}" "${gtexDir}"
	rm -f ${inFile}.tmp
	) &
done
wait)

# ===========================================================================
# === Discovery + replication: Calculate Polygenic Priority Scores (PoPS) ===
# ===========================================================================

# set working directory
cd "/slow/projects/ukb_brainage"
conda activate envs/pops

# run magma
(for trait in gap_gm gap_wm gap_gwm; do (
	mkdir -p results/${trait}/gwama/eur/magma
	zcat results/${trait}/gwama/eur/metal.ivweight.qc.gz > results/${trait}/gwama/eur/magma/sumstats.txt
	for chr in {1..22}; do (
		magma --bfile data/genetics/chr#CHR#/imp_mri_qc_EURjoined/bed/chr#CHR#_mri_qc \
			--batch ${chr} chr \
			--gene-annot /fast/software/pops/data/utils/magma_0kb.genes.annot \
			--pval results/${trait}/gwama/eur/magma/sumstats.txt use=ID,P ncol=N \
			--out results/${trait}/gwama/eur/magma/magma
			) &
	done
	wait
	magma --merge results/${trait}/gwama/eur/magma/magma --out results/${trait}/gwama/eur/magma/magma
	rm -f results/${trait}/gwama/eur/magma/sumstats.txt
	rm -f results/${trait}/gwama/eur/magma/magma.batch*
	chmod 750 results/${trait}/gwama/eur/magma/*) &
done
wait)

# make annotations
annotationFile="/fast/software/pops/data/utils/gene_annot_jun10.txt"
pthreshMapping=5E-8
glist_hg19="data/glist_hg19/glist-hg19-refseq.txt"
threads=100
window=0
humandb="data/annovar/humandb"
refseq="data/refseq/GRCh37_latest_genomic.edit.gff.gz"
clumpingWindow=3000
(for trait in gap_gm gap_wm gap_gwm; do (
	inputFile="results/${trait}/gwama/eur/magma/magma.genes.out"
	outputFile="results/${trait}/gwama/eur/magma/magma.genes.out.annot"
	conditionalFile="results/${trait}/gwama/eur/conditional/conditional.cleaned.tophits.annovar.txt"
	./code/genetics/magma.sh "${inputFile}" "${outputFile}" "${annotationFile}" "${conditionalFile}" "${pthreshMapping}" "${glist_hg19}" "${threads}" "${window}" "${humandb}" "${refseq}" "${clumpingWindow}" 
) &
done
wait)

# calculate polygenic priority score (PoPS)
annotationFile="/fast/software/pops/data/utils/gene_annot_jun10.txt"
featurePrefix="/fast/software/pops/data/features_munged/pops_features"
nChunks=115
controlFeatures="/fast/software/pops/data/utils/control.features"
(for trait in gap_gm gap_wm gap_gwm; do (
	mkdir -p results/${trait}/gwama/eur/pops
	python /fast/software/pops/pops.py \
		--gene_annot_path ${annotationFile} \
		--feature_mat_prefix ${featurePrefix} \
		--num_feature_chunks ${nChunks} \
		--magma_prefix results/${trait}/gwama/eur/magma/magma \
		--control_features_path ${controlFeatures} \
		--out_prefix results/${trait}/gwama/eur/pops/pops
	awk 'NR==1 { next } NR==FNR { symbol[$1]=$2; next }
			 FNR==1 { print "SYMBOL",$0; next }
			 $1 in symbol { print symbol[$1], $0 }' OFS='\t' ${annotationFile} results/${trait}/gwama/eur/pops/pops.preds \
			 > results/${trait}/gwama/eur/pops/pops.preds.symbols
	chmod 750 results/${trait}/gwama/eur/pops/*
	) &
done
wait)

# combine results
conda activate envs/default
traits="gap_gm,gap_wm,gap_gwm"
magmaFiles="results/gap_gm/gwama/eur/magma/magma.genes.out.annot.clumped,results/gap_wm/gwama/eur/magma/magma.genes.out.annot.clumped,results/gap_gwm/gwama/eur/magma/magma.genes.out.annot.clumped"
popsFiles="results/gap_gm/gwama/eur/pops/pops.preds,results/gap_wm/gwama/eur/pops/pops.preds,results/gap_gwm/gwama/eur/pops/pops.preds"
clumpingWindow=3000
outputFile="results/combined/gwama.eur.pops.txt"
./code/genetics/pops.combine.R "${traits}" "${magmaFiles}" "${popsFiles}" "${clumpingWindow}" "${outputFile}"

# prioritize genes based on GWAS index variants
conditionalPthresh=5E-8
annotationFile="/fast/software/pops/data/utils/gene_annot_jun10.txt"
window=500
(for trait in gap_gm gap_wm gap_gwm; do (
	conditionalFile="results/${trait}/gwama/eur/conditional/conditional.cleaned.txt"
	popsFile="results/${trait}/gwama/eur/pops/pops.preds.symbols"
	outputFile="results/${trait}/gwama/eur/pops/pops.prio.txt"
	./code/genetics/pops.prioritize.R "${conditionalFile}" "${conditionalPthresh}" "${popsFile}" "${annotationFile}" "${outputFile}" "${window}"
) &
done)
wait

# =========================================================================================
# === Discovery + replication: identify independent snp-level discoveries across traits ===
# =========================================================================================

# set working directory
cd /slow/projects/ukb_brainage/
conda activate envs/default

# settings
traitlist="gap_gm gap_wm gap_gwm"
conditionalFiles="results/gap_gm/gwama/eur/conditional/conditional.cleaned.tophits.annovar.txt results/gap_wm/gwama/eur/conditional/conditional.cleaned.tophits.annovar.txt results/gap_gwm/gwama/eur/conditional/conditional.cleaned.tophits.annovar.txt"
nonsynFiles="results/gap_gm/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.snpRes.annot.nonsynonymous.summary.txt results/gap_wm/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.snpRes.annot.nonsynonymous.summary.txt results/gap_gwm/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.snpRes.annot.nonsynonymous.summary.txt"
catalogFiles="results/gap_gm/gwama/eur/catalog/catalog.by.locus.txt results/gap_wm/gwama/eur/catalog/catalog.by.locus.txt results/gap_gwm/gwama/eur/catalog/catalog.by.locus.txt"
susieRFiles="results/gap_gm/gwama/eur/susieR/susieR.ls.purity50.genes.txt results/gap_wm/gwama/eur/susieR/susieR.ls.purity50.genes.txt results/gap_gwm/gwama/eur/susieR/susieR.ls.purity50.genes.txt"
sbayesFiles="results/gap_gm/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.joined.genes.txt results/gap_wm/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.joined.genes.txt results/gap_gwm/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.joined.genes.txt"
finemapFiles="results/gap_gm/gwama/eur/finemap/finemap.ls.genes.txt results/gap_wm/gwama/eur/finemap/finemap.ls.genes.txt results/gap_gwm/gwama/eur/finemap/finemap.ls.genes.txt"
smrFiles="results/gap_gm/gwama/eur/smr/smr.summary.txt results/gap_wm/gwama/eur/smr/smr.summary.txt results/gap_gwm/gwama/eur/smr/smr.summary.txt"
eqtlFiles="results/gap_gm/gwama/eur/eqtl/eqtl.summary.txt results/gap_wm/gwama/eur/eqtl/eqtl.summary.txt results/gap_gwm/gwama/eur/eqtl/eqtl.summary.txt"
popsFiles="results/gap_gm/gwama/eur/pops/pops.prio.txt results/gap_wm/gwama/eur/pops/pops.prio.txt results/gap_gwm/gwama/eur/pops/pops.prio.txt"
targetDir="results/combined/gwama.eur.snplevel/"
chrFilehandler='data/genetics/chr${i}/imp_mri_qc_EURjoined/bed/chr${i}_mri_qc'
LDsample="data/genetics/chr1/imp_mri_qc_EURjoined/chr1_mri_qc.psam"
pthresh=1E-6
./code/genetics/snplevel.sh "${traitlist}" "${conditionalFiles}" "${nonsynFiles}" "${catalogFiles}" "${sbayesFiles}" "${susieRFiles}" "${finemapFiles}" "${smrFiles}" "${eqtlFiles}" "${popsFiles}" "${targetDir}" "${chrFilehandler}" "${LDsample}" "${pthresh}" 

# add gene prioritization score
inputFile="results/combined/gwama.eur.snplevel/snplevel.gws.txt"
cols="sbayes_genes nonsyn_customPthresh SMR_eqtl SMR_sqtl GTEx_singleTissue GTEx_multiTissue PoPS"
outputFile="results/combined/gwama.eur.snplevel/snplevel.gws.prio.txt"
./code/genetics/prioritize.R "${inputFile}" "${cols}" "${outputFile}"

# =========================================================================================
# === Discovery + replication: check novelty by comparing own with previous discoveries ===
# =========================================================================================

# set working directory
cd /slow/projects/ukb_brainage/
conda activate envs/default

# settings
ownDiscoveries="results/combined/gwama.eur.snplevel/snplevel.gws.prio.txt"
prevDiscoveries="data/prevDiscoveries/Jonsson_2019_munged.txt,data/prevDiscoveries/Ning_2021_munged.txt,data/prevDiscoveries/Smith_2020_munged.txt,data/prevDiscoveries/Kim_2023_munged.txt,data/prevDiscoveries/Leonardsen_2023_munged.txt,data/prevDiscoveries/Wen_2024_munged.txt"
targetDir="results/combined/gwama.eur.snplevel/"
chrFilehandler='data/genetics/2024/chr${i}/imp_mri/chr${i}_mri'
LDsample="data/genetics/chr1/imp_mri_qc_EURjoined/chr1_mri_qc.psam"
clumpkb=10000
clumpr2=0.1

# run analysis
./code/genetics/snplevel.novelty.sh "${ownDiscoveries}" "${prevDiscoveries}" "${targetDir}" "${chrFilehandler}" "${LDsample}" "${clumpkb}" "${clumpr2}"

# ============================================================================================
# === Discovery + replication: make GWAS discovery table for main article and supplementum ===
# ============================================================================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage
conda activate envs/default

# settings
noveltyFile="results/combined/gwama.eur.snplevel/snplevel.novelty.txt"
replicationFile=""
outputFileMain="results/combined/gwama.eur.snplevel/discoveries.main.txt"
outputFileSuppl="results/combined/gwama.eur.snplevel/discoveries.suppl.txt"
traits="gap_gm,gap_wm,gap_gwm"
traitNamesMain="GM,WM,GWM"
traitNamesSuppl="grey matter,white matter,grey and white matter"
Rscript code/genetics/snplevel.novelty.tables.R "${noveltyFile}" "${replicationFile}" "${outputFileMain}" "${outputFileSuppl}" "${traits}" "${traitNamesMain}" "${traitNamesSuppl}"

# ========================================================================================================
# === Discovery + replication: combine results based on discovered loci across traits for supplementum ===
# ========================================================================================================

# set working director and conda environment
cd /slow/projects/ukb_brainage/
conda activate envs/default

# combine SBayesRC finemapping results across traits
snplevelFile="results/combined/gwama.eur.snplevel/snplevel.gws.txt"
finemapFiles="results/gap_gm/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.snpRes.annot.txt,results/gap_wm/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.snpRes.annot.txt,results/gap_gwm/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.snpRes.annot.txt"
outputFile="results/combined/gwama.eur.credible.sbayesrc.suppl.txt"
traits="gap_gm,gap_wm,gap_gwm"
traitNames="grey matter,white matter,grey and white matter"
Rscript code/genetics/sbayes.cs.combine.R "${snplevelFile}" "${finemapFiles}" "${outputFile}" "${traits}" "${traitNames}"

# combine SusieR finemapping results across traits
# truncate at 100 variants per credible set
snplevelFile="results/combined/gwama.eur.snplevel/snplevel.gws.txt"
finemapFiles="results/gap_gm/gwama/eur/susieR/susieR.df.purity50.annot.txt,results/gap_wm/gwama/eur/susieR/susieR.df.purity50.annot.txt,results/gap_gwm/gwama/eur/susieR/susieR.df.purity50.annot.txt"
outputFile="results/combined/gwama.eur.credible.susieR.suppl.txt"
traits="gap_gm,gap_wm,gap_gwm"
traitNames="grey matter,white matter,grey and white matter"
Rscript code/genetics/susieR.cs.combine.R "${snplevelFile}" "${finemapFiles}" "${outputFile}" "${traits}" "${traitNames}"
awk -F'\t' 'BEGIN { count=0 } seen[$1":"$2":"$3":"$11]++ < 100 { print }' OFS='\t' "results/combined/gwama.eur.credible.susieR.suppl.txt" > "results/combined/gwama.eur.credible.susieR.suppl.truncated.txt"

# combine FINEMAP results across traits
# truncate at 100 variants per credible set
snplevelFile="results/combined/gwama.eur.snplevel/snplevel.gws.txt"
finemapFiles="results/gap_gm/gwama/eur/finemap/finemap.df.annot.txt,results/gap_wm/gwama/eur/finemap/finemap.df.annot.txt,results/gap_gwm/gwama/eur/finemap/finemap.df.annot.txt"
outputFile="results/combined/gwama.eur.credible.finemap.suppl.txt"
traits="gap_gm,gap_wm,gap_gwm"
traitNames="grey matter,white matter,grey and white matter"
Rscript code/genetics/finemap.combine.R "${snplevelFile}" "${finemapFiles}" "${outputFile}" "${traits}" "${traitNames}"
awk -F'\t' 'BEGIN { count=0 } seen[$1":"$2":"$3":"$13]++ < 100 { print }' OFS='\t' "results/combined/gwama.eur.credible.finemap.suppl.txt" > "results/combined/gwama.eur.credible.finemap.suppl.truncated.txt"

# compare credible sets produced by different fine-mapping techniques
sbayesfile="results/combined/gwama.eur.credible.sbayesrc.suppl.txt"
susieRfile="results/combined/gwama.eur.credible.susieR.suppl.txt"
finemapfile="results/combined/gwama.eur.credible.finemap.suppl.txt"
topsetsize=10
topsetoverlap=0
outputfile="results/combined/gwama.eur.credible.overlap"
Rscript code/genetics/cs.compare.R "${sbayesfile}" "${susieRfile}" "${finemapfile}" "${topsetsize}" "${topsetoverlap}" "${outputfile}" 
cat results/combined/gwama.eur.credible.overlap.topsets.txt | awk 'NR>1 && !seen[$1]++' | wc -l

# combine exonic results for supplementum
snplevelFile="results/combined/gwama.eur.snplevel/snplevel.gws.txt"
nsynFiles="results/gap_gm/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.snpRes.annot.nonsynonymous.txt,results/gap_wm/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.snpRes.annot.nonsynonymous.txt,results/gap_gwm/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.snpRes.annot.nonsynonymous.txt"
outputFile="results/combined/gwama.eur.nonsynonymous.suppl.txt"
traits="gap_gm,gap_wm,gap_gwm"
traitNames="grey matter,white matter,grey and white matter"
matchcol="LEAD_SNP"
outcols="LOCUS_COUNT,trait,LEAD_SNP,LEAD_SNP_CHR,LEAD_SNP_BP,Name,Position,A1,A2,A1Frq,A1Effect,SE,PIP,REGION,NEAREST_GENE,NEAREST_GENE_DESCRIPTION,NEAREST_GENE_BIOTYPE,EXONIC_FUNCTION,TRANSCRIPT_CONSEQUENCE,CADD_PHRED,CADD_RAW,CADD_RANK,DANN,DANN_RANK,REVEL,REVEL_RANK"
code/genetics/annotation.exonic.combine.R "${snplevelFile}" "${nsynFiles}" "${outputFile}" "${traits}" "${traitNames}" "${matchcol}" "${outcols}" 
awk -F'\t' 'NR==1 || $14<5E-8' "${outputFile}" > "${outputFile}".tmp
\mv "${outputFile}".tmp "${outputFile}"

# combine SMR eqtl results across traits
snplevelFile="results/combined/gwama.eur.snplevel/snplevel.gws.txt"
smrFiles="results/gap_gm/gwama/eur/smr/smr.eqtl.filtered.assigned.txt,results/gap_wm/gwama/eur/smr/smr.eqtl.filtered.assigned.txt,results/gap_gwm/gwama/eur/smr/smr.eqtl.filtered.assigned.txt"
outputFile="results/combined/gwama.eur.smr.eqtl.suppl.txt"
traits="gap_gm,gap_wm,gap_gwm"
traitNames="grey matter,white matter,grey and white matter"
Rscript code/combine.smr.R "${snplevelFile}" "${smrFiles}" "${outputFile}" "${traits}" "${traitNames}"

# combine SMR sqtl results across traits
snplevelFile="results/combined/gwama.eur.snplevel/snplevel.gws.txt"
smrFiles="results/gap_gm/gwama/eur/smr/smr.sqtl.filtered.assigned.txt,results/gap_wm/gwama/eur/smr/smr.sqtl.filtered.assigned.txt,results/gap_gwm/gwama/eur/smr/smr.sqtl.filtered.assigned.txt"
outputFile="results/combined/gwama.eur.smr.sqtl.suppl.txt"
traits="gap_gm,gap_wm,gap_gwm"
traitNames="grey matter,white matter,grey and white matter"
Rscript code/combine.smr.R "${snplevelFile}" "${smrFiles}" "${outputFile}" "${traits}" "${traitNames}"

# combine eQTL single-tissue results across traits
snplevelFile="results/combined/gwama.eur.snplevel/snplevel.gws.txt"
eqtlFiles="results/gap_gm/gwama/eur/eqtl/eqtl.singleTissue.txt,results/gap_wm/gwama/eur/eqtl/eqtl.singleTissue.txt,results/gap_gwm/gwama/eur/eqtl/eqtl.singleTissue.txt"
outputFile="results/combined/gwama.eur.gtex.singletissue.txt"
traits="gap_gm,gap_wm,gap_gwm"
traitNames="grey matter,white matter,grey and white matter"
Rscript code/genetics/eqtl.single.combine.R "${snplevelFile}" "${eqtlFiles}" "${outputFile}" "${traits}" "${traitNames}"

# combine eQTL multi-tissue results across traits
snplevelFile="results/combined/gwama.eur.snplevel/snplevel.gws.txt"
eqtlFiles="results/gap_gm/gwama/eur/eqtl/eqtl.multiTissue.txt,results/gap_wm/gwama/eur/eqtl/eqtl.multiTissue.txt,results/gap_gwm/gwama/eur/eqtl/eqtl.multiTissue.txt"
outputFile="results/combined/gwama.eur.gtex.multitissue.txt"
traits="gap_gm,gap_wm,gap_gwm"
traitNames="grey matter,white matter,grey and white matter"
Rscript code/genetics/eqtl.multi.combine.R "${snplevelFile}" "${eqtlFiles}" "${outputFile}" "${traits}" "${traitNames}"

# combine NHGRI GWAS-catalog results across traits
snplevelFile="results/combined/gwama.eur.snplevel/snplevel.gws.txt"
catalogFiles="results/gap_gm/gwama/eur/catalog/catalog.txt,results/gap_wm/gwama/eur/catalog/catalog.txt,results/gap_gwm/gwama/eur/catalog/catalog.txt"
outputFile="results/combined/gwama.eur.catalog.txt"
traits="gap_gm,gap_wm,gap_gwm"
traitNames="grey matter,white matter,grey and white matter"
Rscript code/genetics/catalog.combine.R "${snplevelFile}" "${catalogFiles}" "${outputFile}" "${traits}" "${traitNames}"

# ==================================================================
# === Discovery + replication: run gene-based analysis (fastBAT) ===
# ==================================================================

# set working directory and load conda environment
cd "/slow/projects/ukb_brainage"
conda activate envs/default

# run analysis
id="ID"
p="P"
chrFilehandler='data/genetics/chr${i}/imp_mri_qc_EURjoined/bed/chr${i}_mri_qc'
sampleFile="data/genetics/chr1/imp_mri_qc_EURjoined/chr1_mri_qc.psam"
pthreshMapping=5E-8
glist_hg19="data/glist_hg19/glist-hg19-refseq.txt"
threads=20
window=0
humandb="data/annovar/humandb"
refseq="data/refseq/GRCh37_latest_genomic.edit.gff.gz"
clumpingWindow=3000

	# for gcta software: change chromosome strings to numeric values
	for i in X Y MT; do \cp $(eval echo ${chrFilehandler}_numeric.bim) $(eval echo ${chrFilehandler}.bim); done
	i=XY; \cp $(eval echo ${chrFilehandler}_numeric_23.bim) $(eval echo ${chrFilehandler}.bim)

	# run analysis
	(for trait in gap_gm gap_wm gap_gwm; do (
		ulimit -Sv 50000000
		targetDir="results/${trait}/gwama/eur/fastbat"
		sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz"
		conditionalFile="results/${trait}/gwama/eur/conditional/conditional.cleaned.tophits.annovar.txt"
		./code/genetics/fastbat.sh "${targetDir}" "${chrFilehandler}" "${sampleFile}" "${sumstats}" "${id}" "${p}" "${conditionalFile}" "${pthreshMapping}" "${glist_hg19}" "${threads}" "${window}" "${humandb}" "${refseq}" "${clumpingWindow}"
		) &
	done)
	wait

	# replace .bim files with original files (change numeric values to chromosome strings)
	for i in X Y XY MT; do \cp $(eval echo ${chrFilehandler}_original.bim) $(eval echo ${chrFilehandler}.bim); done

# combine results
traits="gap_gm,gap_wm,gap_gwm"
inputFiles="results/gap_gm/gwama/eur/fastbat/fastbat.0kb.gene.fastbat.annovar.clumped,results/gap_wm/gwama/eur/fastbat/fastbat.0kb.gene.fastbat.annovar.clumped,results/gap_gwm/gwama/eur/fastbat/fastbat.0kb.gene.fastbat.annovar.clumped"
clumpingWindow=3000 # cross-trait clumping window size in kb
outputFile="results/combined/gwama.eur.fastbat.txt"
./code/genetics/fastbat.combine.R "${traits}" "${inputFiles}" "${clumpingWindow}" "${outputFile}"

# create gene-based Manhattan plots
annotationThresh=1e-100 # 'fdr' or 'bonferroni' or value
ylim=40 # upper y axis limit
ysteps=10 # y axis breaks
width=11 # plot width
height=4 # plot height
(for trait in gap_gm gap_wm gap_gwm; do (
	targetDir="results/${trait}/gwama/eur/fastbat/"
	fastbatResults="results/${trait}/gwama/eur/fastbat/fastbat.0kb.gene.fastbat.annovar.clumped"
	Rscript code/genetics/fastbat.manhattan.R "${trait}" "${targetDir}" "${fastbatResults}" "${annotationThresh}" "${ylim}" "${ysteps}" "${width}" "${height}"
	) &
done
wait)

# create gene-based qq-plots
pCol="Pvalue" # column containing p-values
nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5
prune=TRUE # ommit overplotting of variants with pval > 0.01
drawCI=FALSE # draw confidence interval?
drawLambda=FALSE # draw LambdaGC?
xend=8 # upper x axis limit
xsteps=2 # x axis breaks
yend=40 # upper y axis limit
ysteps=10 # y axis breaks
width=3 # plot width (4 inch)
height=4 # plot height (4 inch)
(for trait in gap_gm gap_wm gap_gwm; do (
	targetDir="results/${trait}/gwama/eur/fastbat/"
	sumstats="results/${trait}/gwama/eur/fastbat/fastbat.0kb.gene.fastbat.annovar.clumped"
	Rscript code/genetics/qqplot.R "${trait}" "${targetDir}" "${sumstats}" "${pCol}" "${nCriterion}" "${prune}" "${drawCI}" "${drawLambda}" "${xend}" "${xsteps}" "${yend}" "${ysteps}" "${width}" "${height}"
	) &
done
wait)

# combine Manhattan and qq-plots
traits="gap_gm,gap_wm,gap_gwm"
plotTitles="Grey_matter,White_matter,Grey_and_white_matter" # annotation will be a-f if no plot titles are provided
manhattanPlots="results/gap_gm/gwama/eur/fastbat/fastbat.manhattan.png,results/gap_wm/gwama/eur/fastbat/fastbat.manhattan.png,results/gap_gwm/gwama/eur/fastbat/fastbat.manhattan.png"
qqPlots="results/gap_gm/gwama/eur/fastbat/qqplot.png,results/gap_wm/gwama/eur/fastbat/qqplot.png,results/gap_gwm/gwama/eur/fastbat/qqplot.png"
outputFile="results/combined/gwama.eur.fastbat.qqplot.manhattan.png"
width=14
height=12
Rscript code/genetics/manhattan.qq.combine.R "${traits}" "${plotTitles}" "${manhattanPlots}" "${qqPlots}" "${outputFile}" "${width}" "${height}"

# =====================================================
# === Discovery + replication: run pathway analyses ===
# =====================================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage
conda activate envs/gofuncr

# run gene set enrichment test
fastbatFile="results/combined/gwama.eur.fastbat.txt"
fastbatGeneCol="Gene"
fastbatClumpCol="LOCUS_COUNT" # column that indicates locus of gene
crossFWER=TRUE # calculate FWER across the three ontologies (instead of ontology-wise calculation)
(for trait in gap_gm gap_wm gap_gwm; do (
	mkdir -p "results/${trait}/gwama/eur/gofuncr/"
	fastbatpCol="${trait}_Pvalue"
	outFile="results/${trait}/gwama/eur/gofuncr/gofuncr.gsea"
	Rscript code/genetics/gofuncr.gsea.R "${fastbatFile}" "${fastbatpCol}" "${fastbatGeneCol}" "${fastbatClumpCol}" "${crossFWER}" "${outFile}"
	) &
done
wait)

# combine results of gene set enrichment tests
traits="gap_gm,gap_wm,gap_gwm"
overFiles="results/gap_gm/gwama/eur/gofuncr/gofuncr.gsea.results.txt,results/gap_wm/gwama/eur/gofuncr/gofuncr.gsea.results.txt,results/gap_gwm/gwama/eur/gofuncr/gofuncr.gsea.results.txt"
outFile="results/combined/gwama.eur.gofuncr.gsea.txt"
Rscript code/genetics/gofuncr.gsea.combine.R "${traits}" "${overFiles}" "${outFile}"

# =================================================================
# === Discovery + replication: run genetic correlation analysis ===
# =================================================================

# set working directory and load conda environment
cd "/slow/projects/ukb_brainage"
conda activate envs/ldsc

# calculate rg with selected traits
threads=100 # set number of parallel workers
LDchr="data/ldsc/resources/eur_w_ld_chr/"
sumstatsY=$(find data/sumstats/munged/*sumstats.gz)
xprefix=""
yprefix=""
for trait in gap_gm gap_wm gap_gwm; do
	targetDir="results/${trait}/gwama/eur/rgSelection"
	sumstatsX="results/${trait}/gwama/eur/ldsc/ldsc.${trait}.sumstats.gz"
	./code/genetics/rg.sh "${targetDir}" "${sumstatsX}" "${sumstatsY}" "${xprefix}" "${yprefix}" "${LDchr}" "${threads}"
done

	# combine rg results across brain age gap variables
	conda activate envs/default
	traits="gap_gm,gap_wm,gap_gwm"
	inputFiles=$(echo $(for i in $(echo ${traits} | sed 's/,/ /g'); do echo results/${i}/gwama/eur/rgSelection/rg.results; done) | sed 's/ /,/g')
	outputFile="results/combined/gwama.eur.rgSelection.txt"
	./code/genetics/rg.combine.R "${traits}" "${inputFiles}" "${outputFile}"

	# replace trait names by trait labels
	inputFile="results/combined/gwama.eur.rgSelection.txt"
	outFile="results/combined/gwama.eur.rgSelection.labels.txt"
	matchVar="p2"
	./code/genetics/rg.addLabels.R "${inputFile}" "${outFile}" "${matchVar}" 

	# plot rg results
	corrTable="results/combined/gwama.eur.rgSelection.txt"
	traits="gap_gm,gap_wm,gap_gwm"
	traitLabels="Grey_matter,White_matter,Grey_&_White"
	matchVar="p2"
	outFile="results/combined/gwama.eur.rgSelection.png"
	scaleLimit=0.3
	width=2.76
	height=6.85
	./code/genetics/rg.plotSelection.R "${corrTable}" "${traits}" "${traitLabels}" "${matchVar}" "${outFile}" "${scaleLimit}" "${width}" "${height}"

# calculate rg with Neale traits
conda activate envs/ldsc
nealeDir="data/ldsc/resources/neale"
nealeManifest="data/ldsc/resources/neale_manifest_dropbox.tsv"
showcaseCategories="data/ldsc/resources/neale_showcase_categories.txt"
LDchr="data/ldsc/resources/eur_w_ld_chr/"
threads=100
for trait in gap_gm gap_wm gap_gwm; do
    targetDir="results/${trait}/gwama/eur/rg"
    sumstats="results/${trait}/gwama/eur/ldsc/ldsc.${trait}.sumstats.gz"
    ./code/genetics/rgNeale.sh "${trait}" "${targetDir}" "${sumstats}" "${nealeDir}" "${nealeManifest}" "${showcaseCategories}" "${LDchr}" "${threads}" # taskset -c 0-24 
done

	# combine rg results of Neale traits
	conda activate envs/default
	traits="gap_gm,gap_wm,gap_gwm"
	inputFiles=$(echo $(for i in $(echo ${traits} | sed 's/,/ /g'); do echo results/${i}/gwama/eur/rg/rg.results; done) | sed 's/ /,/g')
	outputFile="results/combined/gwama.eur.rgNeale.txt"
	./code/genetics/rgNeale.combine.R "${traits}" "${inputFiles}" "${outputFile}"

	# draw volcano and forest plot of rg Neale results
	inputFile="results/combined/gwama.eur.rgNeale.txt"
	outFile="results/combined/gwama.eur.rgNeale.volcano.forest.png"
	rgCol="gap_gm_rg"
	seCol="gap_gm_se"
	pCol="gap_gm_p"
	multipleTesting='fdr' # 'fdr' or 'bonferroni' or 'both'
	ylim=8 # upper y axis limit
	ysteps=2 # y axis breaks
	xlim=0.5 # upper y axis limit 
	xsteps=0.25 # y axis breaks 
	width=7.56 # plot width
	height=8.58 # 3.85 # plot height
	./code/genetics/rgNeale.volcano.forest.R "${inputFile}" "${outFile}" "${rgCol}" "${seCol}" "${pCol}" "${multipleTesting}" "${ylim}" "${ysteps}" "${xlim}" "${xsteps}" "${width}" "${height}"

	# combine plot of selected traits and Neale traits
	rgSelection="results/combined/gwama.eur.rgSelection.png"
	rgNeale="results/combined/gwama.eur.rgNeale.volcano.forest.png"
	outFile="results/combined/gwama.eur.rgCombined.png"
	./code/genetics/rg.plotCombine.R "${rgSelection}" "${rgNeale}" "${outFile}"

	# draw qq plot of rg Neale results
	inputFile="results/combined/gwama.eur.rgNeale.txt"
	outFile="results/combined/gwama.eur.rgNeale.qq.png"
	rgCol="gap_gm_rg"
	pCol="gap_gm_p"
	ylim=9 # upper y axis limit
	ysteps=3 # y axis breaks
	xlim=3 # upper y axis limit;
	xsteps=1 # y axis breaks
	width=4.79 # plot width
	height=4.57 # plot height
	./code/genetics/rgNeale.qq.R "${inputFile}" "${outFile}" "${rgCol}" "${pCol}" "${ylim}" "${ysteps}" "${xlim}" "${xsteps}" "${width}" "${height}"

# ============================================================
# === Discovery + replication: run Mendelian Randomization ===
# ============================================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage
conda activate envs/mendelian

# run gcta-gsmr
sumstatsCols="ID,A1,A2,A1_FREQ,BETA,SE,P,N"
mrFiles=$(echo $(ls data/sumstats/mr/mr_*.gz) | sed 's/ /,/g')
mrLabels=$(echo $(for i in $(echo $mrFiles | sed 's/,/ /g'); do echo $i | sed 's%.*/%%g' | sed 's%\..*%%g'; done) | sed 's/ /,/g')
chrFilehandler='data/genetics/chr${i}/imp_mri_qc_EURjoined/bed/chr${i}_mri_qc'
clumpr2=0.001
clumpkb=10000
gsmr2beta=0
(for trait in gap_gm gap_wm gap_gwm; do (
	outFile="results/${trait}/gwama/eur/gsmr/gsmr"
	sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz"
	sampleFile="results/mri/iid.EURjoined.txt"
	./code/genetics/gsmr.sh "${trait}" "${outFile}" "${sumstats}" "${sumstatsCols}" "${sampleFile}" "${mrFiles}" "${mrLabels}" "${chrFilehandler}" "${clumpr2}" "${clumpkb}" "${gsmr2beta}"
	) &
done
wait)

# run multiple additional MR analyses on variants pre-selected by gcta-gsmr (without HEIDI outlier removal)
(for trait in gap_gm gap_wm gap_gwm; do (
	gsmrFile="results/${trait}/gwama/eur/gsmr/gsmr.eff_plot.gz"
	gsmrFileNoHEIDI="results/${trait}/gwama/eur/gsmr/gsmr.noHEIDI.eff_plot.gz"
	outFile="results/${trait}/gwama/eur/gsmr/gsmr.multi"
	./code/genetics/gsmr.multi.R "${gsmrFile}" "${gsmrFileNoHEIDI}" "${outFile}"
	) &
done
wait)

# combine mr results
traits="gap_gm,gap_wm,gap_gwm"
inFiles="results/gap_gm/gwama/eur/gsmr/gsmr.multi,results/gap_wm/gwama/eur/gsmr/gsmr.multi,results/gap_gwm/gwama/eur/gsmr/gsmr.multi"
outFile="results/combined/gwama.eur.gsmr.multi.exposure.txt"
direction="exposure"
mInstruments=-1
./code/genetics/gsmr.multi.combine.R "${traits}" "${inFiles}" "${outFile}" "${direction}" "${mInstruments}"
outFile="results/combined/gwama.eur.gsmr.multi.outcome.txt"
direction="outcome"
mInstruments=-1
./code/genetics/gsmr.multi.combine.R "${traits}" "${inFiles}" "${outFile}" "${direction}" "${mInstruments}"

# replace trait names by trait labels
inputFile="results/combined/gwama.eur.gsmr.multi.exposure.txt"
outFile="results/combined/gwama.eur.gsmr.multi.exposure.labels.txt"
matchVar="trait"
./code/genetics/gsmr.addLabels.R "${inputFile}" "${outFile}" "${matchVar}" 
inputFile="results/combined/gwama.eur.gsmr.multi.outcome.txt"
outFile="results/combined/gwama.eur.gsmr.multi.outcome.labels.txt"
matchVar="trait"
./code/genetics/gsmr.addLabels.R "${inputFile}" "${outFile}" "${matchVar}" 

# plot mr results
traits=(gap_gm gap_wm gap_gwm)
outcome_labels=("Grey matter brain age gap" "White matter brain age gap" "Grey and white matter brain age gap")
for (( i=0; i<${#traits[@]}; i++ )); do
	gsmrFile="results/${traits[i]}/gwama/eur/gsmr/gsmr.eff_plot.gz"
	exposure_levels="mr_bmi_locke_2015,mr_whr_shungin_2015,mr_dbp_noUKB_evangelou_2018,mr_sbp_noUKB_evangelou_2018,mr_pp_noUKB_evangelou_2018,mr_ldlc_willer_2013,mr_hdlc_willer_2013,mr_trigl_willer_2013,mr_cad_nikpay_2015,mr_t2d_scott_2017,mr_scz_trubetskoy_2022,mr_eduPruned_okbay_2016"
	exposure_labels="Body-Mass-Index_(BMI),Waist-Hip-Ratio_(BMI adjusted),Diastolic Blood Pressure,Systolic Blood Pressure,Pulse Pressure,LDL-c,HDL-c,Triglyceride,Coronary_Artery_Disease,Type-2-Diabetes,Schizophrenia,Educational_Attainment"
	outcome_levels=${traits[i]}
	titles="Body-Mass-Index (Exposure),Waist-Hip-Ratio_(Exposure),Diastolic Blood Pressure (Exposure),Systolic Blood Pressure (Exposure),Pulse Pressure (Exposure),LDL-c (Exposure),HDL-c (Exposure),Triglyceride (Exposure),Coronary_Artery_Disease (Exposure),Type-2-Diabetes (Exposure),Schizophrenia (Exposure),Educational Attainment (Exposure)"
	outFile="results/${traits[i]}/gwama/eur/gsmr/gsmr.outcome"
	./code/genetics/gsmr.plot.R "${gsmrFile}" "${exposure_levels}" "${exposure_labels}" "${outcome_levels}" "${outcome_labels}" "${titles}" "${outFile}"
done

traits=(gap_gm gap_wm gap_gwm)
exposure_labels=("Grey matter brain age gap" "White matter brain age gap" "Grey and white matter brain age gap")
for (( i=0; i<${#traits[@]}; i++ )); do
	exposure_levels=${traits[i]}
	exposure_labels=${exposure_labels[i]}
	gsmrFile="results/${traits[i]}/gwama/eur/gsmr/gsmr.eff_plot.gz"
	outcome_levels="mr_bmi_locke_2015,mr_whr_shungin_2015,mr_dbp_noUKB_evangelou_2018,mr_sbp_noUKB_evangelou_2018,mr_pp_noUKB_evangelou_2018,mr_ldlc_willer_2013,mr_hdlc_willer_2013,mr_trigl_willer_2013,mr_cad_nikpay_2015,mr_t2d_scott_2017,mr_scz_trubetskoy_2022"
	outcome_labels="Body-Mass-Index_(BMI),Waist-Hip-Ratio_(BMI adjusted),Diastolic Blood Pressure,Systolic Blood Pressure,Pulse Pressure,LDL-c,HDL-c,Triglyceride,Coronary_Artery_Disease,Type-2-Diabetes,Schizophrenia"
	titles="Body-Mass-Index (Outcome),Waist-Hip-Ratio (Outcome),Diastolic Blood Pressure (Outcome),Systolic Blood Pressure (Outcome),Pulse Pressure (Outcome),LDL-c (Outcome),HDL-c (Outcome),Triglyceride (Outcome),Coronary_Artery_Disease (Outcome),Type-2-Diabetes (Outcome),Schizophrenia (Outcome)"
	outFile="results/${traits[i]}/gwama/eur/gsmr/gsmr.exposure"
	./code/genetics/gsmr.plot.R "${gsmrFile}" "${exposure_levels}" "${exposure_labels}" "${outcome_levels}" "${outcome_labels}" "${titles}" "${outFile}"
done

# ========================================================================================
# === Discovery + replication: run genetic effect size distribution analysis (GENESIS) ===
# ========================================================================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage
conda activate envs/genesis

# run analysis for traits of interest
colsIn="ID,BETA,SE,N"
colsOut="SNP,BETA,SE,N"
ncores=20
compFit3=TRUE
(for trait in gap_gm gap_wm gap_gwm; do (
	mkdir -p results/${trait}/genesis/
	sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz"
	outFile="results/${trait}/gwama/eur/genesis/genesis"
	taskset -c 0-59 Rscript code/genetics/genesis.R "${trait}" "${sumstats}" "${colsIn}" "${colsOut}" "${ncores}" "${compFit3}" "${outFile}"
	) &
done
wait)

# run analysis for reference traits
trait="04_neur_baselmans_2019"
mkdir -p data/sumstats/genesis/${trait}/
sumstats="data/sumstats/harmonized/${trait}.gz"
colsIn="ID,Z,N"
colsOut="SNP,Z,N"
ncores=10
compFit3=TRUE
outFile="data/sumstats/genesis/${trait}/genesis"
Rscript code/genetics/genesis.R "${trait}" "${sumstats}" "${colsIn}" "${colsOut}" "${ncores}" "${compFit3}" "${outFile}"

trait="mr_height_wood_2014"
mkdir -p data/sumstats/genesis/${trait}/
sumstats="data/sumstats/harmonized/${trait}.gz"
colsIn="ID,BETA,SE,N"
colsOut="SNP,BETA,SE,N"
ncores=10
compFit3=TRUE
outFile="data/sumstats/genesis/${trait}/genesis"
Rscript code/genetics/genesis.R "${trait}" "${sumstats}" "${colsIn}" "${colsOut}" "${ncores}" "${compFit3}" "${outFile}"

# plot results
traitLabels="Brain age gap"
reftraits="neur,height"
reftraitLabels="Neuroticism,Height"
refgenesisFit="data/sumstats/genesis/04_neur_baselmans_2019/genesis.fit3.Rds,data/sumstats/genesis/mr_height_wood_2014/genesis.fit3.Rds"
for trait in gap_gm gap_wm gap_gwm; do
	genesisFit="results/${trait}/gwama/eur/genesis/genesis.fit3.Rds"
	outFile="results/${trait}/gwama/eur/genesis/genesis"
	Rscript code/genetics/genesis.plot.R "${trait}" "${traitLabels}" "${genesisFit}" "${reftraits}" "${reftraitLabels}" "${refgenesisFit}" "${outFile}"
done

# create stats table
traits="gap_gm,gap_wm,gap_gwm,height,neur"
traitLabels="Grey_matter_brain_age_gap,White_matter_brain_age_gap,Grey_and_white_matter_brain_age_gap,Height_(Wood_et_al._2014),Neuroticism_(Baselmans_et_al._2019)"
genesisFit="results/gap_gm/gwama/eur/genesis/genesis.fit3.Rds,results/gap_wm/gwama/eur/genesis/genesis.fit3.Rds,results/gap_gwm/gwama/eur/genesis/genesis.fit3.Rds,data/sumstats/genesis/mr_height_wood_2014/genesis.fit3.Rds,data/sumstats/genesis/04_neur_baselmans_2019/genesis.fit3.Rds"
outFile="results/combined/gwama.eur.genesis.stats"
Rscript code/genetics/genesis.stats.R "${traits}" "${traitLabels}" "${genesisFit}" "${outFile}"

# ====================================
# === Run polygenic score analysis ===
# ====================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage
conda activate envs/default

# Discovery Sample: Get SbayesR + SBayesRC weights

	# run SBayesRC with eigen-decomposition data and baseline v2.2 annotations obtained from GCTB webpage
	cols="ID,A1,A2,A1_FREQ,BETA,SE,P,N" # set input column names that correspond to header columns SNP A1 A2 freq b se p N
	ldmFolder="/fast/software/gctb/resources/ukbEUR_Imputed"
	annotFile="/fast/software/gctb/resources/annot_baseline2.2.txt"
	threads=15
	sbayes="sbayesrc"
	imputation=1
	trait="gap_gm"; targetDir="results/${trait}/discovery/sbayes"; sumstats="results/${trait}/discovery/gwas/sumstats.txt.gz";
		taskset -c 0-14 ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}"
	trait="gap_wm"; targetDir="results/${trait}/discovery/sbayes"; sumstats="results/${trait}/discovery/gwas/sumstats.txt.gz";
		taskset -c 15-29 ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}"
	trait="gap_gwm"; targetDir="results/${trait}/discovery/sbayes"; sumstats="results/${trait}/discovery/gwas/sumstats.txt.gz";
		taskset -c 30-45 ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}"

	# run SBayesR with eigen-decomposition data and baseline v2.2 annotations obtained from GCTB webpage
	cols="SNP,A1,A2,freq,b,se,p,N" # set input column names that correspond to header columns SNP A1 A2 freq b se p N
	ldmFolder="/fast/software/gctb/resources/ukbEUR_Imputed"
	annotFile="/fast/software/gctb/resources/annot_baseline2.2.txt"
	threads=15
	sbayes="sbayesr"
	imputation=0
	trait="gap_gm"; targetDir="results/${trait}/discovery/sbayes"; sumstats="results/${trait}/discovery/sbayes/sbayesrc.imputed.ma.gz";
		taskset -c 0-14  ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}"
	trait="gap_wm"; targetDir="results/${trait}/discovery/sbayes"; sumstats="results/${trait}/discovery/sbayes/sbayesrc.imputed.ma.gz";
		taskset -c 15-29 ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}"
	trait="gap_gwm"; targetDir="results/${trait}/discovery/sbayes"; sumstats="results/${trait}/discovery/sbayes/sbayesrc.imputed.ma.gz";
		taskset -c 30-45 ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}"

	# run SBayesRC with custom eigen-decomposition data (9M variants from 32k disovery individuals) and baselineLF.v2.2.UKB annotations
	cols="ID,A1,A2,A1_FREQ,BETA,SE,P,N" # set input column names that correspond to header columns SNP A1 A2 freq b se p N
	ldmFolder="/fast/software/gctb/resources/ukbEUR_32k/ldm"
	annotFile="/fast/software/gctb/resources/ukbEUR_32k/ldm/annot_baselineLF_v2.2.UKB_32k.txt"
	threads=50
	sbayes="sbayesrc"
	imputation=1
	trait="gap_gm"; targetDir="results/${trait}/discovery/sbayes_ld32k"; sumstats="results/${trait}/discovery/gwas/sumstats.txt.gz"
		taskset -c 0-49 ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}"
	trait="gap_wm"; targetDir="results/${trait}/discovery/sbayes_ld32k"; sumstats="results/${trait}/discovery/gwas/sumstats.txt.gz"
		taskset -c 0-49 ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}"
	trait="gap_gwm"; targetDir="results/${trait}/discovery/sbayes_ld32k"; sumstats="results/${trait}/discovery/gwas/sumstats.txt.gz"
		taskset -c 56-105 ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}"
	for trait in gap_gm gap_wm gap_gwm; do
		mv "results/${trait}/discovery/sbayes_ld32k/sbayesrc.badSNPlist" "results/${trait}/discovery/sbayes/sbayesrc_ld32k.badSNPlist"
		mv "results/${trait}/discovery/sbayes_ld32k/sbayesrc.imputed.ma.gz" "results/${trait}/discovery/sbayes/sbayesrc_ld32k.imputed.ma.gz"
		mv "results/${trait}/discovery/sbayes_ld32k/sbayesrc.parRes" "results/${trait}/discovery/sbayes/sbayesrc_ld32k.parRes"
		mv "results/${trait}/discovery/sbayes_ld32k/sbayesrc.parSetRes" "results/${trait}/discovery/sbayes/sbayesrc_ld32k.parSetRes"
		mv "results/${trait}/discovery/sbayes_ld32k/sbayesrc.snpRes.gz" "results/${trait}/discovery/sbayes/sbayesrc_ld32k.snpRes.gz"
	done

		# calculate PGS in replication sample based on discovery SBayesRC weights
		pgenFileHandlerBeta='data/genetics/chr${i}/imp_mri_qc_ANCESTRY/chr${i}_mri_qc'
		idCol="Name"
		a1Col="A1"
		effectCol="A1Effect"
		maf=0.01
		threads=100
		for ancestry in EUR; do
			for trait in gap_gm gap_wm gap_gwm; do
				for pgsMethod in sbayesrc_ld32k sbayesrc sbayesr; do
					pgenFileHandler=$(echo ${pgenFileHandlerBeta} | sed "s/ANCESTRY/${ancestry}/g")
					outFile="data/genetics/prs/imp_mri_qc_${ancestry}/${trait}.${pgsMethod}.32k"
					weightFile="results/${trait}/discovery/sbayes/${pgsMethod}.snpRes.gz"
					code/genetics/pgs.predict.sh "${weightFile}" "${pgenFileHandler}" "${idCol}" "${a1Col}" "${effectCol}" "${maf}" "${threads}" "${outFile}"
				done
			done
		done

		# correlate SBayesRC PGS with phenotype in replication sample
		phenoFile="data/traits/replicate.txt"
		joinVar="IID"
		pgsVar="SCORE1_SUM"
		covs="sex,age,age2,ac1,ac2,ac3,TIV,array,PanC1,PanC2,PanC3,PanC4,PanC5,PanC6,PanC7,PanC8,PanC9,PanC10,PanC11,PanC12,PanC13,PanC14,PanC15,PanC16,PanC17,PanC18,PanC19,PanC20"
		for ancestry in EUR; do # AFR AMR CSA EAS EUR MID
			for trait in gap_gm gap_wm gap_gwm; do
				for pgsMethod in sbayesrc_ld32k sbayesrc sbayesr; do
					phenoVar=${trait}
					pgsFile="data/genetics/prs/imp_mri_qc_${ancestry}/${trait}.${pgsMethod}.32k.score"
					outFile="≈${pgsMethod}.base32k.target20k.assoc.txt"
					code/genetics/pgs.corr.R "${pgsMethod}" "${pgsFile}" "${phenoFile}" "${joinVar}" "${pgsVar}" "${phenoVar}" "${covs}" "${outFile}"
				done
			done
		done

# use PRSice-2 to calculate PGS in replication sample
for ancestry in EUR; do
	mkdir -p data/genetics/prs/imp_mri_qc_${ancestry}/
	awk 'NR<3 { print; next } { print $1"_"$2, $1"_"$2, $3, $4 }' data/genetics/chr1/imp_mri_qc_${ancestry}/bgen/chr1_mri_qc.sample > data/genetics/prs/imp_mri_qc_${ancestry}/chr1_mri_qc.sample
	(for trait in gap_gm gap_wm gap_gwm; do (
		Rscript /fast/software/PRSice2/PRSice.R \
			--dir /fast/software/PRSice2/ \
			--prsice /fast/software/PRSice2/PRSice_linux \
			--no-regress T \
			--base "results/${trait}/discovery/gwas/sumstats.txt.gz" \
			--snp ID \
			--base-info INFO:0.8 \
			--target "data/genetics/chr#/imp_mri_qc_${ancestry}/bgen/chr#_mri_qc","data/genetics/prs/imp_mri_qc_${ancestry}/chr1_mri_qc.sample" \
			--maf 0.01 \
			--type bgen \
			--ld-type bgen \
			--allow-inter T \
			--clump-p 1 \
			--clump-r2 0.1 \
			--clump-kb 500 \
			--score std \
			--bar-levels 1.00,0.50,0.20,0.10,0.05,0.01,1E-3,1E-4,1E-6,5E-8 \
			--fastscore \
			--out "data/genetics/prs/imp_mri_qc_${ancestry}/${trait}" \
			--thread 25
			) &
	done)
	wait
done

# Randomly select 2,000 EUR-ancestry replication individuals.
# Re-run the replication GWAS excluding these 2,000 individuals
# Re-run the GWAS meta-analysis across discovery and replication exclduding the same 2,000 individuals
# Re-run SBayesR/RC
# Re-run PRSice-2
# Evaluate polygenic scores (PGS) accuracy in the selected 2,000 individuals using weights from discovery-only and discovery+replication datasets.

	# 1) Randomly select 2,000 individuals
	rows=$(awk -F'\t' '$3=="EUR" {print NR}' data/traits/replicate.txt | shuf --random-source=<(yes 123) | head -2000)
	awk 'NR==FNR { selected[$1]; next } !(FNR in selected) { print }' <(echo "${rows}") data/traits/replicate.txt > data/traits/replicate.excl.2k.txt
	awk 'NR==FNR { selected[$1]; next } FNR==1 { print; next } FNR in selected { print }' <(echo "${rows}") data/traits/replicate.txt > data/traits/replicate.2k.txt

	# 2) Re-run the replication GWAS excluding these 2,000 individuals
	chrFileHandle='data/genetics/chr${i}/imp_mri_qc_${ancestry}/chr${i}_mri_qc'
	data="data/traits/replicate.excl.2k.txt"
	ancestriesCol="pan"
	snpQC="data/genetics/meta/ukb_snp_qc.txt"
	impMFI="data/genetics/ukb_imp_mfi"
	maf=0.01
	threads=100
	ancestries="EUR"
	covs="sex,age,age2,ac1,ac2,ac3,TIV,array,PanC{1..20}"
	for trait in gap_gm gap_wm gap_gwm; do
		#targetDir="results/${trait}/replicate/EUR/EUR.excl.2k"
		#taskset -c 0-99 code/genetics/replicate.sh "${trait}" "${targetDir}" "${chrFileHandle}" "${data}" "${covs}" "${ancestries}" "${ancestriesCol}" "${snpQC}" "${impMFI}" "${maf}" "${threads}"
		mv "results/${trait}/replicate/EUR/EUR.excl.2k/EUR/"* "results/${trait}/replicate/EUR/EUR.excl.2k/"
	done

	# 3a) Re-run the GWAS meta-analysis across discovery and replication excluding the same 2,000 individuals
	cohorts="DISCOV,EUR_excl_2k,LIFE" # set output file handler & make sure that the order matches the list of sumstats (below)
	type="ivweight" # type="nweight"
	idCol="ID"
	carryOver="CHR,BP"
	leaveOneOut=0
	(for trait in gap_gm gap_wm gap_gwm; do (
		targetDir="results/${trait}/gwama/eur/eur.excl.2k/"; mkdir -p ${targetDir}
		sumstats="results/${trait}/discovery/gwas/sumstats.txt.gz,results/${trait}/replicate/EUR/EUR.excl.2k/sumstats.txt.gz,results/${trait}/replicate/LIFE/sumstats.txt.gz"
		code/genetics/metal.sh ${trait} ${targetDir} ${cohorts} ${sumstats} ${type} ${idCol} ${carryOver} ${leaveOneOut}
		) &
	done
	wait)

	# 3b)quality-check METAL results
	nCriterion=TRUE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5
	hetP=1E-6 # heterogeneity p threshold for excluding SNPs
	(for trait in gap_gm gap_wm gap_gwm; do (
		sumstats="results/${trait}/gwama/eur/eur.excl.2k/metal.ivweight.gz"
		outFile="results/${trait}/gwama/eur/eur.excl.2k/metal.ivweight.qc"
		Rscript code/genetics/metal.qc.R ${sumstats} ${outFile} ${nCriterion} ${hetP}
		) &
	done
	wait)

	# 4a) run SBayesRC with eigen-decomposition data and baseline v2.2 annotations obtained from GCTB webpage
	cols="ID,A1,A2,A1_FREQ,BETA,SE,P,N" # set input column names that correspond to header columns SNP A1 A2 freq b se p N
	ldmFolder="/fast/software/gctb/resources/ukbEUR_Imputed"
	annotFile="/fast/software/gctb/resources/annot_baseline2.2.txt"
	threads=18
	sbayes="sbayesrc"
	imputation=1
	trait="gap_gm"; targetDir="results/${trait}/gwama/eur/eur.excl.2k"; sumstats="results/${trait}/gwama/eur/eur.excl.2k/metal.ivweight.qc.gz";
		taskset -c 0-17 ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}"
	trait="gap_wm"; targetDir="results/${trait}/gwama/eur/eur.excl.2k"; sumstats="results/${trait}/gwama/eur/eur.excl.2k/metal.ivweight.qc.gz";
		taskset -c 18-35 ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}"
	trait="gap_gwm"; targetDir="results/${trait}/gwama/eur/eur.excl.2k"; sumstats="results/${trait}/gwama/eur/eur.excl.2k/metal.ivweight.qc.gz";
		taskset -c 36-53 ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}"

	# 4b) run SBayesR with eigen-decomposition data and baseline v2.2 annotations obtained from GCTB webpage
	cols="SNP,A1,A2,freq,b,se,p,N" # set input column names that correspond to header columns SNP A1 A2 freq b se p N
	sbayes="sbayesr"
	imputation=0
	trait="gap_gm"; targetDir="results/${trait}/gwama/eur/eur.excl.2k"; sumstats="results/${trait}/gwama/eur/eur.excl.2k/sbayesrc.imputed.ma.gz";
		taskset -c 0-17 ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}"
	trait="gap_wm"; targetDir="results/${trait}/gwama/eur/eur.excl.2k"; sumstats="results/${trait}/gwama/eur/eur.excl.2k/sbayesrc.imputed.ma.gz";
		taskset -c 18-35 ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}"
	trait="gap_gwm"; targetDir="results/${trait}/gwama/eur/eur.excl.2k"; sumstats="results/${trait}/gwama/eur/eur.excl.2k/sbayesrc.imputed.ma.gz";
		taskset -c 36-53 ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}"

	# 4c) run SBayesRC with custom eigen-decomposition data (9M variants from 32k disovery individuals) and baselineLF.v2.2.UKB annotations
	cols="ID,A1,A2,A1_FREQ,BETA,SE,P,N" # set input column names that correspond to header columns SNP A1 A2 freq b se p N
	ldmFolder="/fast/software/gctb/resources/ukbEUR_32k/ldm"
	annotFile="/fast/software/gctb/resources/ukbEUR_32k/ldm/annot_baselineLF_v2.2.UKB_32k.txt"
	threads=30
	sbayes="sbayesrc"
	imputation=1
	trait="gap_gm"; targetDir="results/${trait}/gwama/eur/eur.excl.2k/sbayes_ld32k"; sumstats="results/${trait}/gwama/eur/eur.excl.2k/metal.ivweight.qc.gz";
		taskset -c 0-29 ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}"
	trait="gap_wm"; targetDir="results/${trait}/gwama/eur/eur.excl.2k/sbayes_ld32k"; sumstats="results/${trait}/gwama/eur/eur.excl.2k/metal.ivweight.qc.gz";
		taskset -c 40-79 ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}"
	trait="gap_gwm"; targetDir="results/${trait}/gwama/eur/eur.excl.2k/sbayes_ld32k"; sumstats="results/${trait}/gwama/eur/eur.excl.2k/metal.ivweight.qc.gz";
		taskset -c 80-109 ./code/genetics/sbayes.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}"
	for trait in gap_gm gap_wm gap_gwm; do
		mv "results/${trait}/gwama/eur/eur.excl.2k/sbayes_ld32k/sbayesrc.badSNPlist" "results/${trait}/gwama/eur/eur.excl.2k/sbayesrc_ld32k.badSNPlist"
		mv "results/${trait}/gwama/eur/eur.excl.2k/sbayes_ld32k/sbayesrc.imputed.ma.gz" "results/${trait}/gwama/eur/eur.excl.2k/sbayesrc_ld32k.imputed.ma.gz"
		mv "results/${trait}/gwama/eur/eur.excl.2k/sbayes_ld32k/sbayesrc.parRes" "results/${trait}/gwama/eur/eur.excl.2k/sbayesrc_ld32k.parRes"
		mv "results/${trait}/gwama/eur/eur.excl.2k/sbayes_ld32k/sbayesrc.parSetRes" "results/${trait}/gwama/eur/eur.excl.2k/sbayesrc_ld32k.parSetRes"
		mv "results/${trait}/gwama/eur/eur.excl.2k/sbayes_ld32k/sbayesrc.snpRes.gz" "results/${trait}/gwama/eur/eur.excl.2k/sbayesrc_ld32k.snpRes.gz"
	done
	rm -rf results/{gap_gm,gap_wm,gap_gwm}/gwama/eur/eur.excl.2k/sbayes_ld32k/

	# 5) Calculate polygenic scores (PGS) in the selected 2,000 individuals using weights from discovery-only and discovery+replication datasets.
	pgenFileHandlerBeta='data/genetics/chr${i}/imp_mri_qc_ANCESTRY/chr${i}_mri_qc'
	idCol="Name"
	a1Col="A1"
	effectCol="A1Effect"
	maf=0.01
	threads=100
	for ancestry in AFR CSA EAS EUR; do # do not use AMR and MID due to low sample sizes (n < 100) 
		for trait in gap_gm gap_wm gap_gwm; do
			for pgsMethod in sbayesrc_ld32k sbayesrc sbayesr; do
				pgenFileHandler=$(echo ${pgenFileHandlerBeta} | sed "s/ANCESTRY/${ancestry}/g")
				outFile="data/genetics/prs/imp_mri_qc_${ancestry}/${trait}.${pgsMethod}.52k"
				weightFile="results/${trait}/gwama/eur/eur.excl.2k/${pgsMethod}.snpRes.gz"
				taskset -c 0-99 code/genetics/pgs.predict.sh "${weightFile}" "${pgenFileHandler}" "${idCol}" "${a1Col}" "${effectCol}" "${maf}" "${threads}" "${outFile}"
			done
		done
	done

	# 6a) correlate SBayesR/RC PGS with phenotype in replication sample - European ancestry individuals
	phenoFile="data/traits/replicate.2k.txt"
	joinVar="IID"
	pgsVar="SCORE1_SUM"
	covs="sex,age,age2,ac1,ac2,ac3,TIV,array,PanC1,PanC2,PanC3,PanC4,PanC5,PanC6,PanC7,PanC8,PanC9,PanC10,PanC11,PanC12,PanC13,PanC14,PanC15,PanC16,PanC17,PanC18,PanC19,PanC20"
	for ancestry in EUR; do
		for trait in gap_gm gap_wm gap_gwm; do
			for pgsMethod in sbayesrc_ld32k sbayesrc sbayesr; do
				phenoVar=${trait}
				pgsFile="data/genetics/prs/imp_mri_qc_${ancestry}/${trait}.${pgsMethod}.32k.score"
				outFile="results/${trait}/replicate/${ancestry}/${pgsMethod}.base32k.target2k.assoc.txt"
				code/genetics/pgs.corr.R "${pgsMethod}" "${pgsFile}" "${phenoFile}" "${joinVar}" "${pgsVar}" "${phenoVar}" "${covs}" "${outFile}"
				pgsFile="data/genetics/prs/imp_mri_qc_${ancestry}/${trait}.${pgsMethod}.52k.score"
				outFile="results/${trait}/replicate/${ancestry}/${pgsMethod}.base52k.target2k.assoc.txt"
				code/genetics/pgs.corr.R "${pgsMethod}" "${pgsFile}" "${phenoFile}" "${joinVar}" "${pgsVar}" "${phenoVar}" "${covs}" "${outFile}"
			done	
		done
	done

	# 6b) correlate SBayesR/RC PGS with phenotype in replication sample - non-European ancestry individuals
	phenoFile="data/traits/replicate.txt"
	joinVar="IID"
	pgsVar="SCORE1_SUM"
	covs="sex,age,age2,ac1,ac2,ac3,TIV,array,PanC1,PanC2,PanC3,PanC4"
	for ancestry in AFR CSA EAS; do
		for trait in gap_gm gap_wm gap_gwm; do
			for pgsMethod in sbayesrc_ld32k sbayesrc sbayesr; do
				phenoVar=${trait}
				pgsFile="data/genetics/prs/imp_mri_qc_${ancestry}/${trait}.${pgsMethod}.52k.score"
				outFile="results/${trait}/replicate/${ancestry}/${pgsMethod}.base52k.assoc.txt"
				code/genetics/pgs.corr.R "${pgsMethod}" "${pgsFile}" "${phenoFile}" "${joinVar}" "${pgsVar}" "${phenoVar}" "${covs}" "${outFile}"
			done
		done
	done

	# 7) use PRSice-2 to calculate PGS in left-out 2k individuals
	for ancestry in AFR CSA EAS EUR; do
		mkdir -p data/genetics/prs/imp_mri_qc_${ancestry}/
		awk 'NR<3 { print; next } { print $1"_"$2, $1"_"$2, $3, $4 }' data/genetics/chr1/imp_mri_qc_${ancestry}/bgen/chr1_mri_qc.sample > data/genetics/prs/imp_mri_qc_${ancestry}/chr1_mri_qc.sample
		(for trait in gap_gm gap_wm gap_gwm; do (
			Rscript /fast/software/PRSice2/PRSice.R \
				--dir /fast/software/PRSice2/ \
				--prsice /fast/software/PRSice2/PRSice_linux \
				--no-regress T \
				--base "results/${trait}/gwama/eur/eur.excl.2k/metal.ivweight.qc.gz" \
				--snp ID \
				--target "data/genetics/chr#/imp_mri_qc_${ancestry}/bgen/chr#_mri_qc","data/genetics/prs/imp_mri_qc_${ancestry}/chr1_mri_qc.sample" \
				--maf 0.01 \
				--type bgen \
				--ld-type bgen \
				--allow-inter T \
				--clump-p 1 \
				--clump-r2 0.1 \
				--clump-kb 500 \
				--score std \
				--bar-levels 1.00,0.50,0.20,0.10,0.05,0.01,1E-3,1E-4,1E-6,5E-8 \
				--fastscore \
				--out "data/genetics/prs/imp_mri_qc_${ancestry}/${trait}.prsice.52k" \
				--thread 25
				) &
		done)
		wait
	done

	# 8) correlate PRSice-2 PGS with phenotype in replication sample
	covs="sex,age,age2,ac1,ac2,ac3,TIV,array,PanC1,PanC2,PanC3,PanC4,PanC5,PanC6,PanC7,PanC8,PanC9,PanC10,PanC11,PanC12,PanC13,PanC14,PanC15,PanC16,PanC17,PanC18,PanC19,PanC20"
	ancestry=EUR
	for trait in gap_gm gap_wm gap_gwm; do
		data="data/traits/replicate.txt"
			prsFile="data/genetics/prs/imp_mri_qc_${ancestry}/${trait}.all_score"
			outFile="results/${trait}/replicate/${ancestry}/prsice.base32k.target20k.assoc.txt"
			code/genetics/prsice.corr.R "${trait}" "${data}" "${prsFile}" "${covs}" "${outFile}"
		data="data/traits/replicate.2k.txt"
			prsFile="data/genetics/prs/imp_mri_qc_${ancestry}/${trait}.all_score"
			outFile="results/${trait}/replicate/${ancestry}/prsice.base32k.target2k.assoc.txt"
			code/genetics/prsice.corr.R "${trait}" "${data}" "${prsFile}" "${covs}" "${outFile}"
			prsFile="data/genetics/prs/imp_mri_qc_${ancestry}/${trait}.prsice.52k.all_score"
			outFile="results/${trait}/replicate/${ancestry}/prsice.base52k.target2k.assoc.txt"
			code/genetics/prsice.corr.R "${trait}" "${data}" "${prsFile}" "${covs}" "${outFile}"
	done
	covs="sex,age,age2,ac1,ac2,ac3,TIV,array,PanC1,PanC2,PanC3,PanC4"
	for ancestry in AFR CSA EAS; do
		for trait in gap_gm gap_wm gap_gwm; do
			data="data/traits/replicate.txt"
			prsFile="data/genetics/prs/imp_mri_qc_${ancestry}/${trait}.prsice.52k.all_score"
			outFile="results/${trait}/replicate/${ancestry}/prsice.base52k.assoc.txt"
			code/genetics/prsice.corr.R "${trait}" "${data}" "${prsFile}" "${covs}" "${outFile}"
		done
	done

	# 9) join PRS association results
		# join results from sbayesrc and prsice for each trait
		ancestry=EUR
		for trait in gap_gm gap_wm gap_gwm; do
			awk 'NR==1 { print; next } FNR==1 { next } { print }' "results/${trait}/replicate/${ancestry}/"{sbayesrc_ld32k,sbayesrc,sbayesr,prsice}".base32k.target20k.assoc.txt"  > "results/${trait}/replicate/${ancestry}/pgs.base32k.target20k.assoc.txt"
			awk 'NR==1 { print; next } FNR==1 { next } { print }' "results/${trait}/replicate/${ancestry}/"{sbayesrc_ld32k,sbayesrc,sbayesr,prsice}".base32k.target2k.assoc.txt"  > "results/${trait}/replicate/${ancestry}/pgs.base32k.target2k.assoc.txt"
			awk 'NR==1 { print; next } FNR==1 { next } { print }' "results/${trait}/replicate/${ancestry}/"{sbayesrc_ld32k,sbayesrc,sbayesr,prsice}".base52k.target2k.assoc.txt" > "results/${trait}/replicate/${ancestry}/pgs.base52k.target2k.assoc.txt"
		done
		for ancestry in AFR CSA EAS; do
			for trait in gap_gm gap_wm gap_gwm; do
				awk 'NR==1 { print; next } FNR==1 { next } { print }' "results/${trait}/replicate/${ancestry}/"{sbayesrc_ld32k,sbayesrc,sbayesr,prsice}".base52k.assoc.txt"  > "results/${trait}/replicate/${ancestry}/pgs.base52k.assoc.txt"
			done
		done

		# join results across traits
		traits="gap_gm,gap_wm,gap_gwm"
		pgsFiles="results/gap_gm/replicate/EUR/pgs.base32k.target20k.assoc.txt,results/gap_wm/replicate/EUR/pgs.base32k.target20k.assoc.txt,results/gap_gwm/replicate/EUR/pgs.base32k.target20k.assoc.txt"
		outFile="results/combined/pgs.base32k.target20k.assoc.txt"
			code/genetics/pgs.combine.R "${traits}" "${pgsFiles}" "${outFile}"
		pgsFiles="results/gap_gm/replicate/EUR/pgs.base32k.target2k.assoc.txt,results/gap_wm/replicate/EUR/pgs.base32k.target2k.assoc.txt,results/gap_gwm/replicate/EUR/pgs.base32k.target2k.assoc.txt"
		outFile="results/combined/pgs.base32k.target2k.assoc.txt"
			code/genetics/pgs.combine.R "${traits}" "${pgsFiles}" "${outFile}"
		pgsFiles="results/gap_gm/replicate/EUR/pgs.base52k.target2k.assoc.txt,results/gap_wm/replicate/EUR/pgs.base52k.target2k.assoc.txt,results/gap_gwm/replicate/EUR/pgs.base52k.target2k.assoc.txt"
		outFile="results/combined/pgs.base52k.target2k.assoc.txt"
			code/genetics/pgs.combine.R "${traits}" "${pgsFiles}" "${outFile}"
		for ancestry in AFR CSA EAS; do
			pgsFiles="results/gap_gm/replicate/${ancestry}/pgs.base52k.assoc.txt,results/gap_wm/replicate/${ancestry}/pgs.base52k.assoc.txt,results/gap_gwm/replicate/${ancestry}/pgs.base52k.assoc.txt"
			outFile="results/combined/pgs.base52k.${ancestry}.assoc.txt"
			code/genetics/pgs.combine.R "${traits}" "${pgsFiles}" "${outFile}"
		done

		# join results from different base / target cohorts
		awk 'NR==1 { print; next } FNR==1 { print ""; next} { print }' results/combined/pgs.base32k.target20k.assoc.txt results/combined/pgs.base{32,52}k.target2k.assoc.txt results/combined/pgs.base52k.{AFR,CSA,EAS}.assoc.txt > results/combined/pgs.assoc.txt
		rm -f results/combined/pgs.base52k* results/combined/pgs.base32k*


# combined discovery and replication dataset: in addition to SBayesRC weights from custom eigen-decomposition data (9M variants from 32k individuals, see fine-mapping section), get SBayesRC weights with eigen-decomposition data and baseline v2.2 annotations obtained from GCTB webpage
cols="ID,A1,A2,A1_FREQ,BETA,SE,P,N" # set input column names that correspond to header columns SNP A1 A2 freq b se p N
ldmFolder="/fast/software/gctb/resources/ukbEUR_Imputed"
annotFile="/fast/software/gctb/resources/annot_baseline2.2.txt"
threads=18
mcmc=
trait="gap_gm"; targetDir="results/${trait}/gwama/eur/sbayesrc"; sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz";
	taskset -c 0-17 ./code/genetics/sbayesrc.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}"
trait="gap_wm"; targetDir="results/${trait}/gwama/eur/sbayesrc"; sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz";
	taskset -c 18-35 ./code/genetics/sbayesrc.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}"
trait="gap_gwm"; targetDir="results/${trait}/gwama/eur/sbayesrc"; sumstats="results/${trait}/gwama/eur/metal.ivweight.qc.gz";
	taskset -c 36-53 ./code/genetics/sbayesrc.sh "${targetDir}" "${sumstats}" "${cols}" "${ldmFolder}" "${annotFile}" "${threads}"

# =================================================================
# === Explorative: compare genetics vs. phenotypic correlations ===
# =================================================================

# set working directory
cd /slow/projects/ukb_brainage
mamba activate envs/default

# make plots
traits="gap_gm,gap_wm,gap_gwm"
rgFile="results/combined/gwama.eur.rgNeale.txt"
phenoFile="results/combined/discovery.phewas.txt"
outFile="results/combined/rgVSrp.png"
width=7.77
height=2.83
ncols=3
./code/genetics/rgVSrp.R "${traits}" "${rgFile}" "${phenoFile}" "${outFile}" "${width}" "${height}" "${ncols}"

# =============================================================================
# === Explorative: genetic correlations between sex-stratified GWAS results ===
# =============================================================================

# set working directory
cd /slow/projects/ukb_brainage
mamba activate envs/default

# discovery gwas in sex-stratified samples
covsFile="data/traits/covs"
chrFileHandle="data/genetics/chr\${i}/imp_mri_qc/chr\${i}_mri_qc"
snpQC="data/genetics/meta/ukb_snp_qc.txt"
impMFI="data/genetics/ukb_imp_mfi"
maf=0.01
threads=100
awk 'NR==1 || $3==1 { out=$1; for(i=2;i<=NF;i++) { if(i!=3) { out=out"\t"$i} }; print out; out=""}' "${covsFile}.txt" > ${covsFile}.female.txt
awk 'NR==1 || $3==2 { out=$1; for(i=2;i<=NF;i++) { if(i!=3) { out=out"\t"$i} }; print out; out=""}' "${covsFile}.txt" > ${covsFile}.male.txt
for trait in gap_gm gap_wm gap_gwm; do
	traitFile="data/traits/${trait}"
	targetDir="results/${trait}/discovery/sex"
	for sex in male female; do
		awk 'NR==FNR { id[$1]; next } $1 in id { print }' ${covsFile}.${sex}.txt "${traitFile}.txt" > "${traitFile}.${sex}.txt"
		code/genetics/gwas.sh "${trait}" "${traitFile}.${sex}.txt" "${covsFile}.${sex}.txt" "${targetDir}/${sex}/gwas/" "${chrFileHandle}" "${snpQC}" "${impMFI}" "${maf}" "${threads}"
	done
done

# replication gwas in sex-stratified samples (European ancestry)
chrFileHandle="data/genetics/chr\${i}/imp_mri_qc_\${ancestry}/chr\${i}_mri_qc"
data="data/traits/replicate"
ancestriesCol="pan"
snpQC="data/genetics/meta/ukb_snp_qc.txt"
impMFI="data/genetics/ukb_imp_mfi"
maf=0.01
threads=100
ancestries="EUR"
covs="age,age2,ac1,ac2,ac3,TIV,array,PanC{1..20}"
awk -F'\t' 'NR==1 || $4==1 { out=$1; for(i=2;i<=NF;i++) { if(i!=4) { out=out"\t"$i} }; print out; out=""}' "${data}.txt" > ${data}.female.txt
awk -F'\t' 'NR==1 || $4==2 { out=$1; for(i=2;i<=NF;i++) { if(i!=4) { out=out"\t"$i} }; print out; out=""}' "${data}.txt" > ${data}.male.txt
for trait in gap_gm gap_wm gap_gwm; do
	targetDir="results/${trait}/replicate/sex"
	for sex in male female; do
		taskset -c 0-99 code/genetics/replicate.sh "${trait}" "${targetDir}/${sex}" "${chrFileHandle}" "${data}.${sex}.txt" "${covs}" "${ancestries}" "${ancestriesCol}" "${snpQC}" "${impMFI}" "${maf}" "${threads}"
	done
done

# meta-analysis across discovery and replication samples stratified by sex (European ancestry)
cohorts="DISCOV,REPLIC" # set output file handler & make sure that the order matches the list of sumstats (below)
type="ivweight" # type="nweight"
idCol="ID"
carryOver="CHR,BP"
leaveOneOut=0
(for trait in gap_gm gap_wm gap_gwm; do (
	for sex in male female; do (
		targetDir="results/${trait}/gwama/eur/sex/${sex}/"; mkdir -p ${targetDir}
		sumstats="results/${trait}/discovery/sex/${sex}/gwas/sumstats.txt.gz,results/${trait}/replicate/sex/${sex}/EUR/sumstats.txt.gz"
		code/genetics/metal.sh ${trait} ${targetDir} ${cohorts} ${sumstats} ${type} ${idCol} ${carryOver} ${leaveOneOut}
		) &
	done
	) &
done
wait)

	# quality-check METAL results
	nCriterion=TRUE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5
	hetP=1E-6 # heterogeneity p threshold for excluding SNPs
	(for trait in gap_gm gap_wm gap_gwm; do (
			for sex in male female; do (
				sumstats="results/${trait}/gwama/eur/sex/${sex}/metal.ivweight.gz"
				outFile="results/${trait}/gwama/eur/sex/${sex}/metal.ivweight.qc"
				Rscript code/genetics/metal.qc.R ${sumstats} ${outFile} ${nCriterion} ${hetP}
			) &
		done
		) & 
	done
	wait)

	# run ld score regression across sex-stratified meta-analysis results
	conda activate envs/ldsc
	LDsnplist="data/ldsc/resources/w_hm3.noMHC.snplist"
	LDchr="data/ldsc/resources/eur_w_ld_chr/"
	LDbaseline="data/ldsc/resources/baselineLD/baselineLD."
	LDweights="data/ldsc/resources/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
	LDfrqfiles="data/ldsc/resources/1000G_Phase3_frq/1000G.EUR.QC."
	partitioned=0
	sumstatsCols="CHR,BP,ID,A1,A2,BETA,SE,P,N"
	sumstatsColNames="CHR,POS,SNP,A1,A2,BETA,SE,P,N"
	(for trait in gap_gm gap_wm gap_gwm; do (
		for sex in male female; do (
	    	targetDir="results/${trait}/gwama/eur/sex/${sex}/ldsc"
	    	sumstats="results/${trait}/gwama/eur/sex/${sex}/metal.ivweight.qc.gz"
			./code/genetics/ldsc.sh "${trait}" "${targetDir}" "${sumstats}" "${sumstatsCols}" "${sumstatsColNames}" "${LDsnplist}" "${LDchr}" "${LDbaselineH2}" "${LDbaselineTau}" "${LDweights}" "${LDfrqfiles}"
			) &
		done
		) &
	done
	wait)

			# combine estimates in one file
			conda activate envs/default
			traitList="gap_gm gap_wm gap_gwm"
			for sex in male female; do
				ldscFiles="results/gap_gm/gwama/eur/sex/${sex}/ldsc/ldsc.h2.results results/gap_wm/gwama/eur/sex/${sex}/ldsc/ldsc.h2.results results/gap_gwm/gwama/eur/sex/${sex}/ldsc/ldsc.h2.results"
				outFile="results/combined/gwama.eur.sex.${sex}.ldsc.h2.txt"
				./code/genetics/ldsc.combine.sh "${ldscFiles}" "${outFile}"
			done

			# calculate rg between male and female results
			conda activate envs/ldsc
			LDchr="data/ldsc/resources/eur_w_ld_chr/"
			xprefix="male."
			yprefix="female."
			threads=100 # set number of parallel workers
			(for trait in gap_gm gap_wm gap_gwm; do (
			    sumstatsX="results/${trait}/gwama/eur/sex/male/ldsc/ldsc.${trait}.sumstats.gz"
			    sumstatsY="results/${trait}/gwama/eur/sex/female/ldsc/ldsc.${trait}.sumstats.gz"
			    targetDir="results/${trait}/gwama/eur/sex/rg"
				./code/genetics/rg.sh "${targetDir}" "${sumstatsX}" "${sumstatsY}" "${xprefix}" "${yprefix}" "${LDchr}" "${threads}"
				) &
			done)

			# combine estimates of rg between male and female results
			awk 'NR==1 { print } FNR==1 { next } { print }' "results/gap_gm/gwama/eur/sex/rg/rg.results" "results/gap_wm/gwama/eur/sex/rg/rg.results" "results/gap_gwm/gwama/eur/sex/rg/rg.results" \
				> results/combined/gwama.eur.sex.rg

			# create supplementum table with h2 and rg
			header=trait$'\t'$'\t'male.snp_count$'\t'male.h2$'\t'male.h2_se$'\t'male.lambda_GC$'\t'male.chi2$'\t'male.intercept$'\t'male.intercept_se$'\t'male.ratio$'\t'male.ratio_se$'\t'$'\t'male.snp_count$'\t'male.h2$'\t'female.h2_se$'\t'female.lambda_GC$'\t'female.chi2$'\t'female.intercept$'\t'female.intercept_se$'\t'female.ratio$'\t'female.ratio_se$'\t'$'\t'rg$'\t'se$'\t'z$'\t'p$'\t'gcov_int$'\t'gcov_int_se
			awk -v header="${header}" 'BEGIN { print header } 
				FNR==1 { filenum++; next}
				filenum==1 { trait[FNR-1]=$1"\t\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 }
				filenum==2 { trait[FNR-1]=trait[FNR-1]"\t\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 }
				filenum==3 { trait[FNR-1]=trait[FNR-1]"\t\t"$3"\t"$4"\t"$5"\t"$6"\t"$11"\t"$12; print trait[FNR-1] }' results/combined/gwama.eur.sex.male.ldsc.h2.txt results/combined/gwama.eur.sex.female.ldsc.h2.txt results/combined/gwama.eur.sex.rg \
				> results/combined/gwama.eur.sex.rg.suppl.txt

# =====================================================
# === Explorative: plot ancestry components vs. BAG ===
# =====================================================

# set working directory and activate conda environment
cd /slow/projects/ukb_brainage
conda activate envs/pointdensity

# make plots for pc1-pc4
traitLabel="BAG (years)"
covsFile="data/traits/covs.txt"
covs="sex,age,age2,array,TIV" # do not use scanner site due to covariance with genetic principal components
PDadjust=0.1
PDsize=0.5
ymin=-10 # argument only sets annotation limit but not axis limits
ymax=10 # argument only sets annotation limit but not axis limits
ysteps=10
width=8.77
height=5
(for trait in gap_gm gap_wm gap_gwm; do (
	outFile="results/${trait}/discovery/ancestry/ancestry.png"
	traitFile="data/traits/${trait}.txt"
	mkdir -p results/${trait}/discovery/ancestry/
	Rscript code/genetics/ancestry.R "${outFile}" "${traitFile}" "${trait}" "${traitLabel}" "${covsFile}" "${covs}" "${PDadjust}" "${PDsize}" "${ymin}" "${ymax}" "${ysteps}" "${width}" "${height}"
	) &
done
wait)

	# combine plots
	plotTitles="Grey matter BAG,White matter BAG,Combined BAG"
	plots="results/gap_gm/discovery/ancestry/ancestry.png,results/gap_wm/discovery/ancestry/ancestry.png,results/gap_gwm/discovery/ancestry/ancestry.png"
	outputFile="results/combined/discovery.ancestry.png"
	width=8.77
	height=15
	ncol=1
	titleSize=13
	titlePosX=0.5
	titlePosY=1.01
	Rscript code/genetics/plot.combine.R "${plotTitles}" "${plots}" "${outputFile}" "${width}" "${height}" "${ncol}" "${titleSize}" "${titlePosX}" "${titlePosY}"

# ==========================================================================
# === Explorative: cross-trait associations of sex-stratified BAG traits ===
# ==========================================================================

# set working directory and activate conda environment
cd /slow/projects/ukb_brainage
conda activate envs/phesant

# settings
covs="age,age2,ac1,ac2,TIV"
phenoFile="data/basket/20210205_2007685/data/ukb45233_imaging_phesant" # phenoFile="data/basket/20200409_2007685/data/ukb41573_imaging_wba_phesant.csv"
phesantDir="/fast/software/PHESANT" # phesantDir="/home/groups/markett/software/PHESANT"
nparts=1
standardise="FALSE"
mawk 'NR==FNR { if($3==1 && $6==1) {id1[$1]}; next } FNR==1 { print; next } { id2=$1; gsub(/"/, "", id2); if(id2 in id1) { print } }' OFS=',' FS="\t" "results/mri/r2024.vars.txt" FS="," "${phenoFile}.csv" > "${phenoFile}.female.csv"
mawk 'NR==FNR { if($3==1 && $6==2) {id1[$1]}; next } FNR==1 { print; next } { id2=$1; gsub(/"/, "", id2); if(id2 in id1) { print } }' OFS=',' FS="\t" "results/mri/r2024.vars.txt" FS="," "${phenoFile}.csv" > "${phenoFile}.male.csv"

# run analysis
trait="gap_gm"; sex="male"; traitFile="data/${trait}/${trait}.txt"; covsFile="data/${trait}/covs.txt"; targetDir="results/${trait}/discovery/sex/${sex}/phesant"
	taskset -c 0-15 ./code/genetics/phesant.sh "${trait}" "${traitFile}" "${covsFile}" "${covs}" "${targetDir}" "${phenoFile}.${sex}.csv" "${phesantDir}" "${nparts}" "${standardise}"
trait="gap_gm"; sex="female"; traitFile="data/${trait}/${trait}.txt"; covsFile="data/${trait}/covs.txt"; targetDir="results/${trait}/discovery/sex/${sex}/phesant"
	taskset -c 16-30 ./code/genetics/phesant.sh "${trait}" "${traitFile}" "${covsFile}" "${covs}" "${targetDir}" "${phenoFile}.${sex}.csv" "${phesantDir}" "${nparts}" "${standardise}"
trait="gap_wm"; sex="male"; traitFile="data/${trait}/${trait}.txt"; covsFile="data/${trait}/covs.txt"; targetDir="results/${trait}/discovery/sex/${sex}/phesant"
	taskset -c 31-45 ./code/genetics/phesant.sh "${trait}" "${traitFile}" "${covsFile}" "${covs}" "${targetDir}" "${phenoFile}.${sex}.csv" "${phesantDir}" "${nparts}" "${standardise}"
trait="gap_wm"; sex="female"; traitFile="data/${trait}/${trait}.txt"; covsFile="data/${trait}/covs.txt"; targetDir="results/${trait}/discovery/sex/${sex}/phesant"
	taskset -c 46-60 ./code/genetics/phesant.sh "${trait}" "${traitFile}" "${covsFile}" "${covs}" "${targetDir}" "${phenoFile}.${sex}.csv" "${phesantDir}" "${nparts}" "${standardise}"
trait="gap_gwm"; sex="male"; traitFile="data/${trait}/${trait}.txt"; covsFile="data/${trait}/covs.txt"; targetDir="results/${trait}/discovery/sex/${sex}/phesant"
	taskset -c 61-75 ./code/genetics/phesant.sh "${trait}" "${traitFile}" "${covsFile}" "${covs}" "${targetDir}" "${phenoFile}.${sex}.csv" "${phesantDir}" "${nparts}" "${standardise}"
trait="gap_gwm"; sex="female"; traitFile="data/${trait}/${trait}.txt"; covsFile="data/${trait}/covs.txt"; targetDir="results/${trait}/discovery/sex/${sex}/phesant"
	taskset -c 76-90 ./code/genetics/phesant.sh "${trait}" "${traitFile}" "${covsFile}" "${covs}" "${targetDir}" "${phenoFile}.${sex}.csv" "${phesantDir}" "${nparts}" "${standardise}"

	# create plots and get result summary
	imaging=FALSE
	ylim=60
	ysteps=10
	multipleTesting=both
	width=8.7
	height=5.0
	repel_nudge_y=10
	(for trait in gap_gm gap_wm gap_gwm; do (
		for sex in male female; do (
			phesantResults="results/${trait}/discovery/sex/${sex}/phesant/phesant.output/results-combined.txt"
			targetDir="results/${trait}/discovery/sex/${sex}/phesant"
			code/genetics/phesant.output.R "${trait}" "${phesantResults}" "${targetDir}" "${imaging}" "${multipleTesting}" "${ylim}" "${ysteps}" "${width}" "${height}" "${repel_nudge_y}"
			) &
		done
		) &
	done)

	# combine phesant plots and result summaries across traits
	conda activate envs/default
	for sex in male female; do
		traits="gap_gm,gap_wm,gap_gwm"
		phesantSummary="results/gap_gm/discovery/sex/${sex}/phesant/phesant.summary.txt,results/gap_wm/discovery/sex/${sex}/phesant/phesant.summary.txt,results/gap_gwm/discovery/sex/${sex}/phesant/phesant.summary.txt"
		phesantPlot="results/gap_gm/discovery/sex/${sex}/phesant/phesant.png,results/gap_wm/discovery/sex/${sex}/phesant/phesant.png,results/gap_gwm/discovery/sex/${sex}/phesant/phesant.png"
		outputFile="results/combined/discovery.sex.${sex}.phewas"
		code/genetics/phesant.combine.R "${traits}" "${phesantSummary}" "${phesantPlot}" "${outputFile}"
	done

	# combine results from male and female samples
	traits="gap_gm,gap_wm,gap_gwm"
	phesantSummary="results/combined/discovery.sex.male.phewas.txt,results/combined/discovery.sex.female.phewas.txt"
	outputFile="results/combined/discovery.sex.phewas"
	code/genetics/phesant.combine.sex.R "${traits}" "${phesantSummary}" "${outputFile}"

		# plot z scores from male vs. female analysis
		conda activate envs/pointdensity
		input1="results/combined/discovery.sex.phewas.txt"
		input2="results/combined/discovery.sex.phewas.txt"
		matchCol="varName"
		conversion="p-to-z"
		PDadjust=0.375
		PDsize=0.5
		xlabel="Males (z-scores)"
		ylabel="Females (z-scores)"
		xmin=-13
		xmax=13
		ymin=-12
		ymax=12
		width=2.9
		height=2.8
		preview=FALSE
		(for trait in gap_gm gap_wm gap_gwm; do (
			outFile="results/${trait}/discovery/sex/phesant.zscores.png"
			input1Col="a_${trait}_pvalue"
			input2Col="b_${trait}_pvalue"
			input1SignCol="a_${trait}_rho"
			input2SignCol="b_${trait}_rho"
			./code/genetics/compareZ.R "${outFile}" "${input1}" "${input2}" "${input1Col}" "${input2Col}" "${input1SignCol}" "${input2SignCol}" "${matchCol}" "${conversion}" "${PDadjust}" "${PDsize}" "${xlabel}" "${ylabel}" "${xmin}" "${xmax}" "${ymin}" "${ymax}" "${width}" "${height}" "${preview}"
			) &
		done
		wait)

			# combine Z score plots
			plotTitles="Grey matter BAG,White matter BAG,Combined BAG"
			plots="results/gap_gm/discovery/sex/phesant.zscores.png,results/gap_wm/discovery/sex/phesant.zscores.png,results/gap_gwm/discovery/sex/phesant.zscores.png"
			outputFile="results/combined/discovery.sex.phewas.zscores.png"
			width=8.7
			height=2.9
			ncol=3
			titleSize=10
			titlePosX=0.53
			titlePosY=1
			Rscript code/genetics/plot.combine.R "${plotTitles}" "${plots}" "${outputFile}" "${width}" "${height}" "${ncol}" "${titleSize}" "${titlePosX}" "${titlePosY}"

		# create qq-plots from male vs. female analysis
		conda activate envs/default
		nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5
		prune=FALSE # ommit overplotting of variants with pval > 0.01
		drawCI=FALSE # draw confidence interval?
		drawLambda=FALSE # draw LambdaGC?
		xend=6 # upper x axis limit
		xsteps=2 # x axis breaks
		yend=8 # upper y axis limit
		ysteps=2 # y axis breaks
		width=3 # plot width (4 inch)
		height=3 # plot height (4 inch)
		(for trait in gap_gm gap_wm gap_gwm; do (
			targetDir="results/${trait}/discovery/sex/"
			pCol="${trait}_deltaP" # column containing p-values
			sumstats="results/combined/discovery.sex.phewas.txt"
			Rscript code/genetics/qqplot.R "${trait}" "${targetDir}" "${sumstats}" "${pCol}" "${nCriterion}" "${prune}" "${drawCI}" "${drawLambda}" "${xend}" "${xsteps}" "${yend}" "${ysteps}" "${width}" "${height}"
			mv "results/${trait}/discovery/sex/qqplot.png" "results/${trait}/discovery/sex/phesant.qqplot.png" 
			mv "results/${trait}/discovery/sex/qqplot.log" "results/${trait}/discovery/sex/phesant.qqplot.log" 
			) &
		done
		wait)

			# combine qq plots
			plotTitles="Grey matter BAG,White matter BAG,Combined BAG"
			plots="results/gap_gm/discovery/sex/phesant.qqplot.png,results/gap_wm/discovery/sex/phesant.qqplot.png,results/gap_gwm/discovery/sex/phesant.qqplot.png"
			outputFile="results/combined/discovery.sex.phewas.qqplot.png"
			width=8.5
			height=2.6
			ncol=3
			titleSize=10
			titlePosX=0.53
			titlePosY=1
			Rscript code/genetics/plot.combine.R "${plotTitles}" "${plots}" "${outputFile}" "${width}" "${height}" "${ncol}" "${titleSize}" "${titlePosX}" "${titlePosY}"

# calculate FreeSurfer correlations
conda activate envs/default
covNames="age,age2,ac1,ac2,TIV"
freesurferFile="results/mri/freesurfer.tsv.gz"
(for trait in gap_gm gap_wm gap_gwm; do (
	for sex in male female; do (
		traitFile="data/traits/${trait}.${sex}.txt"
		covsFile="data/traits/covs.${sex}.txt"
		targetDir="results/${trait}/discovery/sex/${sex}/surfcorr"
		code/mri/surfcorr.R "${trait}" "${traitFile}" "${covsFile}" "${covNames}" "${freesurferFile}" "${targetDir}"
		) &
	done
	) &
done
wait)

	# combine surfcorr results
	for sex in male female; do
		traits="gap_gm,gap_wm,gap_gwm"
		surfcorrFiles="results/gap_gm/discovery/sex/${sex}/surfcorr/surfcorr.txt,results/gap_wm/discovery/sex/${sex}/surfcorr/surfcorr.txt,results/gap_gwm/discovery/sex/${sex}/surfcorr/surfcorr.txt"
		outFile="results/combined/discovery.sex.${sex}.surfcorr.txt"
		code/mri/surfcorr.combine.R "${traits}" "${surfcorrFiles}" "${outFile}"
	done

	# combine results from male and female samples
	traits="gap_gm,gap_wm,gap_gwm"
	surfcorrFiles="results/combined/discovery.sex.male.surfcorr.txt,results/combined/discovery.sex.female.surfcorr.txt"
	outFile="results/combined/discovery.sex.surfcorr.txt"
	code/mri/surfcorr.combine.sex.R "${traits}" "${surfcorrFiles}" "${outFile}"

		# plot z scores from male vs. female analysis
		conda activate envs/pointdensity
		input1="results/combined/discovery.sex.surfcorr.txt"
		input2="results/combined/discovery.sex.surfcorr.txt"
		matchCol="id"
		conversion="none"
		PDadjust=0.01
		PDsize=0.5
		xlabel="Correlations (rho) in males"
		ylabel="Correlations (rho) in females"
		xmin=-0.36
		xmax=0.36
		ymin=-0.36
		ymax=0.36
		width=2.9
		height=2.8
		preview=FALSE
		(for trait in gap_gm gap_wm gap_gwm; do (
			rm -f "results/${trait}/discovery/sex/freesurfer.zscores.png"
			outFile="results/${trait}/discovery/sex/surfcorr.zscores.png"
			input1Col="a_${trait}_rho"
			input2Col="b_${trait}_rho"
			input1SignCol="a_${trait}_rho"
			input2SignCol="b_${trait}_rho"
			./code/genetics/compareZ.R "${outFile}" "${input1}" "${input2}" "${input1Col}" "${input2Col}" "${input1SignCol}" "${input2SignCol}" "${matchCol}" "${conversion}" "${PDadjust}" "${PDsize}" "${xlabel}" "${ylabel}" "${xmin}" "${xmax}" "${ymin}" "${ymax}" "${width}" "${height}" "${preview}"
			) &
		done
		wait)

			# combine Z score plots
			plotTitles="Grey matter BAG,White matter BAG,Combined BAG"
			plots="results/gap_gm/discovery/sex/surfcorr.zscores.png,results/gap_wm/discovery/sex/surfcorr.zscores.png,results/gap_gwm/discovery/sex/surfcorr.zscores.png"
			outputFile="results/combined/discovery.sex.surfcorr.zscores.png"
			width=8.7
			height=2.9
			ncol=3
			titleSize=10
			Rscript code/genetics/compareZ.combine.R "${plotTitles}" "${plots}" "${outputFile}" "${width}" "${height}" "${ncol}" "${titleSize}"

		# create qq-plots from male vs. female analysis
		conda activate envs/default
		nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5
		prune=FALSE # ommit overplotting of variants with pval > 0.01
		drawCI=FALSE # draw confidence interval?
		drawLambda=FALSE # draw LambdaGC?
		xend=3 # upper x axis limit
		xsteps=1 # x axis breaks
		yend=10 # upper y axis limit
		ysteps=2 # y axis breaks
		width=3 # plot width (4 inch)
		height=3 # plot height (4 inch)
		(for trait in gap_gm gap_wm gap_gwm; do (
			targetDir="results/${trait}/discovery/sex/"
			pCol="${trait}_deltaP" # column containing p-values
			sumstats="results/combined/discovery.sex.surfcorr.txt"
			Rscript code/genetics/qqplot.R "${trait}" "${targetDir}" "${sumstats}" "${pCol}" "${nCriterion}" "${prune}" "${drawCI}" "${drawLambda}" "${xend}" "${xsteps}" "${yend}" "${ysteps}" "${width}" "${height}"
			\mv "results/${trait}/discovery/sex/qqplot.png" "results/${trait}/discovery/sex/freesurfer.qqplot.png" 
			\mv "results/${trait}/discovery/sex/qqplot.log" "results/${trait}/discovery/sex/freesurfer.qqplot.log" 
			) &
		done
		wait)

			# combine qq plots
			plotTitles="Grey matter BAG,White matter BAG,Combined BAG"
			plots="results/gap_gm/discovery/sex/freesurfer.qqplot.png,results/gap_wm/discovery/sex/freesurfer.qqplot.png,results/gap_gwm/discovery/sex/freesurfer.qqplot.png"
			outputFile="results/combined/discovery.sex.surfplot.qqplot.png"
			width=8.5
			height=2.6
			ncol=3
			Rscript code/genetics/compareZ.combine.R "${plotTitles}" "${plots}" "${outputFile}" "${width}" "${height}" "${ncol}"

# =============================================
# === prepare summary statistics for zenodo ===
# =============================================

# set working directory
cd /slow/projects/ukb_brainage
conda activate envs/default

# prepare summary statistics for zenodo
targetDir="zenodo"
mkdir -p "${targetDir}"
for tissue in gm wm gwm; do
	echo "Starting with gap_${tissue}"

	# discovery
	pigz -dc "results/gap_${tissue}/discovery/gwas/sumstats.txt.gz" > "${targetDir}/brainage2025.discov.${tissue}"
	chmod 770 "${targetDir}/brainage2025.discov.${tissue}"
	pigz -f "${targetDir}/brainage2025.discov.${tissue}"

	# replication
	pigz -dc "results/gap_${tissue}/replicate/metal.eur/metal.ivweight.qc.gz" > "${targetDir}/brainage2025.replic.eur.${tissue}"
	chmod 770 "${targetDir}/brainage2025.replic.eur.${tissue}"
	pigz -f "${targetDir}/brainage2025.replic.eur.${tissue}"
	pigz -dc "results/gap_${tissue}/replicate/mrmega.all/mrmega.weights.gz" > "${targetDir}/brainage2025.replic.multi.${tissue}"
	chmod 770 "${targetDir}/brainage2025.replic.multi.${tissue}"
	pigz -f "${targetDir}/brainage2025.replic.multi.${tissue}"

	# discovery & replication
	pigz -dc "results/gap_${tissue}/gwama/eur/metal.ivweight.qc.gz" > "${targetDir}/brainage2025.full.eur.${tissue}"
	chmod 770 "${targetDir}/brainage2025.full.eur.${tissue}"
	pigz -f "${targetDir}/brainage2025.full.eur.${tissue}"
	pigz -dc "results/gap_${tissue}/gwama/all/mrmega.weights.gz" > "${targetDir}/brainage2025.full.multi.${tissue}"
	chmod 770 "${targetDir}/brainage2025.full.multi.${tissue}"
	pigz -f "${targetDir}/brainage2025.full.multi.${tissue}"

	# pgsweights
	pigz -dc "results/gap_${tissue}/gwama/eur/sbayes/sbayesrc.snpRes.gz" > "${targetDir}/brainage2025.pgsweights.${tissue}"
	chmod 770 "${targetDir}/brainage2025.pgsweights.${tissue}"
	pigz -f "${targetDir}/brainage2025.pgsweights.${tissue}"
done


