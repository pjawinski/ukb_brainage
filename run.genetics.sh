#!/bin/bash

# ==================================
# === install conda environments ===
# ==================================

for env in default eqtl genesis gofuncr ldsc locuszoom phesant pops power xgb; do
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
	vars="FID,IID,brainage_${trait}_stk"
	varsRename="FID,IID,${trait}"
	outFile="data/traits/${trait}.txt"
	Rscript code/genetics/getvars.R ${subsFile} ${varsFile} ${idCol} ${vars} ${varsRename} ${outFile} 
done

# create discovery covs file
subsFile='data/genetics/chr1/imp_mri_qc/chr1_mri_qc.psam'
varsFile='results/genetics/r2022.vars.pca.txt'
idCol='IID'
vars='FID,IID,sex,t1.age,t1.age2,t1.ac1,t1.ac2,array,t1.TIV,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20'
varsRename='FID,IID,sex,age,age2,ac1,ac2,array,TIV,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20'
outFile="data/traits/covs.txt"
Rscript code/genetics/getvars.R ${subsFile} ${varsFile} ${idCol} ${vars} ${varsRename} ${outFile} 

# create all-in-one replication file
subsFile='results/mri/iid.replication.txt'
varsFile='results/mri/brainage.replication.singlemodel.txt'
idCol='IID'
vars='FID,IID,pan,sex,t1_age,t1_age2,t1_ac1,t1_ac2,t1_ac3,t1_TIV,array,PanC1,PanC2,PanC3,PanC4,PanC5,PanC6,PanC7,PanC8,PanC9,PanC10,PanC11,PanC12,PanC13,PanC14,PanC15,PanC16,PanC17,PanC18,PanC19,PanC20,brainage_gap_gm_stack,brainage_gap_wm_stack,brainage_gap_gwm_stack'
varsRename='FID,IID,pan,sex,age,age2,ac1,ac2,ac3,TIV,array,PanC1,PanC2,PanC3,PanC4,PanC5,PanC6,PanC7,PanC8,PanC9,PanC10,PanC11,PanC12,PanC13,PanC14,PanC15,PanC16,PanC17,PanC18,PanC19,PanC20,gap_gm,gap_wm,gap_gwm'
outFile="data/traits/replicate.txt"
Rscript code/genetics/getvars.R ${subsFile} ${varsFile} ${idCol} ${vars} ${varsRename} ${outFile} 

# ============================================
# === run cross-trait association analysis ===
# ============================================

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
trait="gap_gm"; traitFile="data/${trait}/${trait}.txt"; covsFile="data/${trait}/covs.txt"; targetDir="results/${trait}/phesant"
	taskset -c 0-24 ./code/genetics/phesant.sh "${trait}" "${traitFile}" "${covsFile}" "${covs}" "${targetDir}" "${phenoFile}" "${phesantDir}" "${nparts}" "${standardise}"

trait="gap_wm"; traitFile="data/${trait}/${trait}.txt"; covsFile="data/${trait}/covs.txt"; targetDir="results/${trait}/phesant"
	taskset -c 25-49 ./code/genetics/phesant.sh "${trait}" "${traitFile}" "${covsFile}" "${covs}" "${targetDir}" "${phenoFile}" "${phesantDir}" "${nparts}" "${standardise}"

trait="gap_gwm"; traitFile="data/${trait}/${trait}.txt"; covsFile="data/${trait}/covs.txt"; targetDir="results/${trait}/phesant"
	taskset -c 50-74 ./code/genetics/phesant.sh "${trait}" "${traitFile}" "${covsFile}" "${covs}" "${targetDir}" "${phenoFile}" "${phesantDir}" "${nparts}" "${standardise}"

# create plots and get result summary
imaging=TRUE
ylim=60
ysteps=10
multipleTesting=both
width=8.7
height=5.0
repel_nudge_y=10
for trait in gap_gm gap_wm gap_gwm; do
	phesantResults="results/${trait}/phesant/phesant.output/results-combined.txt"
	targetDir="results/${trait}/phesant_imaging/"
	code/genetics/phesant.output.R "${trait}" "${phesantResults}" "${targetDir}" "${imaging}" "${multipleTesting}" "${ylim}" "${ysteps}" "${width}" "${height}" "${repel_nudge_y}" 
done

# combine phesant plots and result summaries across traits
traits="gap_gm,gap_wm,gap_gwm"
phesantSummary="results/gap_gm/phesant_imaging/phesant.summary.txt,results/gap_wm/phesant_imaging/phesant.summary.txt,results/gap_gwm/phesant_imaging/phesant.summary.txt"
phesantPlot="results/gap_gm/phesant_imaging/phesant.png,results/gap_wm/phesant_imaging/phesant.png,results/gap_gwm/phesant_imaging/phesant.png"
outputFile="results/combined/phewas.imaging"
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
	traitFile="data/${trait}/${trait}.txt"
	covsFile="data/${trait}/covs.txt"
	targetDir="results/${trait}/surfplot"
	code/mri/surfcorr.R "${trait}" "${traitFile}" "${covsFile}" "${covNames}" "${freesurferFile}" "${targetDir}"
	) &
done)

# combine surfcorr results
traits="gap_gm,gap_wm,gap_gwm"
surfcorrFiles="results/gap_gm/surfplot/surfcorr.txt,results/gap_wm/surfplot/surfcorr.txt,results/gap_gwm/surfplot/surfcorr.txt"
outFile="results/combined/surfcorr.txt"
code/mri/surfcorr.combine.R "${traits}" "${surfcorrFiles}" "${outFile}"

# make plots
(for trait in gap_gm gap_wm gap_gwm; do (
	rm -rf "results/${trait}/surfplot"
	mkdir -p "results/${trait}/surfplot"
	matlab -r "\
		workingDir = '$(pwd)';\
		ENIGMAtoolboxPath = '/fast/software/ENIGMA/matlab';\
		surfCorr = 'results/combined/surfcorr.txt';\
		mappingFileSurfArea = 'code/mri/surfplot.mapping.surfarea.txt';\
		mappingFileThickAvg = 'code/mri/surfplot.mapping.thickavg.txt';\
		mappingFileGrayVol = 'code/mri/surfplot.mapping.grayvol.txt';\
		mappingFileSubcortical = 'code/mri/surfplot.mapping.subcortical.txt';\
		rCol = '${trait}_rho';\
		pCol = '${trait}_fdr';\
		barTitle = 'Correlation (r)';\
		outFile = 'results/${trait}/surfplot/surfplot.png';\
		colorBar = 'horizontal';\
		cblim = '0.35';\
		run code/mri/surfplot.m;\
		exit"
		) &
done
wait)

# combine surfplots for supplementum
traits="gap_gm,gap_wm,gap_gwm"
surfPlots="results/gap_gm/surfplot/surfplot.png,results/gap_wm/surfplot/surfplot.png,results/gap_gwm/surfplot/surfplot.png"
surfPlotsTXT="results/gap_gm/surfplot/surfplot.txt,results/gap_wm/surfplot/surfplot.txt,results/gap_gwm/surfplot/surfplot.txt"
outputPrefix="results/combined/surfplot"
width=8
height=10
Rscript code/mri/surfplot.combine.R "${traits}" "${surfPlots}" "${surfPlotsTXT}" "${outputPrefix}" "${width}" "${height}"

# combine gwm PHESANT and Surface plot for main article
phesantplot="results/gap_gwm/phesant/phewas.png"
surfplot="results/gap_gwm/surfplot/surfplot.png"
outFile="results/combined/surfplot.phesant.png"
Rscript code/mri/surfplot.phesant.combine.R "${phesantplot}" "${surfplot}" "${outFile}"

# =================================
# === run heritability analysis ===
# =================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage
conda activate envs/default

# general settings
covsDiscr=sex,ac1,ac2,array
covsQuant=age,age2,TIV,PC{1..20}
grm=data/grm/ukb_merged
threads=50

# run analysis
for trait in gap_gm gap_wm gap_gwm; do
	traitFile="data/${trait}/${trait}.txt"
	covsFile="data/${trait}/covs.txt"
	targetDir="results/${trait}/greml"
	./code/genetics/greml.sh "${trait}" "${traitFile}" "${covsFile}" "${targetDir}" "${covsDiscr}" "${covsQuant}" "${grm}" "${threads}" 
done

# ============================================
# === run genome-wide association analysis ===
# ============================================

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
	targetDir="results/${trait}/gwas/"
	code/genetics/gwas.sh "${trait}" "${traitFile}" "${covsFile}" "${targetDir}" "${chrFileHandle}" "${snpQC}" "${impMFI}" "${maf}" "${threads}" 
done

# ===============================
# === run ld score regression ===
# ===============================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage
conda activate envs/ldsc

# general settings
LDsnplist="data/ldsc/resources/w_hm3.noMHC.snplist"
LDchr="data/ldsc/resources/eur_w_ld_chr/"
LDbaseline="data/ldsc/resources/baselineLD/baselineLD."
LDweights="data/ldsc/resources/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
LDfrqfiles="data/ldsc/resources/1000G_Phase3_frq/1000G.EUR.QC."

# run analysis
for trait in gap_gm gap_wm gap_gwm; do
    targetDir="results/${trait}/ldsc"
    sumstats="results/${trait}/gwas/sumstats.txt.gz"
    ./code/genetics/ldsc/ldsc.sh "${trait}" "${targetDir}" "${sumstats}" "${LDsnplist}" "${LDchr}" "${LDbaseline}" "${LDweights}" "${LDfrqfiles}"
done

# combine estimates in one file
conda activate envs/default
traitList="gap_gm gap_wm gap_gwm"
ldscFiles="results/gap_gm/ldsc/ldsc.h2.results results/gap_wm/ldsc/ldsc.h2.results results/gap_gwm/ldsc/ldsc.h2.results"
outFile="results/combined/ldsc.h2.txt"
./code/genetics/ldsc.combine.sh "${ldscFiles}" "${outFile}"

# =====================================
# === plot partitioned heritability ===
# =====================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage/
conda activate envs/default

# make partitioned heritability plot and get summary table
ymin=-10
ymax=25
ysteps=10
width=8.1
height=3.7
(for trait in gap_gm gap_wm gap_gwm; do (
    ldscFile="results/${trait}/ldsc/ldsc.partitioned.results"
    Rscript code/genetics/ldsc.partitioned.R "${ldscFile}" "${ymin}" "${ymax}" "${ysteps}" "${width}" "${height}"
    ) &
done)
wait

# combine plots and summary tables
traits="gap_gm,gap_wm,gap_gwm"
plotTitles="Grey_matter,White_matter,Grey_and_white_matter"
ldscPlots="results/gap_gm/ldsc/ldsc.partitioned.results.png,results/gap_wm/ldsc/ldsc.partitioned.results.png,results/gap_gwm/ldsc/ldsc.partitioned.results.png"
ldscSummary="results/gap_gm/ldsc/ldsc.partitioned.results.summary.txt,results/gap_wm/ldsc/ldsc.partitioned.results.summary.txt,results/gap_gwm/ldsc/ldsc.partitioned.results.summary.txt"
outputFile="results/combined/ldsc.partitioned"
width=8.1
height=11.1
Rscript code/genetics/ldsc.partitioned.combined.R "${traits}" "${plotTitles}" "${ldscPlots}" "${ldscSummary}" "${outputFile}" "${width}" "${height}"

# ===========================================
# === Run ANNOVAR enrichment test in FUMA ===
# ===========================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage
conda activate envs/default

# prepare summary statistics for FUMA upload
for trait in gap_gm gap_wm gap_gwm; do (
    targetDir="results/${trait}/fuma"
    sumstatsPLINK="results/${trait}/gwas/sumstats_plink.txt.gz"
    ./code/genetics/fuma/fuma.sh "${trait}" "${targetDir}" "${sumstatsPLINK}") &
done
wait

# combine results of ANNOVAR enrichment test
traits="gap_gm,gap_wm,gap_gwm"
inFiles="results/gap_gm/fuma/FUMA_job174874/annov.stats.txt,results/gap_wm/fuma/FUMA_job174875/annov.stats.txt,results/gap_gwm/fuma/FUMA_job174876/annov.stats.txt"
outFile="results/combined/annovar.txt"
./code/genetics/fuma.annov.combine.R "${traits}" "${inFiles}" "${outFile}"

# make annovar enrichment plot for all traits
annovarStats="results/gap_gm/fuma/FUMA_job174874.zip,results/gap_wm/fuma/FUMA_job174875.zip,results/gap_gwm/fuma/FUMA_job174876.zip"
outputFile="results/combined/fuma.annov.stats.png"
width=8.04
height=7.50
Rscript code/genetics/fuma.annov.R "${annovarStats}" "${outputFile}" "${width}" "${height}"

# make enrichment plot of single trait for main article
annovarStats="results/gap_gwm/fuma/FUMA_job174876.zip"
outputFile="results/combined/fuma.annov.stats.gwm.png"
width=8.04
height=2.90
Rscript code/genetics/fuma.annov.R "${annovarStats}" "${outputFile}" "${width}" "${height}"

# ============================
# === run ApoE-e4 analysis ===
# ============================

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

# ================================
# === run conditional analysis ===
# ================================

# set working directory and load conda environment
cd /slow/projectss/ukb_brainage
conda activate envs/default

# run conditional analysis and do positional annotation (!! do not run in parallel across traits)
subsFile='data/genetics/chr1/imp_mri_qc/chr1_mri_qc.psam'
chrFileHandle="data/genetics/chr\$i/imp_mri_qc/bed/chr\${i}_mri_qc"
cols="ID,A1,A2,A1_FREQ,BETA,SE,P,N,CHR,BP,Z,ETA2,INFO,TYPED" # first 10 columns have the following order ID,A1,A2,A1_FREQ,BETA,SE,P,N,CHR,BP
pthresh=1E-6
threads=100
humandb="data/annovar/humandb"
refseq="data/refseq/GRCh37_latest_genomic.edit.gff.gz"

# run analysis
for trait in gap_gm gap_wm gap_gwm; do
	targetDir="results/${trait}/conditional"
	sumstats="results/${trait}/gwas/sumstats.txt.gz"
	sumstatsPLINK="results/${trait}/gwas/sumstats_plink.txt.gz"
	./code/genetics/conditional.sh "${subsFile}" "${targetDir}" "${chrFileHandle}" "${sumstats}" "${cols}" "${sumstatsPLINK}" "${pthresh}" "${threads}" "${humandb}" "${refseq}"
done

# =============================
# === create regional plots ===
# =============================

# set working directory and load conda environment
cd /slow/projectss/ukb_brainage
conda activate envs/locuzsoom

# general settings
geneticsDir="data/genetics"
pthresh=5E-8
flank="500kb"

# gap_gm
trait="gap_gm"
traitDescription="grey matter"
targetDir="results/${trait}/locuszoom/"
sumstats="results/${trait}/gwas/sumstats.txt.gz"
conditionalFile="results/${trait}/conditional/conditional.cleaned.tophits.annovar.txt"

	# optional settings
	snpsManual="rs796226228 rs376899426 rs377158217 rs76122535 rs56072903"
	flankManual=""
	chrManual="1 2 2 2 17"
	startManual="214800000 197500000 200200000 203000000 43000000"
	endManual="216200000 199500000 201600000 205000000 45500000"
	./code/genetics/locuszoom.sh "${trait}" "${traitDescription}" "${targetDir}" "${geneticsDir}" "${sumstats}" "${conditionalFile}" "${pthresh}" "${flank}" "${snpsManual}" "${flankManual}" "${chrManual}" "${startManual}" "${endManual}" &

# gap_wm
trait="gap_wm"
traitDescription="white matter"
targetDir="results/${trait}/locuszoom/"
sumstats="results/${trait}/gwas/sumstats.txt.gz"
conditionalFile="results/${trait}/conditional/conditional.cleaned.tophits.annovar.txt"

	# optional settings
	snpsManual="rs10932008 rs12287076 rs111854640"
	flankManual=""
	chrManual="2 11 17"
	startManual="203000000 46800000 43000000"
	endManual="205000000 48600000 45500000"
	./code/genetics/locuszoom/locuszoom.sh "${trait}" "${traitDescription}" "${targetDir}" "${geneticsDir}" "${sumstats}" "${conditionalFile}" "${pthresh}" "${flank}" "${snpsManual}" "${flankManual}" "${chrManual}" "${startManual}" "${endManual}" &

# gap_gwm
trait="gap_gwm"
traitDescription="grey and white matter"
targetDir="results/${trait}/locuszoom/"
sumstats="results/${trait}/gwas/sumstats.txt.gz"
conditionalFile="results/${trait}/conditional/conditional.cleaned.tophits.annovar.txt"

	# optional settings
	snpsManual="s550778111 rs377158217 rs76122535 rs199737620"
	flankManual=""
	chrManual="1 2 2 17"
	startManual="43300000 200200000 203000000 43000000"
	endManual="445000000 201600000 205000000 45500000"
	./code/genetics/locuszoom/locuszoom.sh "${trait}" "${traitDescription}" "${targetDir}" "${geneticsDir}" "${sumstats}" "${conditionalFile}" "${pthresh}" "${flank}" "${snpsManual}" "${flankManual}" "${chrManual}" "${startManual}" "${endManual}"

# ========================================
# === Get 95% credible set of variants ===
# ========================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage/
conda activate envs/default

# general settings
geneticsDir="data/genetics"
pthresh=5E-8

# run analysis
for trait in gap_gm gap_wm gap_gwm; do (
	targetDir="results/${trait}/credibleSet"
	LDsample="data/${trait}/${trait}.txt"
	conditionalFile="results/${trait}/conditional/conditional.cleaned.tophits.annovar.txt"
	sumstats="results/${trait}/gwas/sumstats.txt.gz"
	sumstatsPLINK="results/${trait}/gwas/sumstats_plink.txt.gz"
	./code/genetics/finemap.sh "${trait}" "${targetDir}" "${geneticsDir}" "${LDsample}" "${conditionalFile}" "${sumstats}" "${sumstatsPLINK}" "${pthresh}"
	) &
done
wait

# do positional gene mapping
chr_bp_ref_alt="chromosome,bp,allele1,allele2"
humandb="data/annovar/humandb"
refseq="data/refseq/GRCh37_latest_genomic.edit.gff.gz"
for trait in gap_gm gap_wm gap_gwm; do (
	inputFile="results/${trait}/credibleSet/credibleSet.df.txt"
	outputFile="results/${trait}/credibleSet/credibleSet.df.annot.txt"
	code/genetics/annotation.sh "${inputFile}" "${outputFile}" "${chr_bp_ref_alt}" "${humandb}" "${refseq}"
) &
done
wait

# add relevant genes (sorted by cumulative PP of variants from finemapping) to credible.ls.txt file
for trait in gap_gm gap_wm gap_gwm; do (
	dfFile="results/${trait}/credibleSet/credibleSet.df.annot.txt"
	lsFile="results/${trait}/credibleSet/credibleSet.ls.txt"
	outputFile="results/${trait}/credibleSet/credibleSet.ls.genes.txt"
	code/genetics/finemap.genes.R "${dfFile}" "${lsFile}" "${outputFile}"
) &
done
wait

# create summary file for exonic variants
locusId="index_rsid"
id="rsid"
beta_se="beta,se"
pthresh=1
for trait in gap_gm gap_wm gap_gwm; do (
	nsynFile="results/${trait}/credibleSet/credibleSet.df.annot.nonsynonymous.txt"
	outputFile="results/${trait}/credibleSet/credibleSet.df.annot.nonsynonymous.summary.txt"
	code/genetics/annotation.exonic.R "${nsynFile}" "${locusId}" "${id}" "${beta_se}" "${pthresh}" "${outputFile}" 
) &
done
wait

# ==================================
# === run eqtl/sqtl SMR analysis ===
# ==================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage/
conda activate envs/default

# run analysis
geneticsDir="data/genetics"
pthresh=1E-6
threads=100
clumpR2=0.8
smrCorrect="fdr"
for trait in gap_gm gap_wm gap_gwm; do
	targetDir="results/${trait}/smr"
	sumstats="results/${trait}/gwas/sumstats.txt.gz"
	conditionalFile="results/${trait}/conditional/conditional.cleaned.tophits.annovar.txt"
	LDsample="data/${trait}/${trait}.txt"
	./code/genetics/smr.sh "${trait}" "${targetDir}" "${geneticsDir}" "${sumstats}" "${conditionalFile}" "${pthresh}" "${LDsample}" "${threads}" "${clumpR2}" "${smrCorrect}" 
done

# =====================================
# ======= run GTEx eQTL lookup ========
# =====================================

# set working directory
cd /slow/projects/ukb_brainage/
conda activate envs/eqtl

# set arguments
gtexDir="data/gtex"
pthresh=5E-8
(for trait in gap_gm gap_wm gap_gwm; do (
	targetDir="results/${trait}/eqtl/"
	conditionalFile="results/${trait}/conditional/conditional.cleaned.tophits.annovar.txt"
	./code/genetics/eqtl.sh "${trait}" "${targetDir}" "${conditionalFile}" "${gtexDir}" "${pthresh}"
	) &
done)
wait

# =============================================
# === run GWAS catalog screening on tophits ===
# =============================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage/
conda activate envs/default

# run catalog screening
catalogFile="data/gwas_catalog/gwas_catalog_v1.0-associations_e105_r2021-12-21.tsv"
sig=5E-8
catalogsig=5E-8
(for trait in gap_gm gap_wm gap_gwm; do (
	targetDir="results/${trait}/catalog/"
	conditionalFile="results/${trait}/conditional/conditional.cleaned.tophits.annovar.txt"
	./code/genetics/catalog/catalog.sh "${trait}" "${targetDir}" "${conditionalFile}" "${catalogFile}" "${sig}" "${catalogsig}"
	) &
done)
wait

# ==================================================
# === Calculate Polygenic Priority Scores (PoPS) ===
# ==================================================

# set working directory
cd "/slow/projects/ukb_brainage"
conda activate envs/pops

# run magma
(for trait in gap_gm gap_wm gap_gwm; do (
	mkdir -p results/${trait}/magma
	zcat results/${trait}/gwas/sumstats.txt.gz > results/${trait}/magma/sumstats.txt
	for chr in {1..22}; do (
		magma --bfile data/genetics/chr#CHR#/imp_mri_qc/bed/chr#CHR#_mri_qc \
			--batch ${chr} chr \
			--gene-annot /fast/software/pops/data/utils/magma_0kb.genes.annot \
			--pval results/${trait}/magma/sumstats.txt use=ID,P ncol=N \
			--out results/${trait}/magma/magma
			) &
	done
	wait
	magma --merge results/${trait}/magma/magma --out results/${trait}/magma/magma
	rm -f results/${trait}/magma/sumstats.txt
	rm -f results/${trait}/magma/magma.batch*
	chmod 750 results/${trait}/magma/*) &
done)
wait

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
	inputFile="results/${trait}/magma/magma.genes.out"
	outputFile="results/${trait}/magma/magma.genes.out.annot"
	conditionalFile="results/${trait}/conditional/conditional.cleaned.tophits.annovar.txt"
	./code/genetics/magma.sh "${inputFile}" "${outputFile}" "${annotationFile}" "${conditionalFile}" "${pthreshMapping}" "${glist_hg19}" "${threads}" "${window}" "${humandb}" "${refseq}" "${clumpingWindow}" 
) &
done)
wait

# calculate polygenic priority score (PoPS)
annotationFile="/fast/software/pops/data/utils/gene_annot_jun10.txt"
featurePrefix="/fast/software/pops/data/features_munged/pops_features"
nChunks=115
controlFeatures="/fast/software/pops/data/utils/control.features"
(for trait in gap_gm gap_wm gap_gwm; do (
	mkdir -p results/${trait}/pops
	python /fast/software/pops/pops.py \
		--gene_annot_path ${annotationFile} \
		--feature_mat_prefix ${featurePrefix} \
		--num_feature_chunks ${nChunks} \
		--magma_prefix results/${trait}/magma/magma \
		--control_features_path ${controlFeatures} \
		--out_prefix results/${trait}/pops/pops
	awk 'NR==1 { next } NR==FNR { symbol[$1]=$2; next }
			 FNR==1 { print "SYMBOL",$0; next }
			 $1 in symbol { print symbol[$1], $0 }' OFS='\t' ${annotationFile} results/${trait}/pops/pops.preds \
			 > results/${trait}/pops/pops.preds.symbols
	chmod 750 results/${trait}/pops/*
	) &
done)
wait

# combine results
traits="gap_gm,gap_wm,gap_gwm"
magmaFiles="results/gap_gm/magma/magma.genes.out.annot.clumped,results/gap_wm/magma/magma.genes.out.annot.clumped,results/gap_gwm/magma/magma.genes.out.annot.clumped"
popsFiles="results/gap_gm/pops/pops.preds,results/gap_wm/pops/pops.preds,results/gap_gwm/pops/pops.preds"
clumpingWindow=3000
outputFile="results/combined/pops.txt"
./code/genetics/pops.combine.R "${traits}" "${magmaFiles}" "${popsFiles}" "${clumpingWindow}" "${outputFile}"

# prioritize genes based on GWAS index variants
conditionalPthresh=5E-8
annotationFile="/fast/software/pops/data/utils/gene_annot_jun10.txt"
window=500
(for trait in gap_gm gap_wm gap_gwm; do (
	conditionalFile="results/${trait}/conditional/conditional.cleaned.txt"
	popsFile="results/${trait}/pops/pops.preds.symbols"
	outputFile="results/${trait}/pops/pops.prio.txt"
	./code/genetics/pops.prioritize.R "${conditionalFile}" "${conditionalPthresh}" "${popsFile}" "${annotationFile}" "${outputFile}" "${window}"
) &
done)
wait

# ========================================================================
# === cross-trait clumping: identify independent snp-level discoveries ===
# ========================================================================

# set working directory
cd /slow/projects/ukb_brainage/
conda activate envs/default

# settings
traitlist="gap_gm gap_wm gap_gwm"
conditionalFiles="results/gap_gm/conditional/conditional.cleaned.tophits.annovar.txt results/gap_wm/conditional/conditional.cleaned.tophits.annovar.txt results/gap_gwm/conditional/conditional.cleaned.tophits.annovar.txt"
nonsynFiles="results/gap_gm/credibleSet/credibleSet.df.annot.nonsynonymous.summary.txt results/gap_wm/credibleSet/credibleSet.df.annot.nonsynonymous.summary.txt results/gap_gwm/credibleSet/credibleSet.df.annot.nonsynonymous.summary.txt"
catalogFiles="results/gap_gm/catalog/catalog.by.locus.txt results/gap_wm/catalog/catalog.by.locus.txt results/gap_gwm/catalog/catalog.by.locus.txt"
finemapFiles="results/gap_gm/credibleSet/credibleSet.ls.genes.txt results/gap_wm/credibleSet/credibleSet.ls.genes.txt results/gap_gwm/credibleSet/credibleSet.ls.genes.txt"
smrFiles="results/gap_gm/smr/smr.summary.txt results/gap_wm/smr/smr.summary.txt results/gap_gwm/smr/smr.summary.txt"
eqtlFiles="results/gap_gm/eqtl/eqtl.summary.txt results/gap_wm/eqtl/eqtl.summary.txt results/gap_gwm/eqtl/eqtl.summary.txt"
popsFiles="results/gap_gm/pops/pops.prio.txt results/gap_wm/pops/pops.prio.txt results/gap_gwm/pops/pops.prio.txt"
targetDir="results/combined"
geneticsDir="data/genetics"
LDsample="data/gap_gm/gap_gm.txt"
pthresh=1E-6
./code/genetics/snplevel.sh "${traitlist}" "${conditionalFiles}" "${nonsynFiles}" "${catalogFiles}" "${finemapFiles}" "${smrFiles}" "${eqtlFiles}" "${popsFiles}" "${targetDir}" "${geneticsDir}" "${LDsample}" "${pthresh}" 

# add gene prioritization score
inputFile="results/combined/snplevel.gws.txt"
cols="FINEMAP_genes nonsyn_customPthresh SMR_eqtl SMR_sqtl GTEx_singleTissue GTEx_multiTissue PoPS"
outputFile="results/combined/snplevel.gws.prio.txt"
./code/genetics/prioritize.R "${inputFile}" "${cols}" "${outputFile}"

# ===============================================================================
# === Novelty: compare own with previous discoveries and create summary table ===
# ===============================================================================

# set working directory
cd /slow/projects/ukb_brainage/
conda activate envs/default

# settings
ownDiscoveries="results/combined/snplevel.gws.prio.txt"
prevDiscoveries="data/prevDiscoveries/Jonsson_2019_munged.txt data/prevDiscoveries/Ning_2020_munged.txt data/prevDiscoveries/Smith_2020_munged.txt data/prevDiscoveries/Leonardsen_2023_munged.txt" # data/prevDiscoveries/Kaufmann_2019_munged.txt
targetDir="results/combined"
geneticsDir="data/genetics"
LDsample="data/gap_gm/gap_gm.txt"

# run analysis
./code/genetics/snplevel.novelty.sh "${ownDiscoveries}" "${prevDiscoveries}" "${targetDir}" "${geneticsDir}" "${LDsample}"

# ================================
# === run replication analyses ===
# ================================

# set working directory
cd "/slow/projects/ukb_brainage"

# calculate expected number of successfull replications (accounting for winner's curse bias)
conda activate envs/power
discoveryFile="results/combined/snplevel.novelty.txt"
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
replicationN=7259
replicationP=0.10 # set 0.1 for one-tailed p = 0.05
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
covs="sex,age,age2,ac1,ac2,ac3,TIV,array,PanC{1..3}"
for trait in gap_gm gap_wm gap_gwm; do
	targetDir="results/${trait}/replicate"
	code/genetics/replicate.sh "${trait}" "${targetDir}" "${chrFileHandle}" "${data}" "${covs}" "${ancestries}" "${ancestriesCol}" "${snpQC}" "${impMFI}" "${maf}" "${threads}" 
done

ancestries="EUR"
covs="sex,age,age2,ac1,ac2,ac3,TIV,array,PanC{1..10}"
for trait in gap_gm gap_wm gap_gwm; do
	targetDir="results/${trait}/replicate"
	code/genetics/replicate.sh "${trait}" "${targetDir}" "${chrFileHandle}" "${data}" "${covs}" "${ancestries}" "${ancestriesCol}" "${snpQC}" "${impMFI}" "${maf}" "${threads}" 
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
done)
wait

# create qq-plots
pCol="P" # column containing p-values
nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5 
prune=FALSE # ommit overplotting of variants with pval > 0.01
drawCI=TRUE # draw confidence interval?
drawLambda=TRUE # draw LambdaGC?
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
done)
wait
echo completed

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
done)
wait

# run METAL
cohorts="AFR,AMR,CSA,EAS,EUR,LIFE,MID" # set output file handler & make sure that the order matches the list of sumstats (below)
type="ivweight" # type="nweight"
idCol="ID"
carryOver="CHR,BP"
leaveOneOut=0
(for trait in gap_gm gap_wm gap_gwm; do (
	targetDir="results/${trait}/replicate"
	sumstats=$(echo $(for i in $(echo $cohorts | sed 's/,/ /g'); do echo ${targetDir}/${i}/sumstats.txt.gz; done) | sed 's/ /,/g')
	code/genetics/metal.sh ${trait} ${targetDir} ${cohorts} ${sumstats} ${type} ${idCol} ${carryOver} ${leaveOneOut}
	) &
done)
wait

# quality-check METAL results
nCriterion=TRUE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5
hetP=1E-6 # heterogeneity p threshold for excluding SNPs
(for trait in gap_gm gap_wm gap_gwm; do (
	sumstats="results/${trait}/replicate/metal.ivweight.gz"
	outFile="results/${trait}/replicate/metal.ivweight.qc"
	Rscript code/genetics/metal.qc.R ${sumstats} ${outFile} ${nCriterion} ${hetP}
	) &
done)
wait

# run conditional analysis and do positional annotation (!! do not run in parallel across traits)
conda activate envs/default
subsFile="data/genetics/chr1/imp_mri_qc_EUR/chr1_mri_qc.psam"
chrFileHandle="data/genetics/chr\$i/imp_mri_qc_EUR/bed/chr\${i}_mri_qc"
cols="ID,A1,A2,A1_FREQ,BETA,SE,P,N,CHR,BP,A1_FREQ_SE,A1_FREQ_MIN,A1_FREQ_MAX,DIRECTION,HET_ISQ,HetChiSq,HET_DF,HET_P" # first 10 columns have the following order ID,A1,A2,A1_FREQ,BETA,SE,P,N,CHR,BP
pthresh=1E-6
threads=100
humandb="data/annovar/humandb"
refseq="data/refseq/GRCh37_latest_genomic.edit.gff.gz"
for trait in gap_gm gap_wm gap_gwm; do
	targetDir="results/${trait}/replicate/conditional"
	sumstats="results/${trait}/replicate/metal.ivweight.qc.gz"
	sumstatsPLINK="results/${trait}/replicate/EUR/sumstats_plink.txt.gz"
	./code/genetics/conditional.sh "${subsFile}" "${targetDir}" "${chrFileHandle}" "${sumstats}" "${cols}" "${sumstatsPLINK}" "${pthresh}" "${threads}" "${humandb}" "${refseq}"
done

# create GWAMA manhattan plots
nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5 
annotationThresh=5E-8 # pthresh=5E-8
sig=5E-8 # sig=5E-8
yend=12 # upper y axis limit
ysteps=4 # y axis breaks
width=11 # plot width (10.67 inch)
height=2.5 # plot height (4 inch)
preview=FALSE # only 1% of all variations is plotted
(for trait in gap_gm gap_wm gap_gwm; do (
	targetDir="results/${trait}/replicate/"
	sumstats="results/${trait}/replicate/metal.ivweight.qc.gz"
	conditionalFile="results/${trait}/replicate/conditional/conditional.cleaned.tophits.annovar.txt"
	Rscript code/genetics/manhattan.R "${trait}" "${targetDir}" "${sumstats}" "${nCriterion}" "${conditionalFile}" "${annotationThresh}" "${sig}" "${yend}" "${ysteps}" "${width}" "${height}" "${preview}"
	) &
done)
wait

# create GWAMA qq-plots
pCol="P" # column containing p-values
nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5 
prune=FALSE # ommit overplotting of variants with pval > 0.01
drawCI=FALSE # draw confidence interval?
drawLambda=TRUE # draw LambdaGC?
xend=8 # upper x axis limit
xsteps=2 # x axis breaks
yend=12 # upper y axis limit
ysteps=4 # y axis breaks
width=3 # plot width (3 inch)
height=2.5 # plot height (4 inch)
(for trait in gap_gm gap_wm gap_gwm; do (
	targetDir="results/${trait}/replicate/"
	sumstats="results/${trait}/replicate/metal.ivweight.qc.gz"
	Rscript code/genetics/qqplot.R "${trait}" "${targetDir}" "${sumstats}" "${pCol}" "${nCriterion}" "${prune}" "${drawCI}" "${drawLambda}" "${xend}" "${xsteps}" "${yend}" "${ysteps}" "${width}" "${height}"
	) &
done)
wait

# combine GWAMA manhattan and qq plots
traitList="gap_gm,gap_wm,gap_gwm"
plotTitles="Grey_matter,White_matter,Grey_and_white_matter" # underscore "_" will be replaced by white space
width=14
height=8
manhattanPlots=$(echo $(for trait in $(echo $traitList | sed 's/,/ /g'); do echo results/${trait}/replicate/manhattan.png; done) | sed 's/ /,/g' )
qqPlots=$(echo $(for trait in $(echo $traitList | sed 's/,/ /g'); do echo results/${trait}/replicate/qqplot.png; done) | sed 's/ /,/g' )
outputFile="results/combined/replicate.qqplot.manhattan.png"
Rscript code/genetics/manhattan.qq.combine.R "${traitList}" "${plotTitles}" "${manhattanPlots}" "${qqPlots}" "${outputFile}" "${width}" "${height}"

# get replication results of discovery index variants
replicationNcovs=8 # sex,age,age2,TIV,PC1-4 (18 in case of ac1-3,PC1-10,array)
(for trait in gap_gm gap_wm gap_gwm; do (
	discoveryFile="results/${trait}/conditional/conditional.cleaned.tophits.annovar.txt"
	replicationFile="results/${trait}/replicate/metal.ivweight.qc.gz"
	outFile="results/${trait}/replicate/replicateDiscoveries"
	Rscript code/genetics/replicateDiscoveries.R "${discoveryFile}" "${replicationFile}" "${replicationNcovs}" "${outFile}"
	) &
done)
wait

# combine results of discovery index variants (locus-wise)
snplevelFile="results/combined/snplevel.txt"
replicateFiles="results/gap_gm/replicate/replicateDiscoveries.txt,results/gap_wm/replicate/replicateDiscoveries.txt,results/gap_gwm/replicate/replicateDiscoveries.txt"
outFile="results/combined/replicateDiscoveries"
traits="gap_gm,gap_wm,gap_gwm"
traitNames="grey matter,white matter,grey and white matter"
Rscript code/genetics/replicateDiscoveries.combine.R "${snplevelFile}" "${replicateFiles}" "${outFile}" "${traits}" "${traitNames}"

	# draw qq plots
	inputFile="results/combined/replicateDiscoveries.txt"
	pCol="REPLIC_P"
	criterionCols="DISCOV_TRAITNAME"
	duplicatedCol=""
	prune=FALSE
	drawCI=TRUE
	xend=2
	xsteps=1
	yend=12
	ysteps=2
	width=2.5
	height=4

	outFile=("results/combined/replicateDiscoveries.gm.png" "results/combined/replicateDiscoveries.wm.png" "results/combined/replicateDiscoveries.gwm.png")
	criterionVals=("grey matter" "white matter" "grey and white matter")
	for (( i=0; i<${#outFile[@]}; i++ )); do
		Rscript code/genetics/replicateDiscoveries.qq.R "${inputFile}" "${outFile[i]}" "${pCol}" "${criterionCols}" "${criterionVals[i]}" "${duplicatedCol}" "${prune}" "${drawCI}" "${xend}" "${xsteps}" "${yend}" "${ysteps}" "${width}" "${height}"
	done

	xend=2.2
	outFile="results/combined/replicateDiscoveries.crosstrait.png"
	criterionCols=""
	criterionVals=""
	duplicatedCol="DISCOV_LOCUS_COUNT"
	Rscript code/genetics/replicateDiscoveries.qq.R "${inputFile}" "${outFile}" "${pCol}" "${criterionCols}" "${criterionVals}" "${duplicatedCol}" "${prune}" "${drawCI}" "${xend}" "${xsteps}" "${yend}" "${ysteps}" "${width}" "${height}"

	# combine plots
	plotTitles="Grey_matter,White_matter,Grey_and_white_matter"
	plots="results/combined/replicateDiscoveries.gm.png,results/combined/replicateDiscoveries.wm.png,results/combined/replicateDiscoveries.gwm.png"
	outFile="results/combined/replicateDiscoveries.qq.png"
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

# calculate PRS
for ancestry in EUR; do
	mkdir -p data/genetics/prs/imp_mri_qc_${ancestry}/
	awk 'NR<3 { print; next } { print $1"_"$2, $1"_"$2, $3, $4 }' data/genetics/chr1/imp_mri_qc_${ancestry}/bgen/chr1_mri_qc.sample > data/genetics/prs/imp_mri_qc_${ancestry}/chr1_mri_qc.sample
	(for trait in gap_gm gap_wm gap_gwm; do (
		Rscript /fast/software/PRSice2/PRSice.R \
			--dir /fast/software/PRSice2/ \
			--prsice /fast/software/PRSice2/PRSice_linux \
			--no-regress T \
			--base "results/${trait}/gwas/sumstats.txt.gz" \
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

# run PRS association
data="data/traits/replicate.txt"
covs="sex,age,age2,ac1,ac2,ac3,TIV,array,PanC1,PanC2,PanC3,PanC4,PanC5,PanC6,PanC7,PanC8,PanC9,PanC10"
ancestry=EUR
for trait in gap_gm gap_wm gap_gwm; do
	prsFile="data/genetics/prs/imp_mri_qc_${ancestry}/${trait}.all_score"
	outFile="results/${trait}/replicate/${ancestry}/prs.assoc.txt"
	code/genetics/replicate.prs.R "${trait}" "${data}" "${prsFile}" "${covs}" "${outFile}"
done

# combine PRS associations
traits="gap_gm,gap_wm,gap_gwm"
prsFiles="results/gap_gm/replicate/EUR/prs.assoc.txt,results/gap_wm/replicate/EUR/prs.assoc.txt,results/gap_gwm/replicate/EUR/prs.assoc.txt"
outFile="results/combined/replicate.prs.txt"
code/genetics/replicate.prs.combine.R "${traits}" "${prsFiles}" "${outFile}"

# ===================================================================
# === make GWAS discovery table for main article and supplementum ===
# ===================================================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage
conda activate envs/default

# settings
noveltyFile="results/combined/snplevel.novelty.txt"
replicationFile="results/combined/replicateDiscoveries.txt"
outputFileMain="results/combined/discoveries.main.txt"
outputFileSuppl="results/combined/discoveries.suppl.txt"
traits="gap_gm,gap_wm,gap_gwm"
traitNamesMain="GM,WM,GWM"
traitNamesSuppl="grey matter,white matter,grey and white matter"
Rscript code/genetics/snplevel.novelty.tables.R "${noveltyFile}" "${replicationFile}" "${outputFileMain}" "${outputFileSuppl}" "${traits}" "${traitNamesMain}" "${traitNamesSuppl}"

# ==========================================
# === create gwas manhattan and qq plots ===
# ==========================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage/
conda activate envs/default

# create manhattan plots
nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5
annotationThresh=5E-8 # pthresh=5E-8
sig=5E-8 # sig=5E-8
yend=35 # upper y axis limit
ysteps=5 # y axis breaks
width=11 # plot width (10.67 inch)
height=4 # plot height (4 inch)
preview=FALSE # only 1% of all variations is plotted
(for trait in gap_gm gap_wm gap_gwm; do (
	targetDir="results/${trait}/manhattan/"
	sumstats="results/${trait}/gwas/sumstats.txt.gz"
	annotationFile="results/combined/snplevel.gws.prio.txt"
	cat ${annotationFile} | awk -F'\t' '{sub("ENSG00000257890","Lnc-NUAK1-1",$NF); print}' OFS='\t' > ${annotationFile}.tmp 
	Rscript code/genetics/manhattan.R "${trait}" "${targetDir}" "${sumstats}" "${nCriterion}" "${annotationFile}.tmp" "${annotationThresh}" "${sig}" "${yend}" "${ysteps}" "${width}" "${height}" "${preview}"
	rm -f ${annotationFile}.tmp
	) &
done
wait)

# create qq-plots
pCol="P" # column containing p-values
nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5
prune=TRUE # ommit overplotting of variants with pval > 0.01
drawCI=FALSE # draw confidence interval?
drawLambda=FALSE # draw LambdaGC?
xend=8 # upper x axis limit
xsteps=2 # x axis breaks
yend=35 # upper y axis limit
ysteps=5 # y axis breaks
width=3 # plot width (4 inch)
height=4 # plot height (4 inch)
(for trait in gap_gm gap_wm gap_gwm; do (
	targetDir="results/${trait}/qqplot/"
	sumstats="results/${trait}/gwas/sumstats.txt.gz"
	Rscript code/genetics/qqplot.R "${trait}" "${targetDir}" "${sumstats}" "${}" "${nCriterion}" "${prune}" "${drawCI}" "${drawLambda}" "${xend}" "${xsteps}" "${yend}" "${ysteps}" "${width}" "${height}"
	) &
done
wait)

# combine Manhattan and qq-plots
traits="gap_gm,gap_wm,gap_gwm"
plotTitles="" # annotation will be a-f if no plot titles are provided
manhattanPlots="results/gap_gm/manhattan/manhattan.png,results/gap_wm/manhattan/manhattan.png,results/gap_gwm/manhattan/manhattan.png"
qqPlots="results/gap_gm/qqplot/qqplot.png,results/gap_wm/qqplot/qqplot.png,results/gap_gwm/qqplot/qqplot.png"
outputFile="results/combined/qqplot.manhattan.png"
width=14
height=12
Rscript code/genetics/manhattan.qq.combine.R "${traits}" "${plotTitles}" "${manhattanPlots}" "${qqPlots}" "${outputFile}" "${width}" "${height}"

# =========================================
# === run gene-based analysis (fastBAT) ===
# =========================================

# set working directory and load conda environment
cd "/slow/projects/ukb_brainage"
conda activate envs/default

# run analysis
geneticsDir="data/genetics"
pthreshMapping=5E-8
glist_hg19="data/glist_hg19/glist-hg19-refseq.txt"
threads=26
window=0
humandb="data/annovar/humandb"
refseq="data/refseq/GRCh37_latest_genomic.edit.gff.gz"
clumpingWindow=3000

	# for gcta software: change chromosome strings to numeric values
	\cp ${geneticsDir}/chrX/imp_mri_qc/bed/chrX_mri_qc_numeric.bim ${geneticsDir}/chrX/imp_mri_qc/bed/chrX_mri_qc.bim 
	\cp ${geneticsDir}/chrY/imp_mri_qc/bed/chrY_mri_qc_numeric.bim ${geneticsDir}/chrY/imp_mri_qc/bed/chrY_mri_qc.bim
	\cp ${geneticsDir}/chrXY/imp_mri_qc/bed/chrXY_mri_qc_numeric_23.bim ${geneticsDir}/chrXY/imp_mri_qc/bed/chrXY_mri_qc.bim
	\cp ${geneticsDir}/chrMT/imp_mri_qc/bed/chrMT_mri_qc_numeric.bim ${geneticsDir}/chrMT/imp_mri_qc/bed/chrMT_mri_qc.bim

	# run analysis
	(for trait in gap_gm gap_wm gap_gwm; do (
		targetDir="results/${trait}/fastbat"
		sumstats="results/${trait}/gwas/sumstats.txt.gz"
		conditionalFile="results/${trait}/conditional/conditional.cleaned.tophits.annovar.txt"
		./code/genetics/fastbat.sh "${trait}" "${targetDir}" "${geneticsDir}" "${sumstats}" "${conditionalFile}" "${pthreshMapping}" "${glist_hg19}" "${threads}" "${window}" "${humandb}" "${refseq}" "${clumpingWindow}"
		) &
	done)

	# replace .bim files with original files (change numeric values to chromosome strings)
	\cp ${geneticsDir}/chrX/imp_mri_qc/bed/chrX_mri_qc_original.bim ${geneticsDir}/chrX/imp_mri_qc/bed/chrX_mri_qc.bim
	\cp ${geneticsDir}/chrY/imp_mri_qc/bed/chrY_mri_qc_original.bim ${geneticsDir}/chrY/imp_mri_qc/bed/chrY_mri_qc.bim
	\cp ${geneticsDir}/chrXY/imp_mri_qc/bed/chrXY_mri_qc_original.bim ${geneticsDir}/chrXY/imp_mri_qc/bed/chrXY_mri_qc.bim
	\cp ${geneticsDir}/chrMT/imp_mri_qc/bed/chrMT_mri_qc_original.bim ${geneticsDir}/chrMT/imp_mri_qc/bed/chrMT_mri_qc.bim

# combine results
traits="gap_gm,gap_wm,gap_gwm"
inputFiles="results/gap_gm/fastbat/fastbat.0kb.gene.fastbat.annovar.clumped,results/gap_wm/fastbat/fastbat.0kb.gene.fastbat.annovar.clumped,results/gap_gwm/fastbat/fastbat.0kb.gene.fastbat.annovar.clumped"
clumpingWindow=3000 # cross-trait clumping window size in kb
outputFile="results/combined/fastbat.txt"
./code/genetics/fastbat.combine.R "${traits}" "${inputFiles}" "${clumpingWindow}" "${outputFile}"

# create gene-based Manhattan plots
annotationThresh='bonferroni' # 'fdr' or 'bonferroni'
ylim=25 # upper y axis limit
ysteps=5 # y axis breaks
width=11 # plot width
height=4 # plot height
for trait in gap_gm gap_wm gap_gwm; do (
	targetDir="results/${trait}/fastbat/"
	fastbatResults="results/${trait}/fastbat/fastbat.0kb.gene.fastbat.annovar.clumped"
	Rscript code/genetics/fastbat.manhattan.R "${trait}" "${targetDir}" "${fastbatResults}" "${annotationThresh}" "${ylim}" "${ysteps}" "${width}" "${height}"
	) &
done
wait

# create gene-based qq-plots
pCol="Pvalue" # column containing p-values
nCriterion=FALSE # LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5
prune=TRUE # ommit overplotting of variants with pval > 0.01
drawCI=FALSE # draw confidence interval?
drawLambda=FALSE # draw LambdaGC?
xend=8 # upper x axis limit
xsteps=2 # x axis breaks
yend=25 # upper y axis limit
ysteps=5 # y axis breaks
width=3 # plot width (4 inch)
height=4 # plot height (4 inch)
(for trait in gap_gm gap_wm gap_gwm; do (
	targetDir="results/${trait}/fastbat/"
	sumstats="results/${trait}/fastbat/fastbat.0kb.gene.fastbat.annovar.clumped"
	Rscript code/genetics/qqplot.R "${trait}" "${targetDir}" "${sumstats}" "${pCol}" "${nCriterion}" "${prune}" "${drawCI}" "${drawLambda}" "${xend}" "${xsteps}" "${yend}" "${ysteps}" "${width}" "${height}"
	) &
done
wait)

# combine Manhattan and qq-plots
traits="gap_gm,gap_wm,gap_gwm"
plotTitles="" # annotation will be a-f if no plot titles are provided
manhattanPlots="results/gap_gm/fastbat/fastbat.manhattan.png,results/gap_wm/fastbat/fastbat.manhattan.png,results/gap_gwm/fastbat/fastbat.manhattan.png"
qqPlots="results/gap_gm/fastbat/qqplot.png,results/gap_wm/fastbat/qqplot.png,results/gap_gwm/fastbat/qqplot.png"
outputFile="results/combined/fastbat.qqplot.manhattan.png"
width=14
height=12
Rscript code/genetics/manhattan.qq.combine.R "${traits}" "${plotTitles}" "${manhattanPlots}" "${qqPlots}" "${outputFile}" "${width}" "${height}"

# ============================
# === run pathway analyses ===
# ============================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage
conda activate envs/gofuncr

# get genes for overrepresentation test
fastbatFile="results/combined/fastbat.txt"
fastbatGeneCol="Gene"
snplevelFile="results/combined/snplevel.novelty.txt"
refseqFile="data/refseq/GRCh37_latest_genomic.edit.gff.gz"
clumpingWindow=3000 # clumping window in kb
for trait in gap_gm gap_wm gap_gwm; do
	mkdir "results/${trait}/gofuncr/"
	fastbatpCol="${trait}_FDR" # column that indicates gene significance
	outFile="results/${trait}/gofuncr/gofuncr.over.input.txt" # output file
	Rscript code/genetics/gofuncr.over.input.R "${trait}" "${fastbatFile}" "${fastbatpCol}" "${fastbatGeneCol}" "${snplevelFile}" "${refseqFile}" "${clumpingWindow}" "${outFile}"
done

# run gene set enrichment test
fastbatFile="results/combined/fastbat.txt"
fastbatGeneCol="Gene"
fastbatClumpCol="LOCUS_COUNT" # column that indicates locus of gene
crossFWER=TRUE # calculate FWER across the three ontologies (instead of ontology-wise calculation)
(for trait in gap_gm gap_wm gap_gwm; do (
	fastbatpCol="${trait}_Pvalue"
	outFile="results/${trait}/gofuncr/gofuncr.gsea"
	Rscript code/genetics/gofuncr.gsea.R "${fastbatFile}" "${fastbatpCol}" "${fastbatGeneCol}" "${fastbatClumpCol}" "${crossFWER}" "${outFile}"
	) &
done
wait)

# combine results of gene set enrichment tests
traits="gap_gm,gap_wm,gap_gwm"
overFiles="results/gap_gm/gofuncr/gofuncr.gsea.results.txt,results/gap_wm/gofuncr/gofuncr.gsea.results.txt,results/gap_gwm/gofuncr/gofuncr.gsea.results.txt"
outFile="results/combined/gofuncr.gsea.txt"
Rscript code/genetics/gofuncr.gsea.combine.R "${traits}" "${overFiles}" "${outFile}"

# ========================================
# === run genetic correlation analysis ===
# ========================================

# set working directory and load conda environment
cd "/slow/projects/ukb_brainage"
conda activate envs/ldsc

# calculate rg with selected traits
threads=100 # set number of parallel workers
LDchr="data/ldsc/resources/eur_w_ld_chr/"
for trait in gap_gm gap_wm gap_gwm; do
	targetDir="results/${trait}/rgSelection"
	sumstatsX="results/${trait}/ldsc/ldsc.${trait}.sumstats.gz"
	sumstatsY=$(find data/sumstats/munged/*sumstats.gz)
	xprefix=""
	yprefix=""
	./code/genetics/rg.sh "${targetDir}" "${sumstatsX}" "${sumstatsY}" "${xprefix}" "${yprefix}" "${LDchr}" "${threads}"
done

	# combine rg results across brain age gap variables
	conda activate envs/default
	traits="gap_gm,gap_wm,gap_gwm"
	inputFiles=$(echo $(for i in $(echo ${traits} | sed 's/,/ /g'); do echo results/${i}/rgSelection/rg.results; done) | sed 's/ /,/g')
	outputFile="results/combined/rgSelection.txt"
	./code/genetics/rg.combine.R "${traits}" "${inputFiles}" "${outputFile}"

	# replace trait names by trait labels
	inputFile="results/combined/rgSelection.txt"
	outFile="results/combined/rgSelection.labels.txt"
	matchVar="p2"
	./code/genetics/rg.addLabels.R "${inputFile}" "${outFile}" "${matchVar}" 

	# plot rg results
	corrTable="results/combined/rgSelection.txt"
	traits="gap_gm,gap_wm,gap_gwm"
	traitLabels="Grey_matter,White_matter,Grey_&_White"
	matchVar="p2"
	outFile="results/combined/rgSelection.png"
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
    targetDir="results/${trait}/rg"
    sumstats="results/${trait}/ldsc/ldsc.${trait}.sumstats.gz"
    ./code/genetics/rgNeale.sh "${trait}" "${targetDir}" "${sumstats}" "${nealeDir}" "${nealeManifest}" "${showcaseCategories}" "${LDchr}" "${threads}" # taskset -c 0-24 
done

	# combine rg results of Neale traits
	conda activate envs/default
	traits="gap_gm,gap_wm,gap_gwm"
	inputFiles=$(echo $(for i in $(echo ${traits} | sed 's/,/ /g'); do echo results/${i}/rg/rg.results; done) | sed 's/ /,/g')
	outputFile="results/combined/rgNeale.txt"
	./code/genetics/rgNeale.combine.R "${traits}" "${inputFiles}" "${outputFile}"

	# draw volcano and forest plot of rg Neale results
	inputFile="results/combined/rgNeale.txt"
	outFile="results/combined/rgNeale.volcano.forest.png"
	rgCol="gap_gm_rg"
	seCol="gap_gm_se"
	pCol="gap_gm_p"
	multipleTesting='fdr' # 'fdr' or 'bonferroni' or 'both'
	ylim=10 # upper y axis limit
	ysteps=2 # y axis breaks
	xlim=0.5 # upper y axis limit 
	xsteps=0.25 # y axis breaks 
	width=7.56 # plot width
	height=8.58 # 3.85 # plot height
	./code/genetics/rgNeale.volcano.forest.R "${inputFile}" "${outFile}" "${rgCol}" "${seCol}" "${pCol}" "${multipleTesting}" "${ylim}" "${ysteps}" "${xlim}" "${xsteps}" "${width}" "${height}"

	# combine plot of selected traits and Neale traits
	rgSelection="results/combined/rgSelection.png"
	rgNeale="results/combined/rgNeale.volcano.forest.png"
	outFile="results/combined/rgCombined.png"
	./code/genetics/rg.plotCombine.R "${rgSelection}" "${rgNeale}" "${outFile}"

	# draw qq plot of rg Neale results
	inputFile="results/combined/rgNeale.txt"
	outFile="results/combined/rgNeale.qq.png"
	rgCol="gap_gm_rg"
	pCol="gap_gm_p"
	ylim=9 # upper y axis limit
	ysteps=3 # y axis breaks
	xlim=3 # upper y axis limit;
	xsteps=1 # y axis breaks
	width=4.79 # plot width
	height=4.57 # plot height
	./code/genetics/rgNeale.qq.R "${inputFile}" "${outFile}" "${rgCol}" "${pCol}" "${ylim}" "${ysteps}" "${xlim}" "${xsteps}" "${width}" "${height}"

# prepare summary statistics for LDhub
snplist="data/ldsc/resources/w_hm3.noMHC.snplist"
for trait in gap_gm gap_wm gap_gwm; do
    targetDir="results/${trait}/ldhub"
    sumstats="results/${trait}/gwas/sumstats.txt.gz"
    ./code/genetics/ldhub.sh "${trait}" "${targetDir}" "${sumstats}" "${snplist}"
done

# ===================================
# === run Mendelian Randomization ===
# ===================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage
conda activate envs/default

# settings
sumstatsCols="ID,A1,A2,A1_FREQ,BETA,SE,P,N"
mrFiles=$(echo $(ls data/sumstats/mr/mr_*.gz) | sed 's/ /,/g')
mrLabels=$(echo $(for i in $(echo $mrFiles | sed 's/,/ /g'); do echo $i | sed 's%.*/%%g' | sed 's%\..*%%g'; done) | sed 's/ /,/g')
chrFilehandler='data/genetics/chr${i}/imp_mri_qc/bed/chr${i}_mri_qc'

(for trait in gap_gm gap_wm gap_gwm; do (
	targetDir="results/${trait}/gsmr"
	sumstats="results/${trait}/gwas/sumstats.txt.gz"
	sampleFile=data/${trait}/${trait}.txt
	./code/genetics/gsmr.sh "${trait}" "${targetDir}" "${sumstats}" "${sumstatsCols}" "${sampleFile}" "${mrFiles}" "${mrLabels}" "${chrFilehandler}"
	) &
done)
wait

# combine mr results
traits="gap_gm,gap_wm,gap_gwm"
inputFiles=$(echo $(for i in $(echo ${traits} | sed 's/,/ /g'); do echo results/${i}/gsmr/gsmr.gsmr; done) | sed 's/ /,/g')
outputFile="results/combined/gsmr.txt"
direction="outcome"
mInstruments=-1
./code/genetics/gsmr.combine.R "${traits}" "${inputFiles}" "${outputFile}" "${direction}" "${mInstruments}"

# replace trait names by trait labels
inputFile="results/combined/gsmr.txt"
outFile="results/combined/gsmr.labels.txt"
matchVar="trait"
./code/genetics/gsmr.addLabels.R "${inputFile}" "${outFile}" "${matchVar}" 

# plot mr results
conda activate envs/default
traits=(gap_gm gap_wm gap_gwm)
outcome_labels=("Grey matter brain age gap" "White matter brain age gap" "Grey and white matter brain age gap")
for (( i=0; i<${#traits[@]}; i++ )); do
	gsmrFile="results/${traits[i]}/gsmr/gsmr.eff_plot.gz"
	exposure_levels="mr_bmi_locke_2015,mr_whr_shungin_2015,mr_dbp_noUKB_evangelou_2018,mr_sbp_noUKB_evangelou_2018,mr_pp_noUKB_evangelou_2018,mr_ldlc_willer_2013,mr_hdlc_willer_2013,mr_trigl_willer_2013,mr_cad_nikpay_2015,mr_t2d_scott_2017,mr_scz_trubetskoy_2022,mr_eduPruned_okbay_2016"
	exposure_labels="Body-Mass-Index_(BMI),Waist-Hip-Ratio_(BMI adjusted),Diastolic Blood Pressure,Systolic Blood Pressure,Pulse Pressure,LDL-c,HDL-c,Triglyceride,Coronary_Artery_Disease,Type-2-Diabetes,Schizophrenia,Educational_Attainment"
	outcome_level=${traits[i]}
	outcome_label=${outcome_labels[i]}
	targetDir="results/${traits[i]}/gsmr/"
	./code/genetics/gsmr.plot.R "${gsmrFile}" "${exposure_levels}" "${exposure_labels}" "${outcome_level}" "${outcome_label}" "${targetDir}"
done

# ===============================================================
# === run genetic effect size distribution analysis (GENESIS) ===
# ===============================================================

# set working directory and load conda environment
cd /slow/projects/ukb_brainage
conda activate envs/genesis

# run analysis for traits of interest
ncores=10
compFit3=TRUE
(for trait in gap_gm gap_wm gap_gwm; do (
	mkdir -p results/${trait}/genesis/
	sumstats="results/${trait}/gwas/sumstats.txt.gz"
	outFile="results/${trait}/genesis/genesis"
	Rscript code/genetics/genesis.R "${trait}" "${sumstats}" "${ncores}" "${compFit3}" "${outFile}"
	) &
done)
wait

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
	genesisFit="results/${trait}/genesis/genesis.fit3.Rds"
	outFile="results/${trait}/genesis/genesis"
	Rscript code/genetics/genesis.plot.R "${trait}" "${traitLabels}" "${genesisFit}" "${reftraits}" "${reftraitLabels}" "${refgenesisFit}" "${outFile}"
done

# create stats table
traits="gap_gm,gap_wm,gap_gwm,height,neur"
traitLabels="Grey_matter_brain_age_gap,White_matter_brain_age_gap,Grey_and_white_matter_brain_age_gap,Height_(Wood_et_al._2014),Neuroticism_(Baselmans_et_al._2019)"
genesisFit="results/gap_gm/genesis/genesis.fit3.Rds,results/gap_wm/genesis/genesis.fit3.Rds,results/gap_gwm/genesis/genesis.fit3.Rds,data/sumstats/genesis/mr_height_wood_2014/genesis.fit3.Rds,data/sumstats/genesis/04_neur_baselmans_2019/genesis.fit3.Rds"
outFile="results/combined/genesis.stats"
Rscript code/genetics/genesis.stats.R "${traits}" "${traitLabels}" "${genesisFit}" "${outFile}"

# =====================================
# === combine results across traits ===
# =====================================

# set working director and conda environment
cd /slow/projects/ukb_brainage/
conda activate envs/default

# combine finemap results across traits
snplevelFile="results/combined/snplevel.gws.txt"
finemapFiles="results/gap_gm/credibleSet/credibleSet.df.annot.txt,results/gap_wm/credibleSet/credibleSet.df.annot.txt,results/gap_gwm/credibleSet/credibleSet.df.annot.txt"
outputFile="results/combined/finemap.suppl.txt"
traits="gap_gm,gap_wm,gap_gwm"
traitNames="grey matter,white matter,grey and white matter"
Rscript code/genetics/finemap.combine.R "${snplevelFile}" "${finemapFiles}" "${outputFile}" "${traits}" "${traitNames}"

# combine exonic results for supplementum
snplevelFile="results/combined/snplevel.gws.txt"
nsynFiles="results/gap_gm/credibleSet/credibleSet.df.annot.nonsynonymous.txt,results/gap_wm/credibleSet/credibleSet.df.annot.nonsynonymous.txt,results/gap_gwm/credibleSet/credibleSet.df.annot.nonsynonymous.txt"
outputFile="results/combined/nonsynonymous.suppl.txt"
traits="gap_gm,gap_wm,gap_gwm"
traitNames="grey matter,white matter,grey and white matter"
code/genetics/annotation.exonic.combine.R "${snplevelFile}" "${nsynFiles}" "${outputFile}" "${traits}" "${traitNames}" 

# combine SMR eqtl results across traits
snplevelFile="results/combined/snplevel.gws.txt"
smrFiles="results/gap_gm/smr/smr.eqtl.filtered.assigned.txt,results/gap_wm/smr/smr.eqtl.filtered.assigned.txt,results/gap_gwm/smr/smr.eqtl.filtered.assigned.txt"
outputFile="results/combined/smr.eqtl.suppl.txt"
traits="gap_gm,gap_wm,gap_gwm"
traitNames="grey matter,white matter,grey and white matter"
Rscript code/combine.smr.R "${snplevelFile}" "${smrFiles}" "${outputFile}" "${traits}" "${traitNames}"

# combine SMR sqtl results across traits
snplevelFile="results/combined/snplevel.gws.txt"
smrFiles="results/gap_gm/smr/smr.sqtl.filtered.assigned.txt,results/gap_wm/smr/smr.sqtl.filtered.assigned.txt,results/gap_gwm/smr/smr.sqtl.filtered.assigned.txt"
outputFile="results/combined/smr.sqtl.suppl.txt"
traits="gap_gm,gap_wm,gap_gwm"
traitNames="grey matter,white matter,grey and white matter"
Rscript code/combine.smr.R "${snplevelFile}" "${smrFiles}" "${outputFile}" "${traits}" "${traitNames}"

# combine eQTL single-tissue results across traits
snplevelFile="results/combined/snplevel.gws.txt"
eqtlFiles="results/gap_gm/eqtl/eqtl.singleTissue.txt,results/gap_wm/eqtl/eqtl.singleTissue.txt,results/gap_gwm/eqtl/eqtl.singleTissue.txt"
outputFile="results/combined/gtex.singletissue.txt"
traits="gap_gm,gap_wm,gap_gwm"
traitNames="grey matter,white matter,grey and white matter"
Rscript code/genetics/eqtl.single.combine.R "${snplevelFile}" "${eqtlFiles}" "${outputFile}" "${traits}" "${traitNames}"

# combine eQTL multi-tissue results across traits
snplevelFile="results/combined/snplevel.gws.txt"
eqtlFiles="results/gap_gm/eqtl/eqtl.multiTissue.txt,results/gap_wm/eqtl/eqtl.multiTissue.txt,results/gap_gwm/eqtl/eqtl.multiTissue.txt"
outputFile="results/combined/gtex.multitissue.txt"
traits="gap_gm,gap_wm,gap_gwm"
traitNames="grey matter,white matter,grey and white matter"
Rscript code/genetics/eqtl.multi.combine.R "${snplevelFile}" "${eqtlFiles}" "${outputFile}" "${traits}" "${traitNames}"
