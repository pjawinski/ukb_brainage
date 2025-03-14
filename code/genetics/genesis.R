#!/usr/bin/env Rscript

# =====================================================
# === run genetic effect size distribution analysis ===
# =====================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=7) {
  stop(paste0('expected 7 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}                                                                                                                                                    
            
# set arguments
trait = args[1] # trait="gap_gm"
sumstats = args[2] # sumstats="results/gap_gm/gwas/sumstats.txt.gz"
colsIn = args[3] # cols="ID,Z,N"
colsOut = args[4] # cols="SNP,Z,N"
ncores = as.numeric(args[5]) # ncores=100
compFit3 = args[6] # compFit3=FALSE
outFile = args[7] # outFile="results/gap_gm/genesis/genesis"

logInfo = paste0('\n--- run genetic effect size distribution analysis ---',
               '\ntrait: ', trait,
               '\nsumstats: ', sumstats,
               '\ncolsIn: ', colsIn,
               '\ncolsOut: ', colsOut,
               '\nncores: ', ncores,
               '\ncompFit3: ', compFit3,
               '\noutFile: ', outFile,'\n')
message(logInfo)

# check if package genesis is available and conda environment name is 'genesis'
is_genesis_available <- suppressWarnings(require('GENESIS'))
env = system('echo $CONDA_DEFAULT_ENV', intern = T)
if (!is_genesis_available ) {
  message(' - GENESIS not available. Checking conda environment.')
  if (env != 'genesis') {
    message(' - conda environment name is not genesis. Please load environment genesis.'); stop()
  } else {
    message(' - conda environment name is genesis. Downloading missing package genesis.')
    options(download.file.method = "curl")
    devtools::install_github('yandorazhang/GENESIS', upgrade = 'never')
  }
}

# attach required packages
for (pkg in c('dplyr','stringr','foreach','iterators','parallel','doParallel','data.table','GENESIS')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# transform variables
colsIn = str_split(colsIn, ',')[[1]]
colsOut = str_split(colsOut, ',')[[1]]

# import data
message(paste0(' - importing data'))
if (stringr::str_sub(sumstats,-3,-1) == '.gz') {
  df = data.frame(data.table::fread(cmd=sprintf("gzip -dc %s", sumstats), tmpdir = getwd(), select = colsIn, header=T, sep = '\t', stringsAsFactors=FALSE)) 
} else {
  df = data.frame(data.table::fread(file = sumstats, tmpdir = getwd(), select = colsIn, header=T, sep = '\t', stringsAsFactors=FALSE))
}

# The input GWAS summary statistics should contain 3 columns:
# SNP rsID,
# original z-statistics got from GWAS study (uncorrected by genomic control factor),
# effective sample size of GWAS study, which can vary for different SNPs; for disease traits, the sample size should be effective sample size, i.e., (# of cases)*(# of controls)/(total # of cases & controls).
names(df) = colsOut

  # calculate Z if necessary
  if (sum(c('BETA','SE') %in% names(df)) == 2) {
    df$Z = df$BETA/df$SE
    df = df[,c('SNP','Z','N')]
  }

  # calculate effective sample size (as defined by Zhang et al. as n1*n0/(n1+n0), which basically is 2/(1/n1+1/n0)/2, i.e., effective sample size per group divided by 2)
  if (sum(c('BETA','SE') %in% names(df)) == 2) {
    df$Z = df$BETA/df$SE
    df = df[,c('SNP','Z','N')]
  }

# The input GWAS summary statistics are strongly recommended to do filtering before fitting to the model:
# If sample size is different for different SNPs, remove SNPs with an effective sample size less than 0.67 times the 90th percentile of sample size.
# Filter SNPs to Hapmap3 SNPs which are not in the major histocompatibility complex (MHC) region. For Hapmap3 SNP list without MHC region, type data(w_hm3.noMHC.snplist) in R.
# Remove SNPs with extremely large effect sizes (z^2 > 80).
message(' - preprocessing sumstats')
df.prep = preprocessing(df, filter = TRUE)
df.prep = df.prep[,c('SNP','Z','N')]

# ssh time outs when script is run with Rscript: use ForkCluster instead of PSOCKcluster
body(genesis)[8] <- quote((cl <- makeForkCluster(cores))())

# fit the 2-component model
set.seed(4924)
message(' - fitting the 2-component model')
fit2 = genesis(df.prep, filter=F, modelcomponents = 2, cores = ncores, LDcutoff = 0.1, LDwindow = 1, c0=10, startingpic=0.005, herit.liability = TRUE, qqplot = FALSE)

# get starting values of 3-component model from 2-component model
est = fit2$estimates$`Parameter (pic, sigmasq, a) estimates` # the model parameter estimates
v = fit2$estimates$`Covariance matrix of parameter estimates` # the covariance matrix of model parameter estimtaes
starting = rep(0,5)                                                                                                                                                                                          
starting[1] = est[1]                                                                                                                                                                                         
starting[2] = 1/9                                                                                                                                                                                            
starting[3] = est[2]*5                                                                                                                                                                                       
starting[4] = starting[3]/10                                                                                                                                                                                 
starting[5] = est[3]    

# fit the 3-component model
if (compFit3 == TRUE) {
  message(' - fitting the 3-component model.')
  fit3 = genesis(df.prep, filter=F, modelcomponents = 3, cores = ncores, LDcutoff = 0.1, LDwindow = 1, c0=10, starting = starting) # summaryGWASdata.save=TRUE, qqplotdata.save=TRUE)
}

# create outtable for fit2 model
tblfit2 = data.frame(model = 'fit2',
  h2 = fit2$estimates$`Total heritability in log-odds-ratio scale (sd)`,
  sSnps = fit2$estimates$`Number of sSNPs (sd)`,
  logLikelihood = fit2$estimates$`Composite log-likelihood of fitted model`,
  bic = fit2$estimates$`Model selection related`$BIC)
tblfit2$h2_se = tblfit2$h2 %>% gsub(pattern = '.*[(]', replacement = '') %>% gsub(pattern = '[)].*', replacement = '')
tblfit2$h2 = tblfit2$h2 %>% gsub(pattern = '[ ].*', replacement = '')
tblfit2$sSnps_se = tblfit2$sSnps %>% gsub(pattern = '.*[(]', replacement = '') %>% gsub(pattern = '[)].*', replacement = '')
tblfit2$sSnps = tblfit2$sSnps %>% gsub(pattern = '[ ].*', replacement = '')
tblfit2 = tblfit2[,c('model','h2','h2_se','sSnps','sSnps_se','logLikelihood','bic')]

# create outtable for fit3 model
if (compFit3 == TRUE) {
  tblfit3 = data.frame(model = 'fit3',
    h2 = fit3$estimates$`Total heritability in log-odds-ratio scale (sd)`,
    sSnps = fit3$estimates$`Number of sSNPs (sd)`,
    logLikelihood = fit3$estimates$`Composite log-likelihood of fitted model`,
    bic = fit3$estimates$`Model selection related`$BIC)
  tblfit3$h2_se = tblfit3$h2 %>% gsub(pattern = '.*[(]', replacement = '') %>% gsub(pattern = '[)].*', replacement = '')
  tblfit3$h2 = tblfit3$h2 %>% gsub(pattern = '[ ].*', replacement = '')
  tblfit3$sSnps_se = tblfit3$sSnps %>% gsub(pattern = '.*[(]', replacement = '') %>% gsub(pattern = '[)].*', replacement = '')
  tblfit3$sSnps = tblfit3$sSnps %>% gsub(pattern = '[ ].*', replacement = '')
  tblfit3 = tblfit3[,c('model','h2','h2_se','sSnps','sSnps_se','logLikelihood','bic')]
}

# combine tables
if (compFit3 == TRUE) {
  tblfit = rbind(tblfit2,tblfit3)
  } else {
  tblfit = tblfit2
}

# write results
message(sprintf(' - writing %s.stats.txt',outFile))
write.table(tblfit, file = sprintf('%s.stats.txt',outFile), row.names = F, quote = F, sep = '\t', na = "-")
message(sprintf(' - writing %s.fit2.Rds',outFile))
saveRDS(file = sprintf('%s.fit2.Rds',outFile),fit2)
if (compFit3 == TRUE) { 
  message(sprintf(' - writing %s.fit3.Rds',outFile))
  saveRDS(file = sprintf('%s.fit3.Rds',outFile),fit3)
}
system(sprintf('chmod 770 %s*', outFile))
message('-- Completed: run genetic effect size distribution analysis ---')
