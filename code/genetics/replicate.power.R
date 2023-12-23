#!/usr/bin/env Rscript

# =================================
# === Compute replication power ===
# =================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=14) {
  stop(paste0('expected 14 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
discoveryFile=args[1] # discoveryFile="results/combined/snplevel.novelty.txt" # ="results/combined/discoveries.suppl.txt"
traits=args[2] # traits="gap_gm,gap_wm,gap_gwm"
traitNames=args[3] # traitNames="grey matter,white matter,grey and white matter"
colTrait=args[4] # colTrait="TRAIT"
colBETA=args[5] # colBETA="BETA"
colSE=args[6] # colSE="SE"
colN=args[7] # colN="N"
colFREQ=args[8] # colFREQ="A1_FREQ"
cols2keep=args[9] # cols2keep="LOCUS_COUNT,TRAIT,CHR,BP,ID,NEAREST_GENE,A1,A2,A1_FREQ,BETA,SE,Z,P,N"
filterCol=args[10] # filterCol="DISCOV_COUNT"
filterVal=as.numeric(args[11]) # filterVal=1
replicationN=as.numeric(args[12]) # replicationN=7259
replicationP=as.numeric(args[13]) # replicationP=0.10 # set 0.1 for one-tailed p = 0.05
outFile=args[14] # outFile="results/combined/replication.power"

logInfo = paste0('\n--- Compute replication power ---',
               '\ndiscoveryFile: ', discoveryFile,
               '\ntraits: ', traits,
               '\ntraitNames: ', traitNames,
               '\ncolTrait: ', colTrait,
               '\ncolBETA: ', colBETA,
               '\ncolSE: ', colSE,
               '\ncolN: ', colN,
               '\ncolFREQ: ', colFREQ,
               '\ncols2keep: ', cols2keep,
               '\nfilterCol: ', filterCol,
               '\nfilterVal: ', filterVal,
               '\nreplicationN: ', replicationN,
               '\nreplicationP: ', replicationP,
               '\noutFile: ', outFile, '\n')
message(logInfo)

# check if package gwas.winners.curse is available and conda environment name is 'gwaspower'
is_winnercurse_available <- suppressWarnings(require('gwas.winners.curse'))
env = system('echo $CONDA_DEFAULT_ENV', intern = T)
if (!is_winnercurse_available ) {
  message(' - package gwas.winners.curse not available. Checking conda environment.')
  if (env != 'power') {
    message(' - conda environment name is not power. Please load environment gwaspower.'); stop()
  } else {
    message(' - conda environment name is power. Downloading missing package gwas.winners.curse.')
    options(download.file.method = "curl")
    devtools::install_github('lightning-auriga/gwas-winners-curse', ref = 'default', upgrade = 'never')
  }
}

# attach packages to current R session and load required functions
for (pkg in c('dplyr', 'stringr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }
source('code/genetics/power_calc_functions.R')

# transform variables
cols2keep = str_split(cols2keep, ',')[[1]]
traits = str_split(traits, ',')[[1]]
traitNames = str_split(traitNames, ',')[[1]]

# load data and apply filter
message('Loading Files.')
df = read.delim(discoveryFile, sep = '\t', header = T)
df = df[df[[filterCol]]==filterVal,]

# adding traitlabels
message('Adding traitNames.')
df[[colTrait]] = factor(df[[colTrait]], levels = traits, labels = traitNames)

# standardise BETA and SE
message('Standardizing BETA and SE.')
df$BETA_STD = (df[[colBETA]]/df[[colSE]])/sqrt(2*df[[colFREQ]]*(1-df[[colFREQ]])*(df[[colN]]+(df[[colBETA]]/df[[colSE]])^2))
df$SE_STD = 1/sqrt(2*df[[colFREQ]]*(1-df[[colFREQ]])*(df[[colN]]+(df[[colBETA]]/df[[colSE]])^2))

# correct for winner's curse
message('Correcting beta estimates for winners curse.')
input = data.frame(DISCOV_BETA = df$BETA_STD, DISCOV_SE = df$SE_STD, DISCOV_N = df[[colN]], DISCOV_FREQ = df[[colFREQ]], REPL_BETA = rep(NA,nrow(df)), REPL_SE = rep(NA,nrow(df)), REPL_N = rep(NA,nrow(df)), REPL_FREQ = rep(NA,nrow(df)))
write.table(input, 'tmp.wcinput.txt', sep = '\t', row.names = F, col.names = T)
gwas.winners.curse::correct.winners.curse('tmp.wcinput.txt','tmp.wcoutput.txt',0,5e-8,header = TRUE,sep = "\t")
wcurse = read.table('tmp.wcoutput.txt', sep = '\t', header = T)
system('rm -f tmp.wcinput.txt; rm -f tmp.wcinput.txt')
wcurse$z = abs(wcurse$discovery.beta / wcurse$discovery.se)
wcurse$p = pnorm(wcurse$z,lower.tail = F)/2
wcurse$maf = wcurse$discovery.freq
wcurse$maf[wcurse$discovery.freq>0.5] = 1-wcurse$discovery.freq[wcurse$discovery.freq>0.5]

# calculate power for winner's curse-corrected effects
message('Calculating power for winners curse-corrected effects.')
pow = c()
for (i in 1:nrow(wcurse)) {
  tmp = power_beta_maf(beta = wcurse$debiased.beta.mle[i], maf = wcurse$maf[i], n = replicationN, pval=replicationP)
  pow = c(pow,tmp)
}

# store expected number of successfull discoveries in logInfo
logInfo = paste0(logInfo,
  '\n--- calculations ---',
  '\ndiscoveries: ', nrow(df),
  '\nminimum replication power: ', min(pow),
  '\nmaximum replication power: ', max(pow),
  '\naverage replication power (mean): ', mean(pow),
  '\naverage replication power (median): ', median(pow),
  '\nexpected replications: ', sum(pow))

# create output
output = df[,c(cols2keep)]
output$debiased.beta.mle = wcurse$debiased.beta.mle
output$debiased.beta.mle.l95 = wcurse$l95.mle
output$debiased.beta.mle.u95 = wcurse$u95.mle
output$replication_power = pow
output[nrow(output)+1,] = NA
output[nrow(output),ncol(output)] = sum(output$replication_power, na.rm = T)

# writing 
message(sprintf('Writing %s.txt', outFile))
write.table(output, sprintf('%s.txt', outFile), sep = '\t', row.names = F, col.names = T, quote = F, na = c(""))

# save log file
message(sprintf('Writing %s.log', outFile))
sink(sprintf('%s.log', outFile)) 
sink(stdout(), type = "message")
message(logInfo)
sink()
system(sprintf('chmod -R 770 %s*', outFile))
message(paste0('--- Completed: Compute replication power ---\n'))
