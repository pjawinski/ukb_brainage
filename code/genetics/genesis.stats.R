#!/usr/bin/env Rscript

# ======================================================================
# === get stats results of genetic effect size distribution analysis ===
# ======================================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop(paste0('expected 4 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}                                                                                                                                                    

# set arguments
traits = args[1] # traits="gap_gm,gap_wm,gap_gwm,height,neur"
traitLabels = args[2] # traitLabels="Grey_matter,White_matter,Grey_and_white,Height,Neuroticism"
genesisFit = args[3] # genesisFit="results/gap_gm/genesis/genesis.fit3.Rds,results/gap_wm/genesis/genesis.fit3.Rds,results/gap_gwm/genesis/genesis.fit3.Rds,data/sumstats/genesis/mr_height_wood_2014/genesis.fit3.Rds,data/sumstats/genesis/04_neur_baselmans_2019/genesis.fit3.Rds"
outFile = args[4] # outFile="results/combined/genesis.stats"

logInfo = paste0('\n--- get stats of genetic effect size distribution analysis ---',
                 '\ntraits: ', traits,
                 '\ntraitLabels: ', traitLabels,
                 '\ngenesisFit: ', genesisFit,
                 '\noutFile: ', outFile,'\n')
message(logInfo)

# attach required packages
for (pkg in c('dplyr','GENESIS','MASS','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# transform variables
traits = str_split(traits, ',')[[1]]
traitLabels = str_split(traitLabels, ',')[[1]]
traitLabels = str_replace_all(traitLabels, "_", " ")
genesisFit = str_split(genesisFit, ',')[[1]]

# create outtable
out = matrix(NA,nrow = length(traits),ncol = 17)
colnames(out) = c('trait','nsnps','model','h2','h2_se','sSnps','sSnps_se','reqsample','reqsnps','h2_large','h2_large_se','h2_small','h2_small_se','sSnps_large','sSnps_large_se','logLikelihood','bic')

k = 0
for (i in 1:length(traits)) {

  # load data
  message(sprintf('[%d/%d] Loading %s',i,length(traits),genesisFit[i]))
  k = k+1
  tmp = readRDS(genesisFit[i])
  if (is.null(tmp$estimates$`Total time in fitting the 2-component model (in hours)`)) {
    modeltype = '3-component'
  } else {
    modeltype = '2-component'
  }

  # get estimates
  nsnps = tmp$estimates$`Total number of SNPs in the GWAS study after quality control`
  h2 = tmp$estimates$`Total heritability in log-odds-ratio scale (sd)`
  h2_se = h2 %>% gsub(pattern = '.*[(]', replacement = '') %>% gsub(pattern = '[)].*', replacement = '')
  h2 = h2 %>% gsub(pattern = '[ ].*', replacement = '')
  sSnps = tmp$estimates$`Number of sSNPs (sd)`
  sSnps_se = sSnps %>% gsub(pattern = '.*[(]', replacement = '') %>% gsub(pattern = '[)].*', replacement = '')
  sSnps = sSnps %>% gsub(pattern = '[ ].*', replacement = '')
  logLikelihood = tmp$estimates$`Composite log-likelihood of fitted model`
  bic = tmp$estimates$`Model selection related`$BIC

  # specific for 3-component model
  if (modeltype=='3-component') {
    h2_large = tmp$estimates$`Heritability explained by the cluster with larger variance component (sd)`
    h2_large_se = h2_large %>% gsub(pattern = '.*[(]', replacement = '') %>% gsub(pattern = '[)].*', replacement = '')
    h2_large = h2_large %>% gsub(pattern = '[ ].*', replacement = '')
    h2_small = tmp$estimates$`Heritability explained by the cluster with samller variance component`
    h2_small_se = h2_small %>% gsub(pattern = '.*[(]', replacement = '') %>% gsub(pattern = '[)].*', replacement = '')
    h2_small = h2_small %>% gsub(pattern = '[ ].*', replacement = '')
    sSnps_large = tmp$estimates$`Number of sSNPs in the cluster with larger variance component (sd)`
    sSnps_large_se = sSnps_large %>% gsub(pattern = '.*[(]', replacement = '') %>% gsub(pattern = '[)].*', replacement = '')
    sSnps_large = sSnps_large %>% gsub(pattern = '[ ].*', replacement = '')
  } else {
    h2_large = h2_large_se = h2_small = h2_small_se = sSnps_large = sSnps_large_se = NA
  }

  # project number of discoveries and explained variance
  GVpercentage = Numdiscoveries = matrix(NA,nrow=1001,ncol=2)
  l = 0
  for (j in seq(0, 10000000, length.out = 1001)) {
    l = l+1
    GVpercentage[l,1] = Numdiscoveries[l,1] = j
    if (modeltype == '3-component') {
    est = tmp$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates`
    } else {
    est = tmp$estimates$`Parameter (pic, sigmasq, a) estimates`
    }
    v = tmp$estimates$`Covariance matrix of parameter estimates`
    project = projection(est,v,n=j,CI=FALSE);
      GVpercentage[l,2] = project$GVpercentage[1]
      Numdiscoveries[l,2] = project$Numdiscoveries[1]
  }
  colnames(GVpercentage) = c('n','percentage')
  colnames(Numdiscoveries) = c('n','nsnps')
  GVpercentage = data.frame(GVpercentage)
  Numdiscoveries = data.frame(Numdiscoveries)
  reqsample = GVpercentage$n[(abs(GVpercentage$percentage - 80))==min(abs(GVpercentage$percentage - 80))]
  reqsnps = Numdiscoveries$nsnps[Numdiscoveries$n == reqsample]

  # add variables to out table
  out[i,] = c(traitLabels[i],nsnps,modeltype,h2,h2_se,sSnps,sSnps_se,reqsample,reqsnps,h2_large,h2_large_se,h2_small,h2_small_se,sSnps_large,sSnps_large_se,logLikelihood,bic)

  # remove variables
  rm(nsnps,h2,h2_se,h2_large,h2_large_se,h2_small,h2_small_se,sSnps,sSnps_se,reqsample,reqsnps,GVpercentage,Numdiscoveries,sSnps_large,sSnps_large_se,logLikelihood,bic,modeltype)
}
out = data.frame(out)

# write results
message(sprintf(' - writing %s.txt',outFile))
write.table(out, file = sprintf('%s.txt',outFile), row.names = F, quote = F, sep = '\t', na = "-")

# save log file
message(sprintf(' - writing %s.log', outFile))
sink(sprintf('%s.log', outFile))
sink(stdout(), type = "message")
message(logInfo)
sink()
system(sprintf('chmod -R 770 %s*', outFile))
message('--- Completed: get stats of genetic effect size distribution analysis ---')

