#!/usr/bin/env Rscript

# ==================================
# === Calcualte prs associations ===
# ==================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop(paste0('expected 5 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
trait = args[1] # trait="gap_wm"
dataFile = args[2] # dataFile="data/traits/replicate.txt"
prsFile = args[3] # prsFile="data/genetics/prs/imp_mri_qc_EUR/gap_wm.all_score"
covs = args[4] # covs="sex,age,age2,ac1,ac2,ac3,TIV,array,PanC1,PanC2,PanC3,PanC4,PanC5,PanC6,PanC7,PanC8,PanC9,PanC10"
outFile = args[5] # outFile="results/gap_gm/replicate/EUR/prs.assoc.txt"

message(paste0('\n--- Calcualte prs associations ---',
               '\ntrait: ', trait,
               '\ndataFile: ', dataFile,
               '\nprsFile: ', prsFile,
               '\ncovs: ', covs,
               '\noutFile: ', outFile,'\n'))

# attach required packages
for (pkg in c('dplyr','ppcor','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# transform variables
covs = str_split(covs, ',')[[1]]

# load data
message('Loading data.')
prs = read.delim(prsFile, header = T, sep = ' ')
prs$IID = as.numeric(sub('_.*','', prs$IID))
data = read.delim(dataFile, header = T, sep = '\t')
data = data[,c('IID',trait,covs)]

# keep only covariates with more than one value
message('Only keeping covariates with more than one distinct value.')
idx = sapply(data[,covs],function(x) length(unique(x))) > 1
if (sum(!idx) > 1) {
  message(sprintf(' - removing the following variables: %s', paste0(covs[!idx],collapse = ', ')))
  covs = covs[idx] } else {
  message(' - all covariates kept.')
}

# merge datasets
df = inner_join(prs,data,by = 'IID')

# run prs analysis
message('Running prs analysis.')
for (i in ncol(prs):3) {
  options(warn = -1)
  tmp = data.frame(pgsMethod = names(prs)[i], pcor.test(x = df[,names(prs)[i]], y = df[,trait], z = df[,covs], method = 'pearson'))
  tmp$R2 = tmp$estimate^2
  tmp$df = tmp$n-2-tmp$gp
  tmp = tmp[,c('pgsMethod','estimate','R2','statistic','df','p.value','n','gp','Method')]
  if (i == ncol(prs)) { results = tmp } else { results = rbind(results,tmp)}
  options(warn = 0)
}

# save file
message(sprintf('Writing results to file %s.',outFile))
write.table(results, outFile, sep = '\t', row.names = F, quote = F)
message('--- Completed: Calcualte prs associations ---')
