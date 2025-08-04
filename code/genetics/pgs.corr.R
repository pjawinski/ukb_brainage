#!/usr/bin/env Rscript

# ==================================
# === Calculate PGS associations ===
# ==================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=8) {
  stop(paste0('expected 8 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
pgsMethod = args[1] # pgsMethod="sbayesrc"
pgsFile = args[2] # pgsFile="data/genetics/prs/imp_mri_qc_EUR/sbayesrc/gap_gm/sbayesrc.score"
phenoFile = args[3] # phenoFile="data/traits/replicate.txt"
joinVar = args[4] # joinVar="IID"
pgsVar = args[5] # pgsVar="SCORE1_SUM"
phenoVar = args[6] # phenoVar="gap_gm"
covs = args[7] # covs="sex,age,age2,ac1,ac2,ac3,TIV,array,PanC1,PanC2,PanC3,PanC4,PanC5,PanC6,PanC7,PanC8,PanC9,PanC10"
outFile = args[8] # outFile="results/gap_gm/replicate/EUR/sbayesrc.pgs.assoc.txt"

message(paste0('\n--- Calcualte PGS associations ---',
               '\npgsMethod: ', pgsMethod,
               '\npgsFile: ', pgsFile,
               '\nphenoFile: ', phenoFile,
               '\njoinVar: ', joinVar,
               '\npgsVar: ', pgsVar,
               '\nphenoVar: ', phenoVar,
               '\ncovs: ', covs,
               '\noutFile: ', outFile,'\n'))

# attach required packages
for (pkg in c('dplyr','psych','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# transform variables
covs = str_split(covs, ',')[[1]]

# load data
message('Loading data.')
pgs = read.delim(pgsFile, header = T, sep = '\t')
pgs = pgs[,c(joinVar,pgsVar)]
data = read.delim(phenoFile, header = T, sep = '\t')
data = data[,c(joinVar,phenoVar,covs)]

# merge datasets
df = inner_join(pgs,data,by = joinVar)
df = df[complete.cases(df),]

# keep only covariates with more than one value
message('Only keeping covariates with more than one distinct value.')
idx = sapply(df[,covs],function(x) length(unique(x))) > 1
if (sum(!idx) > 0) {
  message(sprintf(' - removing the following variables: %s', paste0(covs[!idx],collapse = ', ')))
  covs = covs[idx] } else {
  message(' - all covariates kept.')
}

# run pgs analysis
message('Running pgs analysis.\n')
rho = partial.r(data = df, x = c(phenoVar,pgsVar), y = covs, use = 'pairwise', method = 'pearson')[[2]]
n = nrow(df)
p = corr.p(rho,n-length(covs),ci = F,adjust="none")$p
t = rho * sqrt(n-2)/sqrt(1-rho^2) 
results = data.frame(pgsMethod = pgsMethod, estimate = rho, R2 = rho^2, statistic = t, df = n-2-length(covs), p.value = p, n = n, gp = length(covs), Method = 'pearson')
print(results)

# save file
message(sprintf('\nWriting %s',outFile))
write.table(results, outFile, sep = '\t', row.names = F, quote = F)
system(sprintf('chmod 770 %s',outFile))
message('--- Completed: Calcualte PGS associations ---')
