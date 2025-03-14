#!/usr/bin/env Rscript

# ===================================
# === Quality-check GWAMA results ===
# ===================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop(paste0('expected 4 arguments, but ', length(args), ' arguments provided.'), call.=FALSE)
}

# get arguments from command line
sumstats = args[1] # sumstats="results/gap_gm/replicate/metal.ivweight.gz"
outFile = args[2] # outFile="results/gap_gm/replicate/metal.ivweight.qc"
nCriterion = args[3] # nCriterion=TRUE | LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5
hetP = as.numeric(args[4]) # hetP=1E-6 | heterogeneity p threshold for excluding SNPs

logInfo = paste0('\n--- Quality-check GWAMA results ---',
               '\nsumstats: ', sumstats,
               '\noutFile: ', outFile,
               '\nnCriterion: ', nCriterion,
               '\nhetP: ', hetP,'\n')
message(logInfo)

# attach required packages
for (pkg in c('data.table')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# import data
message(paste0(' - importing data.'))
df = data.frame(fread(cmd=paste0("gzip -dc ", sumstats), tmpdir = getwd(), header=T, stringsAsFactors=FALSE))

# exclude snps that do not meet N criterion
if (nCriterion == TRUE) {
  message(' - excluding variants that do not meet N criterion (except for Y and MT).')
  maxN = max(df$N)
  minN = quantile(df$N, probs = 0.9)[[1]] / 1.5
  nExcl = sum(df$N < minN & df$CHR != "Y" & df$CHR != "MT")
  nIncl = sum(df$N >= minN | df$CHR == "Y" | df$CHR == "MT")
  df = df[df$N >= minN | df$CHR == "Y" | df$CHR == "MT",]
  logInfo = paste0(logInfo,
    '\n--- calculations ---',
    '\nmaxN: ', maxN,
    '\nminN: ', minN,
    '\nnExcl: ', nExcl,
    '\nnIncl: ', nIncl)

  message(' - excluding variants that do not meet N criterion (Y chromosome).')
  if (sum(df$CHR=="Y") > 0) {
    maxNY = max(df$N[df$CHR=="Y"])
    minNY = quantile(df$N[df$CHR=="Y"], probs = 0.9)[[1]] / 1.5
    nExclY = sum(df$N[df$CHR=="Y"] < minNY)
    nInclY = sum(df$N[df$CHR=="Y"] >= minNY)
    df = df[(df$CHR == "Y" & df$N >= minNY) | df$CHR != "Y",]
    logInfo = paste0(logInfo,
      '\nmaxNY: ', maxNY,
      '\nminNY: ', minNY,
      '\nnExclY: ', nExclY,
      '\nnInclY: ', nInclY)
  }

  message(' - excluding variants that do not meet N criterion (MT chromosome).')
  if (sum(df$CHR=="MT") > 0) {
    maxNMT = max(df$N[df$CHR=="MT"])
    minNMT = quantile(df$N[df$CHR=="MT"], probs = 0.9)[[1]] / 1.5
    nExclMT = sum(df$N[df$CHR=="MT"] < minNMT)
    nInclMT = sum(df$N[df$CHR=="MT"] >= minNMT)
    df = df[(df$CHR == "MT" & df$N >= minNMT) | df$CHR != "MT",]
    logInfo = paste0(logInfo,
      '\nmaxNY: ', maxNMT,
      '\nminNY: ', minNMT,
      '\nnExclY: ', nExclMT,
      '\nnInclY: ', nInclMT)
  }
}

# exclude snps that do not meet N criterion
message(' - excluding variants that do not meet hetP criterion.')
hetExcl = sum(df$HET_P <= hetP)
hetIncl = sum(df$HET_P > hetP)
df = df[df$HET_P > hetP,]
logInfo = paste0(logInfo,
  '\nhetExcl: ', hetExcl,
  '\nhetIncl: ', hetIncl)

# write output
message(' - writing file.')
data.table::fwrite(df, file = sprintf('%s.gz', outFile), sep = ' ', compress = 'gzip')
system(sprintf('chmod -R 770 %s.gz', outFile))

# save log file
sink(sprintf('%s.log', outFile))
sink(stdout(), type = "message")
message(logInfo)
sink()
system(sprintf('chmod -R 770 %s.log', outFile))
message('--- Completed: Quality-check GWAMA results ---')

