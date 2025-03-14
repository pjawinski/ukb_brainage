#!/usr/bin/env Rscript

# ========================
# === get credible set ===
# ========================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 11) {
  stop(paste0('expected 10 arguments, but ', length(args), ' arguments provided.'), call.=FALSE)
}

# get arguments from command line
finemapCred = args[1] # finemapCred='1_215139887/finemap.cred2'
finemapSNP = args[2] # finemapSNP='1_215139887/finemap.snp'
indexLD = args[3] # indexLD='1_215139887/flankingSNPs.ld'
snpTotal = as.numeric(args[4]) # snpTotal=618
snpCausal = as.numeric(args[5]) # snpCausal=2.24
bestK = as.numeric(args[6]) # bestK=2
regionalh2 = args[7] # regionalh2="0.00204 [0.00136,0.00284]"
minBP = as.numeric(args[8]) # minBP=215134041
maxBP = as.numeric(args[9]) # maxBP=215289079
outputDF = args[10] # outputDF = '1_215139887/credibleDF.txt'
outputLS = args[11] # outputLS = '1_215139887/credibleLS.txt'

message(paste0('\n--- Creating credible set output files ---',
               '\nfinemapCred: ', finemapCred,
               '\nfinemapSNP: ', finemapSNP,
               '\nindexLD: ', indexLD,
               '\nsnpTotal: ', snpTotal,
               '\nsnpCausal: ', snpCausal,
               '\nbestK: ', bestK,
               '\nregionalh2: ', regionalh2,
               '\nminBP: ', minBP,
               '\nmaxBP: ', maxBP,
               '\noutputDF: ', outputDF,
               '\noutputLS: ', outputLS,'\n'))

# load packages
for (pkg in c('dplyr')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# load data
cred = read.table(finemapCred, header = T)
snp = read.table(finemapSNP, header = T)
ld = read.table(indexLD, header = T)

# create empty data frame with as many rows as SNPs in credible sets
cols = grep('cred',names(cred))
rownum = 0; for (i in cols) { rownum = rownum + sum(!is.na(cred[i])) }
df = data.frame(matrix(NA, nrow = rownum, ncol = 4))
names(df) = c('signalCount','rsid','prob','cumprob')

signal = 0
k = 0
for (i in cols) {
  signal = signal + 1
  for (j in 1:sum(!is.na(cred[i]))) {
    k = k + 1
    df$signalCount[k] = signal
    df$rsid[k] = cred[j,i]
    df$prob[k] = cred[j,i+1]
  }
  df$cumprob[df$signalCount==signal & !is.na(df$signalCount)] = cumsum(df$prob[df$signalCount==signal & !is.na(df$signalCount)])
}

# add R2 with index snp
ld = ld[,-which(names(ld) %in% c('CHR_B'))]
names(ld) = c('chromosome', 'index_bp', 'index_rsid', 'bp', 'rsid', 'index_r2')
df = left_join(df, ld, by = 'rsid')

# add size of credible set
df$cs.size = 0
setSize = as.data.frame(table(df$signalCount))
for (i in 1:nrow(df)) { 
  df$cs.size[i] = setSize$Freq[setSize$Var == df$signalCount[i]]
}

# add snpTotal, snpCausal, regionalh2, minBP, maxBP
df$snpTotal = snpTotal
df$snpCausal = snpCausal
df$bestK = bestK
df$regionalh2 = regionalh2
df$minBP = minBP
df$maxBP = maxBP
df$rangeBP = maxBP - minBP

# add .snp data
snp = snp[,-which(names(snp) %in% c('index', 'chromosome', 'position'))]
names(snp)[which(names(snp)=='prob')] = 'marginal_prob'
df = left_join(df, snp, by = 'rsid')

# get strings of concatenated credible set sizes, and credible snps
cssize = paste(df$cs.size[!duplicated(df$signalCount)],collapse = ' | ')
cssnps = df %>%
  group_by(signalCount) %>%
  summarise(y = paste0(rsid, collapse = ",")) %>%
  summarise(y = paste(y, collapse = " | ")) %>%
  pull(y)

# sort columns
DF = df[,c('index_rsid', 'chromosome', 'index_bp',  
           'regionalh2', 'minBP', 'maxBP', 'rangeBP', 'snpTotal', 'snpCausal', 'bestK',
           'signalCount', 'cs.size', 'rsid', 'bp', 'allele1', 'allele2', 'maf', 'beta', 'se', 'z', 'prob', 'cumprob',
           'marginal_prob', 'log10bf', 'mean', 'sd', 'mean_incl', 'sd_incl', 'index_r2')]
LS = cbind(df[1,c('index_rsid','chromosome','index_bp','regionalh2', 'minBP','maxBP','rangeBP','snpTotal','snpCausal','bestK')], data.frame(cssize = cssize, cssnps = cssnps))

# output
write.table(DF, outputDF, sep = "\t", row.names = F, col.names = T, quote = F)
write.table(LS, outputLS, sep = "\t", row.names = F, col.names = T, quote = F)
message(paste0('--- Creating credible set output files completed ---\n'))
