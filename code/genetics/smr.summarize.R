#!/usr/bin/env Rscript

# ===========================================
# === Summarize SMR eQTL and sQTL results ===
# ===========================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop(paste0('expected 3 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
eqtlFile=args[1] # eqtlFile="results/gap_gm/smr/smr.eqtl.filtered.assigned.txt"
sqtlFile=args[2] # sqtlFile="results/gap_gm/smr/smr.sqtl.filtered.assigned.txt"
outputFile=args[3] # outputFile="results/gap_gm/smr/smr.summary.txt"
message(paste0('\n--- Summarize GTEx single- and multi-tissue eQTL hits | Settings ---',
               '\neqtlFile: ', eqtlFile,
               '\nsqtlFile: ', sqtlFile,
               '\noutputFile: ', outputFile, '\n'))

# attach packages to current R session
for (pkg in c('dplyr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# load data and add ensembl gene variable without gene version
message('Loading eqtlFile.')
df = read.delim(eqtlFile, sep = '\t', header = T)

# get eqtl locus summary
message('Creating eqtl locus summary.')
locus.summary = data.frame(matrix(NA, nrow = length(unique(df$COND_LEAD_SNP)), ncol = 3))
names(locus.summary) = c('LOCUS_COUNT', 'id', 'SMR_eqtl')
i = 0
for (leadsnp in unique(df$COND_LEAD_SNP)) {
  i = i + 1
  df.tmp = df[df$COND_LEAD_SNP == leadsnp,]
  df.tmp = df.tmp[order(df.tmp$p_SMR),]
  
  gene.tmp = data.frame(matrix(NA, nrow = length(unique(df.tmp$gene)), ncol = 5))
  names(gene.tmp) = c('gene', 'nprobes', 'minp', 'rsid_top', 'concat')
  j = 0
  for (gene in unique(df.tmp$Gene)) {
    j = j + 1
    gene.tmp[j,1] = gene
    gene.tmp[j,2] = nrow(df.tmp[df.tmp$Gene == gene,])
    gene.tmp[j,3] = df.tmp$p_SMR[df.tmp$Gene == gene][1]
    gene.tmp[j,4] = df.tmp$topSNP[df.tmp$Gene == gene][1]
    gene.tmp[j,5] = paste0(gene.tmp[j,1], ' (', gene.tmp[j,4], ';', formatC(gene.tmp[j,3], format = "e", digits = 0), ')')
  }
  gene.tmp = gene.tmp[order(gene.tmp$minp, decreasing = F),]
  locus.summary[i,1:2] = df[df$COND_LEAD_SNP == leadsnp & !duplicated(df$COND_LEAD_SNP), c('COND_LOCUS_COUNT', 'COND_LEAD_SNP')]
  locus.summary[i,3] = paste(gene.tmp$concat, collapse = ' | ')
}
eqtl.summary = locus.summary

# load sqtlFile
message('Loading sqtlFile.')
df = read.delim(sqtlFile, sep = '\t', header = T)

# get sqtl locus summary
message('Creating sqtl locus summary.')
locus.summary = data.frame(matrix(NA, nrow = length(unique(df$COND_LEAD_SNP)), ncol = 3))
names(locus.summary) = c('LOCUS_COUNT', 'id', 'SMR_sqtl')
i = 0
for (leadsnp in unique(df$COND_LEAD_SNP)) {
  i = i + 1
  df.tmp = df[df$COND_LEAD_SNP == leadsnp,]
  df.tmp = df.tmp[order(df.tmp$p_SMR),]
  
  gene.tmp = data.frame(matrix(NA, nrow = length(unique(df.tmp$gene)), ncol = 5))
  names(gene.tmp) = c('gene', 'nprobes', 'minp', 'rsid_top', 'concat')
  j = 0
  for (gene in unique(df.tmp$Gene)) {
    j = j + 1
    gene.tmp[j,1] = gene
    gene.tmp[j,2] = nrow(df.tmp[df.tmp$Gene == gene,])
    gene.tmp[j,3] = df.tmp$p_SMR[df.tmp$Gene == gene][1]
    gene.tmp[j,4] = df.tmp$topSNP[df.tmp$Gene == gene][1]
    gene.tmp[j,5] = paste0(gene.tmp[j,1], ' (', gene.tmp[j,4], ';', formatC(gene.tmp[j,3], format = "e", digits = 0), ')')
  }
  gene.tmp = gene.tmp[order(gene.tmp$minp, decreasing = F),]
  locus.summary[i,1:2] = df[df$COND_LEAD_SNP == leadsnp & !duplicated(df$COND_LEAD_SNP), c('COND_LOCUS_COUNT', 'COND_LEAD_SNP')]
  locus.summary[i,3] = paste(gene.tmp$concat, collapse = ' | ')
}
sqtl.summary = locus.summary

# write results
message('Writing summary file.')
output = full_join(eqtl.summary, sqtl.summary, by = 'id')
output$LOCUS_COUNT = output$LOCUS_COUNT.x; output$LOCUS_COUNT[is.na(output$LOCUS_COUNT)] = output$LOCUS_COUNT.y[is.na(output$LOCUS_COUNT)]
output = output[,c('LOCUS_COUNT', 'id', 'SMR_eqtl', 'SMR_sqtl')]
output = output[order(output$LOCUS_COUNT),]
write.table(output, file = outputFile, row.names = F, quote = F, sep = '\t', na = "NA")
