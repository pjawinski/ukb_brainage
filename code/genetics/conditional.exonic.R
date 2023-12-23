#!/usr/bin/env Rscript

# ======================================================================
# === Summarize results of nonsynonymous exotic variations per locus ===
# ======================================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop(paste0('expected 3 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
nsynFile=args[1] # nsynFile="results/gap_gm/conditional/conditional.cleaned.tophits.annovar.nonsynonymous.txt"
pthresh=as.numeric(args[2]) # pthresh=1E-6
outputFile=args[3] # outputFile="results/gap_gm/conditional/exonic.summary.txt"
message(paste0('\n--- Summarizing results for exonic variants ---',
               '\nnsynFile: ', nsynFile,
               '\npthresh: ', pthresh,
               '\noutputFile: ', outputFile, '\n'))

# attach packages to current R session
for (pkg in c('dplyr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# load data and add ensembl gene variable without gene version
message('Loading nsynFile.')
df = read.delim(nsynFile, sep = '\t', header = T)

# get summary for nonsynonymous variations with p < pthresh and those reaching genome-wide significance
message('Creating nsyn locus summary.')

for (pval in c(pthresh,5E-8)) {
  
  locus.summary = data.frame(matrix(NA, nrow = length(unique(df$LEAD_SNP)), ncol = 3))
  names(locus.summary) = c('LOCUS_COUNT', 'id', 'nonsyn')
  
  i = 0
  for (leadsnp in unique(df$LEAD_SNP)) {
    i = i + 1
    locus.summary[i,1:2] = df[df$LEAD_SNP == leadsnp & !duplicated(df$LEAD_SNP), c('LOCUS_CNT', 'LEAD_SNP')]
    if (sum(df$LEAD_SNP == leadsnp & df$P < pval) == 0) { next }
    df.tmp = df[df$LEAD_SNP == leadsnp & df$P < pval,]
    df.tmp = df.tmp[order(df.tmp$P),]
    
    gene.tmp = data.frame(matrix(NA, nrow = length(unique(df.tmp$GENE)), ncol = 5))
    names(gene.tmp) = c('gene', 'nsnps', 'minp', 'rsid_top', 'concat')
    j = 0
    for (gene in unique(df.tmp$GENE)) {
      j = j + 1
      gene.tmp[j,1] = gene
      gene.tmp[j,2] = nrow(df.tmp[df.tmp$GENE == gene,])
      gene.tmp[j,3] = df.tmp$P[df.tmp$GENE == gene][1]
      gene.tmp[j,4] = df.tmp$ID[df.tmp$GENE == gene][1]
      gene.tmp[j,5] = paste0(gene.tmp[j,1], ' (', gene.tmp[j,2], ';', gene.tmp[j,4], ';', formatC(gene.tmp[j,3], format = "e", digits = 0), ')')
    }
    
    gene.tmp = gene.tmp[order(gene.tmp$minp, decreasing = F),]
    locus.summary[i,3] = paste(gene.tmp$concat, collapse = ' | ')
  }
  assign(paste0('locus.summary.',formatC(pval)),locus.summary)
}

# write results
message('Writing summary file.')
output = full_join(get(paste0('locus.summary.',formatC(pthresh))), `locus.summary.5e-08`, by = 'id')
output$LOCUS_COUNT = output$LOCUS_COUNT.x; output$LOCUS_COUNT[is.na(output$LOCUS_COUNT)] = output$LOCUS_COUNT.y[is.na(output$LOCUS_COUNT)]
output = output[,c('LOCUS_COUNT.x','id','nonsyn.x','nonsyn.y')]
names(output) = c('LOCUS_COUNT','id','nonsyn.customPthresh','nonsyn.gwasPthresh')
output = output[order(output$LOCUS_COUNT),]
write.table(output, file = outputFile, row.names = F, quote = F, sep = '\t', na = "NA")
