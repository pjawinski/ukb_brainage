#!/usr/bin/env Rscript

# ======================================================================
# === Summarize results of nonsynonymous exotic variations per locus ===
# ======================================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=6) {
  stop(paste0('expected 6 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
nsynFile=args[1] # nsynFile="results/gap_gm/credibleSet/credibleSet.df.annot.nonsynonymous.txt"
locusId=args[2] # locusId="index_rsid"
id=args[3] # id="rsid"
beta_se=args[4] # beta_se="beta,se"
pthresh=as.numeric(args[5]) # pthresh=1
outputFile=args[6] # outputFile="results/gap_gm/credibleSet/credibleSet.df.annot.nonsynonymous.summary.txt"
message(paste0('\n--- Summarizing results for exonic variants ---',
               '\nsynFile: ', nsynFile,
               '\nlocusId: ', locusId,
               '\nid: ', id,
               '\nbeta_se: ', beta_se,
               '\npthresh: ', pthresh,
               '\noutputFile: ', outputFile, '\n'))

# attach packages to current R session
for (pkg in c('dplyr','stringr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# transform variables
beta_se = str_split(beta_se, ',')[[1]]
beta = beta_se[1]
se = beta_se[2]

# load data and add ensembl gene variable without gene version
message('Loading nsynFile.')
df = read.delim(nsynFile, sep = '\t', header = T)

# calculate p
df$P = pnorm(abs(df[[beta]]/df[[se]]), lower.tail = F)*2

# get summary for nonsynonymous variations with p < pthresh and those reaching genome-wide significance
message('Creating nsyn locus summary.')

for (pval in c(pthresh,5E-8)) {
  
  locus.summary = data.frame(matrix(NA, nrow = length(unique(df[[locusId]])), ncol = 2))
  names(locus.summary) = c('id', 'nonsyn')
  
  i = 0
  for (leadsnp in unique(df[[locusId]])) {
    i = i + 1
    locus.summary[i,'id'] = df[df[[locusId]] == leadsnp & !duplicated(df[[locusId]]), locusId]
    if (sum(df[[locusId]] == leadsnp & df$P < pval) == 0) { next }
    df.tmp = df[df[[locusId]] == leadsnp & df$P < pval,]
    df.tmp = df.tmp[order(df.tmp$P),]
    
    gene.tmp = data.frame(matrix(NA, nrow = length(unique(df.tmp$GENE)), ncol = 5))
    names(gene.tmp) = c('gene', 'nsnps', 'minCADD_PHRED', 'rsid_top', 'concat')
    j = 0
    for (gene in unique(df.tmp$GENE)) {
      j = j + 1
      gene.tmp[j,1] = gene
      gene.tmp[j,2] = nrow(df.tmp[df.tmp$GENE == gene,])
      gene.tmp[j,3] = df.tmp$CADD_PHRED[df.tmp$GENE == gene][1]
      gene.tmp[j,4] = df.tmp[df.tmp$GENE == gene,id][1]
      gene.tmp[j,5] = sprintf('%s (%d;%s;%0.3f)', gene.tmp[j,1], gene.tmp[j,2], gene.tmp[j,4], gene.tmp[j,3])
    }
    
    gene.tmp = gene.tmp[order(gene.tmp$minCADD_PHRED, decreasing = T),]
    locus.summary[i,2] = paste(gene.tmp$concat, collapse = ' | ')
  }
  assign(paste0('locus.summary.',formatC(pval)),locus.summary)
}

# write results
message(sprintf('Writing %s',outputFile))
output = full_join(get(paste0('locus.summary.',formatC(pthresh))), `locus.summary.5e-08`, by = 'id')
names(output) = c('id','nonsyn.customPthresh','nonsyn.gwasPthresh')
write.table(output, file = outputFile, row.names = F, quote = F, sep = '\t', na = "NA")
system(sprintf('chmod 750 %s',outputFile))
