#!/usr/bin/env Rscript

# ===================================================
# === Get relevant genes from finemapping results ===
# ===================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop(paste0('expected 3 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
dfFile=args[1] # dfFile="results/gap_gm/credibleSet/credibleSet.df.annot.txt"
lsFile=args[2] # lsFile="results/gap_gm/credibleSet/credibleSet.ls.txt"
outputFile=args[3] # outputFile="results/gap_gm/credibleSet/credibleSet.ls.genes.txt"
message(paste0('\n--- Get list of genes sorted by cumulative PP of variants from finemapping | Settings ---',
               '\ndfFile: ', dfFile,
               '\nlsFile: ', lsFile,
               '\noutputFile: ', outputFile, '\n'))

# attach packages to current R session
for (pkg in c('dplyr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# load data and add ensembl gene variable without gene version
message(' - Loading files')
df = read.delim(dfFile, sep = '\t', header = T, quote = "")
ls = read.delim(lsFile, sep = '\t', header = T, quote = "")

# get credible set locus summary
message(' - creating list of genes (sorted by variant PP) for each locus')
locus.summary = data.frame(matrix(NA, nrow = length(unique(df$index_rsid)), ncol = 2))
names(locus.summary) = c('index_rsid', 'genes')
i = 0
for (leadsnp in unique(df$index_rsid)) {
  i = i + 1
  df.tmp = df[df$index_rsid == leadsnp,]
  df.tmp = df.tmp[order(df.tmp$signalCount,df.tmp$cumprob),]

  gene.tmp = data.frame(matrix(NA, nrow = length(unique(df.tmp$NEAREST_GENE)), ncol = 3))
  names(gene.tmp) = c('gene', 'cumprob', 'concat')
  j = 0
  for (gene in unique(df.tmp$NEAREST_GENE)) {
    j = j + 1
    gene.tmp[j,1] = gene
    gene.tmp[j,2] = sum(df.tmp$prob[df.tmp$NEAREST_GENE == gene])
    gene.tmp[j,3] = sprintf('%s (%0.3f)', gene.tmp[j,1], gene.tmp[j,2])
  }
  gene.tmp = gene.tmp[order(gene.tmp$cumprob, decreasing = T),]
  locus.summary[i,1] = leadsnp
  locus.summary[i,2] = paste(gene.tmp$concat, collapse = ' | ')
}

# write results
message(sprintf(' - writing %s',outputFile))
output = left_join(ls,locus.summary,by='index_rsid')
write.table(output, file = outputFile, row.names = F, quote = F, sep = '\t', na = "NA")
system(sprintf('chmod 750 %s',outputFile))
message('--- Completed: Get list of genes sorted by cumulative PP of variants from finemapping')
