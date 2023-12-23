#!/usr/bin/env Rscript

# =========================================================
# === Summarize GTEx single- and multi-tissue eqtl hits ===
# =========================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop(paste0('expected 3 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
singleTissueFile=args[1] # singleTissueFile="results/gap_gm/eqtl/eqtl.singleTissue.txt"
multiTissueFile=args[2] # multiTissueFile="results/gap_gm/eqtl/eqtl.multiTissue.txt"
outputFile=args[3] # outputFile="results/gap_gm/eqtl/eqtl.summary.txt"
message(paste0('\n--- Summarize GTEx single- and multi-tissue eQTL hits | Settings ---',
               '\nsingleTissueFile: ', singleTissueFile,
               '\nmultiTissueFile: ', multiTissueFile,
               '\noutputFile: ', outputFile, '\n'))

# attach packages to current R session
for (pkg in c('dplyr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# load data and add ensembl gene variable without gene version
message('Loading singleTissueFile.')
df = read.delim(singleTissueFile, sep = '\t', header = T)

# get singleTissue locus summary
message('Creating singleTissue locus summary.')
locus.summary = data.frame(matrix(NA, nrow = length(unique(df$LEAD_SNP)), ncol = 3))
names(locus.summary) = c('LOCUS_COUNT', 'LEAD_SNP', 'GTEx_singleTissue')
i = 0
for (leadsnp in unique(df$LEAD_SNP)) {
  i = i + 1
  df.tmp = df[df$LEAD_SNP == leadsnp,]
  df.tmp = df.tmp[order(df.tmp$gene_id),]
  df.tmp$gene = df.tmp$hgnc_symbol; df.tmp$gene[df.tmp$gene == ""] = df.tmp$ensembl_gene_id[df.tmp$gene == ""]; df.tmp = df.tmp[order(df.tmp$gene),]
  
  gene.tmp = data.frame(matrix(NA, nrow = length(unique(df.tmp$gene)), ncol = 4))
  names(gene.tmp) = c('gene', 'ntissue', 'minp', 'concat')
  j = 0
  for (gene in unique(df.tmp$gene)) {
    j = j + 1
    gene.tmp[j,1] = gene
    gene.tmp[j,2] = nrow(df.tmp[df.tmp$gene == gene,])
    gene.tmp[j,3] = min(df.tmp$pval_nominal[df.tmp$gene == gene])
    gene.tmp[j,4] = paste0(gene.tmp[j,1], ' (', gene.tmp[j,2], ';', formatC(gene.tmp[j,3], format = "e", digits = 0), ')')
  }
  gene.tmp = gene.tmp[order(gene.tmp$ntissue, gene.tmp$minp, decreasing = c(T,F)),]
  locus.summary[i,1:2] = df[df$LEAD_SNP == leadsnp & !duplicated(df$LEAD_SNP), c('LOCUS_COUNT', 'LEAD_SNP')]
  locus.summary[i,3] = paste(gene.tmp$concat, collapse = ' | ')
}
singleTissue.summary = locus.summary

# load multiTissueFile
message('Loading multiTissueFile.')
df = read.delim(multiTissueFile, sep = '\t', header = T)

# get multiTissue locus summary
message('Creating multiTissue locus summary.')
locus.summary = data.frame(matrix(NA, nrow = length(unique(df$LEAD_SNP)), ncol = 3))
names(locus.summary) = c('LOCUS_COUNT', 'LEAD_SNP', 'GTEx_multiTissue')
i = 0
for (leadsnp in unique(df$LEAD_SNP)) {
  i = i + 1
  df.tmp = df[df$LEAD_SNP == leadsnp,]
  df.tmp = df.tmp[order(df.tmp$gene_id),]
  df.tmp$gene = df.tmp$hgnc_symbol; df.tmp$gene[df.tmp$gene == ""] = df.tmp$ensembl_gene_id[df.tmp$gene == ""]; df.tmp = df.tmp[order(df.tmp$gene),]
  
  gene.tmp = data.frame(matrix(NA, nrow = length(df.tmp$gene), ncol = 5))
  names(gene.tmp) = c('gene', 'mval0.9', 'ntissue', 'PVALUE_RE2', 'concat')
  j = 0
  for (gene in df.tmp$gene) {
    j = j + 1
    gene.tmp[j,1] = gene
    gene.tmp[j,2] = sum(df.tmp[df.tmp$gene == gene, which(grepl('^mval_', names(df.tmp)))] >= 0.9, na.rm = T)
    gene.tmp[j,3] = sum(!is.na(df.tmp[df.tmp$gene == gene, which(grepl('^mval_', names(df.tmp)))]))
    gene.tmp[j,4] = df.tmp$PVALUE_RE2[df.tmp$gene == gene]
    gene.tmp[j,5] = paste0(gene.tmp[j,1], ' (', gene.tmp[j,2], '/', gene.tmp[j,3], ';', formatC(gene.tmp[j,4], format = "e", digits = 0), ')')
  }
  gene.tmp = gene.tmp[order(gene.tmp$mval0.9/gene.tmp$ntissue, gene.tmp$mval0.9, gene.tmp$PVALUE_RE2, decreasing = c(T,T,F)),]
  locus.summary[i,1:2] = df[df$LEAD_SNP == leadsnp & !duplicated(df$LEAD_SNP), c('LOCUS_COUNT', 'LEAD_SNP')]
  locus.summary[i,3] = paste(gene.tmp$concat, collapse = ' | ')
}
multiTissue.summary = locus.summary

# write results
message('Writing summary file.')
output = full_join(singleTissue.summary, multiTissue.summary, by = 'LEAD_SNP')
output$LOCUS_COUNT = output$LOCUS_COUNT.x; output$LOCUS_COUNT[is.na(output$LOCUS_COUNT)] = output$LOCUS_COUNT.y[is.na(output$LOCUS_COUNT)]
output = output[,c('LOCUS_COUNT', 'LEAD_SNP', 'GTEx_singleTissue', 'GTEx_multiTissue')]
output = output[order(output$LOCUS_COUNT),]
write.table(output, file = outputFile, row.names = F, quote = F, sep = '\t', na = "NA")
