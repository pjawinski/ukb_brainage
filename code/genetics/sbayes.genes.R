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
annotFile=args[1] # annotFile="results/gap_gm/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.snpRes.annot.txt"
summaryFile=args[2] # summaryFile="results/gap_gm/gwama/eur/sbayes/sbayesrc.lcs.leadsnps"
outputFile=args[3] # outputFile="results/gap_gm/gwama/eur/sbayes/sbayesrc.lcs.leadsnps.genes"
message(paste0('\n--- Get list of genes sorted by cumulative PP of variants from finemapping | Settings ---',
               '\nannotFile: ', annotFile,
               '\nsummaryFile: ', summaryFile,
               '\noutputFile: ', outputFile, '\n'))

# attach packages to current R session
for (pkg in c('dplyr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# load data and add ensembl gene variable without gene version
message(' - Loading files')
annot = read.delim(annotFile, sep = '\t', header = T, quote = "")
smry = read.delim(summaryFile, sep = '\t', header = T, quote = "")

# get credible set locus summary
message(' - creating list of genes (sorted by variant PP) for each locus')
locus.summary = data.frame(matrix(NA, nrow = length(unique(annot$index_rsid)), ncol = 2))
names(locus.summary) = c('LEAD_SNP', 'genes')
i = 0
for (leadsnp in unique(annot$LEAD_SNP)) {
  i = i + 1
  annot.tmp = annot[annot$LEAD_SNP == leadsnp,]
  annot.tmp = annot.tmp[order(annot.tmp$csCount,annot.tmp$cumPIP),]

  gene.tmp = data.frame(matrix(NA, nrow = length(unique(annot.tmp$NEAREST_GENE)), ncol = 3))
  names(gene.tmp) = c('gene', 'cumpip', 'concat')
  j = 0
  for (gene in unique(annot.tmp$NEAREST_GENE)) {
    j = j + 1
    gene.tmp[j,1] = gene
    gene.tmp[j,2] = sum(annot.tmp$PIP[annot.tmp$NEAREST_GENE == gene])
    gene.tmp[j,3] = sprintf('%s (%0.3f)', gene.tmp[j,1], gene.tmp[j,2])
  }
  gene.tmp = gene.tmp[order(gene.tmp$cumpip, decreasing = T),]
  locus.summary[i,1] = leadsnp
  locus.summary[i,2] = paste(gene.tmp$concat, collapse = ' | ')
}

# write results
message(sprintf(' - writing %s',outputFile))
output = left_join(smry,locus.summary,by='LEAD_SNP')
write.table(output, file = outputFile, row.names = F, quote = F, sep = '\t', na = "NA")
system(sprintf('chmod 750 %s',outputFile))
message('--- Completed: Get list of genes sorted by cumulative PP of variants from finemapping')
