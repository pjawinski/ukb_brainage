#!/usr/bin/env Rscript

# ========================
# === Prioritize genes ===
# ========================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop(paste0('expected 3 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
inputFile=args[1] # inputFile="results/combined/snplevel.txt"
cols=args[2] # cols="FINEMAP_genes nonsyn_customPthresh SMR_eqtl SMR_sqtl GTEx_singleTissue GTEx_multiTissue PoPS"
outputFile=args[3] # outputFile="results/combined/snplevel.prio.txt"
message(paste0('\n--- Prioritize genes | Settings ---',
               '\ninputFile: ', inputFile,
               '\ncols: ', cols,
               '\noutputFile: ', outputFile, '\n'))

# attach packages to current R session
for (pkg in c('dplyr','stringr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# transform variables
cols = str_split(cols, ' ')[[1]]

# load data and add ensembl gene variable without gene version
message(sprintf(' - Loading %s',inputFile))
df = read.delim(inputFile, sep = '\t', header = T, quote = "")

# get credible set locus summary
message(' - scoring genes')
df$prioritized = NA
for (i in 1:nrow(df)) {
  k = 0
  for (j in 1:length(cols)) {
    if(!is.na(df[i,cols[j]])) {
      k = k + 1

      # extract genes and score them
      tmp.genes = df[i,cols[j]]
      tmp.genes = str_split(tmp.genes, ' \\| ')[[1]]
      if (cols[j] == 'FINEMAP_genes') {
        tmp.scores = sub('[)].*','', tmp.genes)
        tmp.scores = sub('.*[(]','', tmp.scores)
        tmp.genes = sub(' .*','', tmp.genes)
        tmp.df = data.frame(genes = tmp.genes, scores = as.numeric(tmp.scores))
      } else {
        tmp.genes = sub(' .*','', tmp.genes)
        n = length(tmp.genes)
        slices = (n*(n+1))/2
        tmp.df = data.frame(genes = tmp.genes, scores = (n:1)/slices)
      }

      # add scores across culomns
      if (k==1) { prio = tmp.df } else {
        prio = full_join(prio,tmp.df,'genes')
        prio$scores = rowSums(prio[,c('scores.x','scores.y')], na.rm = T)
        prio = prio[,c('genes','scores')]
      }
    }
  }
  if(exists('prio')) {
    prio = prio[order(-prio$scores),]
    df$prioritized[i] = prio$genes[1]
    #df$prioritized[i] = paste0(sprintf('%s (%0.2f)', prio$genes,prio$scores), collapse = ' | ')
    rm(prio)
  }
}

# write results
message(sprintf(' - writing %s',outputFile))
write.table(df, file = outputFile, row.names = F, quote = F, sep = '\t', na = "NA")
system(sprintf('chmod 750 %s',outputFile))
message('--- Completed: Prioritize genes')
