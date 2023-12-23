#!/usr/bin/env Rscript

# =========================================
# === determine significant SMR results ===
# =========================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop(paste0('expected 3 arguments, but ', length(args), ' arguments provided.'), call.=FALSE)
}

# get arguments from command line
smrInput = args[1] # smrInput='results/gap_gm/smr/smr.txt'
smrOutput = args[2] # smrOutput = 'results/gap_gm/smr/smr.filtered.txt'
multipleTesting = args[3] # multipleTesting = 'fdr'

message(paste0('\n--- SMR result filtering ---',
               '\nsmrInput: ', smrInput,
               '\nsmrOutput: ', smrOutput,'\n'))

# load packages
for (pkg in c('data.table', 'dplyr', 'ggplot2', 'ggrepel')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# import data
message(paste0('Importing data.'))
df = read.delim(smrInput, sep='\t', header=T, stringsAsFactors=FALSE)

# calculate FDR and Bonferroni
df$FDR = p.adjust(df$p_SMR, method = 'BH')
sigBonf = 0.05/nrow(df)

# keep significant variant-gene-trait associations
if (multipleTesting == 'fdr') {
  df = df[df$FDR < 0.05,]
} else if (multipleTesting == 'bonferroni') {
  df = df[df$p_SMR < sigBonf,]
}

# remove HEIDI outliers
df = df[df$p_HEIDI > 0.01 & !is.na(df$p_HEIDI),]

# save plot
message('Saving filtered smr results.')
write.table(df, file = smrOutput, sep = '\t', quote = F, row.names = F) 
system(paste0('chmod -R 770 ', smrOutput))
message(paste0('--- Filtering SMR results finished ---\n'))
