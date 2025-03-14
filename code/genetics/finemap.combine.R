#!/usr/bin/env Rscript

# ========================================================
# === Create table of FINEMAP results for supplementum ===
# ========================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop(paste0('expected 5 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
snplevelFile=args[1] # snplevelFile="results/combined/snplevel.txt"
finemapFiles=args[2] # finemapFiles="results/gap_gm/credibleSet/credibleSet.df.txt,results/gap_wm/credibleSet/credibleSet.df.txt,results/gap_gwm/credibleSet/credibleSet.df.txt"
outputFile=args[3] # outputFile="results/combined/finemap.suppl.txt"
traits=args[4] # traits="gap_gm,gap_wm,gap_gwm"
traitNames=args[5] # traitNames="grey matter,white matter,grey and white matter"

message(paste0('\n--- Creating table of FINEMAP results for supplementum ---',
               '\nsnplevelFile: ', snplevelFile,
               '\nfinemapFiles: ', finemapFiles,
               '\noutputFile: ', outputFile,
               '\ntraits: ', traits,
               '\ntraitNames: ', traitNames, '\n'))

# attach packages to current R session
for (pkg in c('dplyr', 'stringr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# transform variables
traits = str_split(traits, ',')[[1]]
traitNames = str_split(traitNames, ',')[[1]]
finemapFiles = str_split(finemapFiles, ',')[[1]]

# load data
message('Loading Files.')
snplevel = read.delim(snplevelFile, sep = '\t', header = T)
snplevel$key = paste0(snplevel$TRAIT,"_", snplevel$ID)
snplevel = snplevel[,c('LOCUS_COUNT','DISCOV_COUNT','key')]

# add traitNames and identifier
for (i in 1:length(finemapFiles)) {
  tmp = read.table(finemapFiles[i], sep = '\t', head = T, quote = "")
  tmp$trait = traitNames[i]
  tmp$key = paste0(traits[i],"_", tmp$index_rsid)
  if (i == 1) { finemap = tmp } else { finemap = rbind(finemap, tmp)}
}

# combine
message('Assign finemap to snplevel results.')
df = left_join(finemap, snplevel, by = 'key')
df = df[order(df$LOCUS_COUNT, df$DISCOV_COUNT, df$signalCount, df$cumprob),]
df$credibleSet_size = 0
setSize = as.data.frame(table(paste0(df$trait,df$index_rsid)))
for (i in 1:nrow(df)) { 
  df$credibleSet_size[i] = setSize$Freq[setSize$Var == paste0(df$trait[i],df$index_rsid[i])]
}

# create supplementum data.frame
suppl = df[,c('LOCUS_COUNT','trait','index_rsid','chromosome',
              'index_bp','regionalh2','minBP','maxBP','rangeBP','snpTotal','snpCausal','bestK','signalCount','cs.size',
              'rsid','bp','allele1','allele2','maf','beta','se','z','prob','cumprob',
              'marginal_prob','log10bf','mean','sd','mean_incl','sd_incl','REGION','NEAREST_GENE','DISTANCE')]

# write results
message('Writing output file.')
write.table(suppl, file = outputFile, row.names = F, quote = F, sep = '\t', na = "NA")
system(paste0('chmod 770 ', outputFile))
message('-- Creating suppl. table of finemap results completed. ---')
