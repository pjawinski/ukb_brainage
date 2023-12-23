#!/usr/bin/env Rscript

# ==================================================
# === Add literature matches to snplevel results ===
# ==================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop(paste0('expected 3 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
ownDiscoveryFile=args[1] # ownDiscoveryFile="results/combined/snplevel.txt"
noveltyFile=args[2] # noveltyFile="results/combined/novelty.clumping/clumped"
outputFile=args[3] # outputFile="results/combined/novelty.txt"
message(paste0('\n--- previous vs. literature discoveries ---',
               '\nownDiscoveryFile: ', ownDiscoveryFile,
               '\nnoveltyFile: ', noveltyFile,
               '\noutputFile: ', outputFile, '\n'))

# attach packages to current R session
for (pkg in c('dplyr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# load data
message('Loading Files.')
df = read.delim(noveltyFile, sep = '\t', header = T)
discov = read.delim(ownDiscoveryFile, sep = '\t', header = T)

# add studies with overlapping discoveries
message('Adding studies with overlapping discoveries.')
discov$LITERATURE_SHORT = discov$LITERATURE_LONG = ""
i = 0
for (discovery in discov$ID) {
  i = i + 1
  df.tmp = df[df$LOCUS_COUNT == unique(df$LOCUS_COUNT[df$ID == discovery]),]
  df.tmp = df.tmp[!duplicated(df.tmp$STUDY) & df.tmp$STUDY != 'ownDiscoveries',]
  
  if (nrow(df.tmp) > 0) {
    df.tmp$concat_long = df.tmp$concat_short = ""
    for (j in 1:nrow(df.tmp)) {
      df.tmp$concat_long[j] = paste0(df.tmp$STUDY[j], ' (', df.tmp$ID[j], ')')
      df.tmp$concat_short[j] = substr(df.tmp$STUDY[j], 1,1)
    }
    discov$LITERATURE_LONG[i] = paste(df.tmp$concat_long, collapse = ', ')
    discov$LITERATURE_SHORT[i] = paste(df.tmp$concat_short, collapse = ',')
  } else {
    discov$LITERATURE_LONG[i] = 'Novel'
    discov$LITERATURE_SHORT[i] = '-'
  }
}

# write results
message('Writing output file.')
write.table(discov, file = outputFile, row.names = F, quote = F, sep = '\t', na = "NA")
