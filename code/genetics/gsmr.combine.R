#!/usr/bin/env Rscript

# ===============================================
# === combine Mendelian Randomization results ===
# ===============================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop(paste0('expected 5 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
traits = args[1] # traits="gap_gm,gap_wm,gap_gwm"
inputFiles = args[2] # inputFiles="results/gap_gm/gsmr/gsmr.gsmr,results/gap_wm/gsmr/gsmr.gsmr,results/gap_gwm/gsmr/gsmr.gsmr"
outputFile = args[3] # outputFile="results/combined/gsmr.txt"
direction = args[4] # direction="exposure" | only export results of 'exposure','outcome', or 'both'
mInstruments = as.numeric(args[5]) # mInstruments=30 # set criterion for minimum number of instruments

message(paste0('\n--- Combine results from Mendelian Randomization ---',
               '\ntraits: ', traits,
               '\ninputFiles: ', inputFiles,
               '\noutputFile: ', outputFile,
               '\ndirection: ', direction,
               '\nmInstruments: ', mInstruments,'\n'))

# attach required packages
for (pkg in c('dplyr','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# transform variables
traits = str_split(traits, ',')[[1]]
inputFiles = str_split(inputFiles, ',')[[1]]

# open and merge datasets
for (i in 1:length(traits)) {
  
  # open datasets and keep relevant variables
  message(sprintf('Loading %s', inputFiles[i]))
  tmp = read.delim(inputFiles[i], sep = '\t', header = T)
  tmp_exposure = tmp[tmp$Exposure == traits[i],c('Outcome','bxy','se','p','nsnp')]
  tmp_outcome = tmp[tmp$Outcome == traits[i],c('Exposure','bxy','se','p','nsnp')]
  tmp_exposure$FDR = p.adjust(tmp_exposure$p, method = 'BH')
  tmp_outcome$FDR = p.adjust(tmp_outcome$p, method = 'BH')
  names(tmp_exposure) = c('trait',paste0(traits[i],c('_exposure_bxy','_exposure_se','_exposure_p','_exposure_nsnp','_exposure_FDR')))
  names(tmp_outcome) = c('trait',paste0(traits[i],c('_outcome_bxy','_outcome_se','_outcome_p','_outcome_nsnp','_outcome_FDR')))
  tmp = full_join(tmp_exposure,tmp_outcome, by = 'trait')
  
  # merge datasets
  if (i ==1) { df = tmp } else { df = left_join(df,tmp, by = 'trait') }
}

# create output
message(' - creating output data frame.')
df$`NA` = ""
df$top_p = NA
df$top_FDR = NA
cols = c('trait','top_p','top_FDR')
for (i in 1:length(traits)) {
  cols = c(cols,'NA',paste0(traits[i],c('_exposure_bxy','_exposure_se','_exposure_p','_exposure_nsnp','_exposure_FDR')))
  cols = c(cols,'NA',paste0(traits[i],c('_outcome_bxy','_outcome_se','_outcome_p','_outcome_nsnp','_outcome_FDR')))
}
output = df[,cols]

# exclude columns
if (direction == 'exposure') {
  message(' - keeping gsmr effects where trait is exposure.')
  output = output[,-grep('outcome',names(output))]
} else if (direction == 'outcome') {
  message(' - keeping gsmr effects where trait is outcome.')
  output = output[,-grep('exposure',names(output))]
}

idx = c()
for (i in 2:ncol(output)) {
  if (str_starts(names(output)[[i]],'NA') & str_starts(names(output)[[i-1]],'NA')) { idx = c(idx,i) }
}
if (length(idx > 0)) { output = output[,-idx] }

# exclude rows containing only NAs
message(' - excluding rows containing only NAs.')
output = output[!(rowSums(is.na(output[,grep('_p',names(output))])) == ncol(output[,grep('_p',names(output))])),]

# get top p-value and top FDR
message(' - getting top p and top FDR.')
output$top_p = output[,grep('_p',names(output))] %>% apply(1, FUN = min, na.rm = T)
output$top_FDR = output[,grep('_FDR',names(output))] %>% apply(1, FUN = min, na.rm = T)
output = output[order(output$top_FDR,output$top_p),]

# instrument criterion
if (mInstruments > 0) {
  message(sprintf(' - applying instrument criterion (m < %d).',mInstruments))
  idx = !(output[,grep('_nsnp',names(output))] %>% apply(1, FUN = min, na.rm = T) < mInstruments)
  output = output[idx,]
}

# write file
message(sprintf('Writing %s',outputFile))
write.table(output, outputFile, sep = '\t', quote = F, row.names = F)
message('--- Completed: Combine results from Mendelian Randomization ---')


