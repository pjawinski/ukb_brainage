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
inFiles = args[2] # inFiles="results/gap_gm/gwama/eur/gsmr/gsmr.multi,results/gap_wm/gwama/eur/gsmr/gsmr.multi,results/gap_gwm/gwama/eur/gsmr/gsmr.multi"
outFile = args[3] # outFile="results/combined/gwama.eur.gsmr.multi.txt"
direction = args[4] # direction="both" | only export results of 'exposure','outcome', or 'both'
mInstruments = as.numeric(args[5]) # mInstruments=-1 # set criterion for minimum number of instruments

message(paste0('\n--- Combine results from Mendelian Randomization ---',
               '\ntraits: ', traits,
               '\ninFiles: ', inFiles,
               '\noutFile: ', outFile,
               '\ndirection: ', direction,
               '\nmInstruments: ', mInstruments,'\n'))

# attach required packages
for (pkg in c('dplyr','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# transform variables
traits = str_split(traits, ',')[[1]]
inFiles = str_split(inFiles, ',')[[1]]

# open and merge datasets
for (i in 1:length(traits)) {
  
  # open datasets and keep relevant variables
  message(sprintf('Loading %s', inFiles[i]))
  tmp = read.delim(inFiles[i], sep = '\t', header = T)
  tmp_exposure = tmp[tmp$exposure == traits[i],c('outcome','n_snps_total','n_snps_HEIDI','gsmr_beta','gsmr_se','gsmr_p','ivw_p','divw_p','pivw_p','median_p','egger_p','ml_p','mbe_p','cm_p','lasso_p','avg_p','sum_p05')]
  tmp_outcome = tmp[tmp$outcome == traits[i],c('exposure','n_snps_total','n_snps_HEIDI','gsmr_beta','gsmr_se','gsmr_p','ivw_p','divw_p','pivw_p','median_p','egger_p','ml_p','mbe_p','cm_p','lasso_p','avg_p','sum_p05')]
  tmp_exposure$gsmr_fdr = p.adjust(tmp_exposure$gsmr_p, method = 'BH')
  tmp_outcome$gsmr_fdr = p.adjust(tmp_outcome$gsmr_p, method = 'BH')
  names(tmp_exposure) = c('trait',paste0(paste0(traits[i],'_exposure_'),c('n_snps_total','n_snps_HEIDI','gsmr_beta','gsmr_se','gsmr_p','ivw_p','divw_p','pivw_p','median_p','egger_p','ml_p','mbe_p','cm_p','lasso_p','avg_p','sum_p05','gsmr_fdr')))
  names(tmp_outcome) = c('trait',paste0(paste0(traits[i],'_outcome_'),c('n_snps_total','n_snps_HEIDI','gsmr_beta','gsmr_se','gsmr_p','ivw_p','divw_p','pivw_p','median_p','egger_p','ml_p','mbe_p','cm_p','lasso_p','avg_p','sum_p05','gsmr_fdr')))
  tmp = full_join(tmp_exposure,tmp_outcome, by = 'trait')
  
  # merge datasets
  if (i ==1) { df = tmp } else { df = left_join(df,tmp, by = 'trait') }
}

# create output
message(' - creating output data frame.')
df$`NA` = ""
df$top_gsmr_p = NA
df$top_gsmr_fdr = NA
cols = c('trait','top_gsmr_p','top_gsmr_fdr')
for (i in 1:length(traits)) {
  cols = c(cols,'NA',paste0(paste0(traits[i],'_exposure_'),c('n_snps_total','n_snps_HEIDI','gsmr_beta','gsmr_se','gsmr_p','gsmr_fdr','ivw_p','divw_p','pivw_p','median_p','egger_p','ml_p','mbe_p','cm_p','lasso_p','sum_p05')))
  cols = c(cols,'NA',paste0(paste0(traits[i],'_outcome_'),c('n_snps_total','n_snps_HEIDI','gsmr_beta','gsmr_se','gsmr_p','gsmr_fdr','ivw_p','divw_p','pivw_p','median_p','egger_p','ml_p','mbe_p','cm_p','lasso_p','sum_p05')))
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
output = output[!(rowSums(is.na(output[,grep('_p$',names(output))])) == ncol(output[,grep('_p$',names(output))])),]

# get top p-value and top FDR
message(' - getting top p and top FDR.')
output$top_gsmr_p = output[,grep('_gsmr_p',names(output))] %>% apply(1, FUN = min, na.rm = T)
output$top_gsmr_fdr = output[,grep('_gsmr_fdr',names(output))] %>% apply(1, FUN = min, na.rm = T)
output = output[order(output$top_gsmr_fdr,output$top_gsmr_p),]

# instrument criterion
if (mInstruments > 0) {
  message(sprintf(' - applying instrument criterion (m < %d).',mInstruments))
  idx = !(output[,grep('_n_snps',names(output))] %>% apply(1, FUN = min, na.rm = T) < mInstruments)
  output = output[idx,]
}

# write file
message(sprintf('Writing %s',outFile))
write.table(output, outFile, sep = '\t', quote = F, row.names = F)
message('--- Completed: Combine results from Mendelian Randomization ---')


