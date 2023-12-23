#!/usr/bin/env Rscript

# ==========================
# === combine rg.results ===
# ==========================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop(paste0('expected 3 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
traits = args[1] # traits="pwr_cz_alpha,pwr_cz_beta,pwr_cz_delta,pwr_cz_theta,pwr_cz_broadband,pwr_occ_alpha,pwr_occ_alphapeakfreq"
inputFiles = args[2] # inputFiles="results/pwr_cz_alpha/rgBrain/rg.results,results/pwr_cz_beta/rgBrain/rg.results,results/pwr_cz_delta/rgBrain/rg.results,results/pwr_cz_theta/rgBrain/rg.results,results/pwr_cz_broadband/rgBrain/rg.results,results/pwr_occ_alpha/rgBrain/rg.results,results/pwr_occ_alphapeakfreq/rgBrain/rg.results"
outputFile = args[3] # outputFile="results/combined/rgBrain.txt"

message(paste0('\n--- Combine results from genetic correlation analysis ---',
               '\ntraits: ', traits,
               '\ninputFiles: ', inputFiles,
               '\noutputFile: ', outputFile,'\n'))

# attach required packages
for (pkg in c('dplyr','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# transform variables
traits = str_split(traits, ',')[[1]]
inputFiles = str_split(inputFiles, ',')[[1]]

# open and merge datasets
for (i in 1:length(traits)) {
  
  # open datasets and keep relevant variables
  message(paste0('Loading ', inputFiles[i]))
  tmp = read.delim(inputFiles[i], sep = '\t', header = T)
  if (i ==1) {
    tmp = tmp[,c('p2','h2_obs','h2_obs_se','h2_int','h2_int_se','gcov_int','gcov_int_se','rg','se','z','p')]
  } else { 
    tmp = tmp[,c('p2','rg','se','z','p')]
  }
  
  # calculate FDR and abs(rg)
  tmp$FDR = p.adjust(tmp$p, method = 'BH')
  tmp$rgAbs = abs(tmp$rg)

  # create unique variable names
  names(tmp)[which(names(tmp) %in% c('rg','se','z','p','FDR','rgAbs'))] = 
    paste0(traits[i],c('_rg','_se','_z','_p','_FDR','_rgAbs'))
  
  # merge datasets
  if (i ==1) { df = tmp } else { df = left_join(df,tmp, by = 'p2') }
}
  
# get top p-value and top |rho|
df$top_rgAbs =  df[,grep('_rgAbs',names(df))] %>% apply(1, FUN = max)
df$top_p = df[,grep('_p',names(df))] %>% apply(1, FUN = min)
df$top_FDR = df[,grep('_FDR',names(df))] %>% apply(1, FUN = min)
df = df[order(df$top_FDR,df$top_p,-df$top_rg),]

# create output
message(sprintf('Writing %s',outputFile))
df$`NA` = ""
cols = c('p2','top_rgAbs','top_p','top_FDR')
for (i in 1:length(traits)) {
  cols = c(cols,'NA',paste0(traits[i],c('_rg','_se','_z','_p','_FDR')))
}
cols = c(cols,'NA','h2_obs','h2_obs_se','h2_int','h2_int_se','gcov_int','gcov_int_se')
output = df[,cols]
write.table(output, paste0(outputFile), sep = '\t', quote = F, row.names = F)
message('--- Completed: Combine results from genetic correlation analysis ---')

