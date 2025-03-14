#!/usr/bin/env Rscript
# ===========================================================================
# === create phesant summary text file and combined plot for supplementum ===
# ===========================================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop(paste0('expected 3 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
traits=args[1] # traits="gap_gm,gap_wm,gap_gwm"
surfcorrFiles=args[2] # surfcorrFiles="results/gap_gm/surfplot/surfcorr.txt,results/gap_wm/surfplot/surfcorr.txt,results/gap_gwm/surfplot/surfcorr.txt"
outFile=args[3] # outFile="results/combined/surfcorr.txt"

message(paste0('\n--- Combine surfcorr results ---',
               '\ntraits: ', traits,
               '\nsurfcorrFiles: ', surfcorrFiles,
               '\noutFile: ', outFile,'\n'))

# attach packages to current R session
for (pkg in c('data.table', 'dplyr', 'stringr')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# transform variables
traits = str_split(traits, ',')[[1]]
surfcorrFiles= str_split(surfcorrFiles, ',')[[1]]

# open and merge datasets
for (i in 1:length(traits)) {
  
  # open datasets and keep relevant variables
  # set p to 1E-307 if observed is lower
  tmp = read.delim(surfcorrFiles[i], sep = '\t', header = T)
  tmp$pval[tmp$pval < 1E-307] = 1E-307
  tmp$fdr = p.adjust(tmp$pval, method = 'BH')
  tmp$rhoabs = abs(tmp$rho)
  if (i == 1) {
    tmp = tmp[,c('id', 'n', 'rho', 'rhoabs', 'pval', 'fdr')]
  } else { 
    tmp = tmp[,c('id', 'rho', 'rhoabs', 'pval', 'fdr')]
  }

  # create unique variable names
  names(tmp)[which(names(tmp) %in% c('rho', 'rhoabs','pval','fdr'))] = 
    paste0(traits[i],c('_rho', '_rhoabs', '_pval', '_fdr'))
  
  # merge datasets
  if (i == 1) { df = tmp } else { df = left_join(df,tmp, by = 'id') }
}
  
# get top p-value and top |rho|
df$top_pval =  df[,grep('pval',names(df))] %>% apply(1, FUN = min)
df$top_fdr = df[,grep('fdr',names(df))] %>% apply(1, FUN = min)
df$top_rhoabs = df[,grep('rhoabs',names(df))] %>% apply(1, FUN = max)

# get type (cortical / subcortical)'
df$category = 'subcortical volume'
df$category[grep('SurfArea',df$id)] = 'cortical surface area'
df$category[grep('ThickAvg',df$id)] = 'cortical thickness'
df$category[grep('GrayVol',df$id)] = 'cortical volume'

# create structure column
df$structure = df$id
df$structure[str_starts(df$structure,'Right')] = paste(df$structure[str_starts(df$structure,'Right')],'(right)')
df$structure[str_starts(df$structure,'Left')] = paste(df$structure[str_starts(df$structure,'Left')],'(left)')
df$structure[str_starts(df$structure,'RH')] = paste(df$structure[str_starts(df$structure,'RH')],'(right)')
df$structure[str_starts(df$structure,'LH')] = paste(df$structure[str_starts(df$structure,'LH')],'(left)')
df$structure = df$structure %>% 
  str_replace_all('\\.', ' ') %>%
  str_replace_all('RH_', '') %>%
  str_replace_all('LH_', '') %>%
  str_replace_all('LH_', '') %>%
  str_replace_all('Right ', '') %>%
  str_replace_all('Left ', '') %>%
  str_replace_all('SurfArea_', '') %>%
  str_replace_all('GrayVol_', '') %>%
  str_replace_all('ThickAvg_', '') %>%
  tolower() %>%
  str_replace('^\\w{1}', toupper)

# create output
message(sprintf('Writing %s',outFile))
df$`NA` = ""
cols = c('id','structure','category', 'n','top_pval','top_fdr','top_rhoabs')
for (i in 1:length(traits)) {
  cols = c(cols,'NA',paste0(traits[i],c('_rho', '_pval', '_fdr')))
}
output = df[,cols]
output = output[order(output$top_fdr,output$top_pval,-output$top_rhoabs),]
write.table(output, outFile, sep = '\t', quote = F, row.names = F)

message('-- Completed: Combine surfcorr results---')
