#!/usr/bin/env Rscript

# ==========================================================
# === create phesant summary text file and combined plot ===
# ==========================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop(paste0('expected 3 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
traits=args[1] # traits="gap_gm,gap_wm,gap_gwm"
surfcorrFiles=args[2] # surfcorrFiles="results/combined/discovery.sex.male.surfcorr.txt,results/combined/discovery.sex.female.surfcorr.txt"
outFile=args[3] # outputFile="results/combined/discovery.sex.surfcorr.txt"

message(paste0('\n--- Create sex-stratified suppl. table of FreeSurfer association result ---',
               '\ntraits: ', traits,
               '\nsurfcorrFiles: ', surfcorrFiles,
               '\noutFile: ', outFile,'\n'))

# attach packages to current R session
for (pkg in c('data.table', 'cocor', 'dplyr', 'stringr', 'tidyr')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# transform variables
surfcorrFiles = str_split(surfcorrFiles, ',')[[1]]
traits = str_split(traits, ',')[[1]]

# open and merge datasets
for (i in 1:2) {
  
  # open datasets and keep relevant variables
  tmp = read.delim(surfcorrFiles[i], sep = '\t', header = T)
  if (i ==1) {
    tmp = tmp[,c('id','structure','category', 'n', do.call(paste0,crossing(traits,c('_rho', '_pval', '_fdr'))))]
    startcol=4
  } else { 
    tmp = tmp[,c('id', 'n', do.call(paste0,crossing(traits,c('_rho', '_pval', '_fdr'))))]
    startcol=2
  }
   
  # create unique variable names
  names(tmp)[startcol:ncol(tmp)] = 
    paste(letters[i],names(tmp)[startcol:ncol(tmp)], sep='_')
  
  # merge datasets
  if (i==1) { df = tmp } else { df = inner_join(df,tmp, by = 'id') }
}

# compare correlations between independent groups
for (i in 1:length(traits)) {
    df[paste(traits[i],c('deltaZ'), sep = '_')] = NA
    df[paste(traits[i],c('deltaP'), sep = '_')] = NA
  for (j in 1:nrow(df)) {
    corr = cocor.indep.groups(df[j,paste('a',traits[i],'rho', sep = '_')], df[j, paste('b',traits[i],'rho', sep = '_')], df$a_n[j], df$b_n[j])
    df[j,paste(traits[i],c('deltaZ'), sep = '_')] = corr@fisher1925$statistic
    df[j,paste(traits[i],c('deltaP'), sep = '_')] = corr@fisher1925$p.value
  }
}

# get top p-value and top |rho|
df$top_pvalue =  df[,grep('pval',names(df))] %>% apply(1, FUN = min)
df$top_fdr = df[,grep('fdr',names(df))] %>% apply(1, FUN = min)
df$top_rhoAbs = df[,grep('rho',names(df))] %>% abs() %>% apply(1, FUN = max)
df$top_deltaP = df[,grep('deltaP',names(df))] %>% apply(1, FUN = min)

# create output
message('Writing txt file.')
df$`NA` = ""
cols = c('id', 'structure', 'category','a_n','b_n','top_deltaP', 'top_rhoAbs','top_pvalue','top_fdr')
for (i in 1:length(traits)) {
  for (j in 1:2) {
    cols = c(cols,'NA',paste(letters[j],traits[i],c('rho', 'pval', 'fdr'), sep = '_'))
  }
  cols = c(cols,'NA',paste0(traits[i],c('_deltaP')))
}
output = df[,cols]
output = output[order(output$top_fdr,output$top_pvalue),]
write.table(output, sprintf('%s',outFile), sep = '\t', quote = F, row.names = F)
system(paste0('chmod 770 ', outFile, '*'))
message('-- Completed: Create sex-stratified suppl. table of FreeSurfer association results ---')
