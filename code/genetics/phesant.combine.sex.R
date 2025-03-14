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
phesantSummary=args[2] # phesantSummary="results/combined/discovery.sex.male.phewas.txt,results/combined/discovery.sex.female.phewas.txt"
outputFile=args[3] # outputFile="results/combined/discovery.sex.phewas.txt"

message(paste0('\n--- Create sex-stratified suppl. table of PheWAS results ---',
               '\ntraits: ', traits,
               '\nphesantSummary: ', phesantSummary,
               '\noutputFile: ', outputFile,'\n'))

# attach packages to current R session
for (pkg in c('data.table', 'dplyr', 'stringr', 'tidyr')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# transform variables
phesantSummary = str_split(phesantSummary, ',')[[1]]
traits = str_split(traits, ',')[[1]]

# open and merge datasets
for (i in 1:2) {
  
  # open datasets and keep relevant variables
  tmp = read.delim(phesantSummary[i], sep = '\t', header = T)
  if (i ==1) {
    tmp = tmp[,c('varName','description','custom_category', 'Path', 'varType', 'resType', 'n', 'ntotal', do.call(paste0,crossing(traits,c('_beta', '_se', '_rho', '_pvalue', '_fdr'))))]
    startcol=5
  } else { 
    tmp = tmp[,c('varName', 'varType', 'resType',  'n', 'ntotal', do.call(paste0,crossing(traits,c('_beta', '_se', '_rho', '_pvalue', '_fdr'))))]
    startcol=2
  }
   
  # create unique variable names
  names(tmp)[startcol:ncol(tmp)] = 
    paste(letters[i],names(tmp)[startcol:ncol(tmp)], sep='_')
  
  # merge datasets
  if (i==1) { df = tmp } else { 
    df = inner_join(df,tmp, by = 'varName')
    df = df[df[paste0(letters[i],'_varType')]==df[paste0(letters[i-1],'_varType')],]
    df = df[df[paste0(letters[i],'_resType')]==df[paste0(letters[i-1],'_resType')],]
  }
}

# compare beta coefficients between groups (see Clogg et al. 1995, p. 1276, doi: 10.1086/230638)
for (i in 1:length(traits)) {
  abeta = df[paste('a',traits[i],'beta', sep = '_')]
  ase = df[paste('a',traits[i],'se', sep = '_')]
  bbeta = df[paste('b',traits[i],'beta', sep = '_')]
  bse = df[paste('b',traits[i],'se', sep = '_')]
  z = (abeta-bbeta)/sqrt(ase^2+bse^2)
  p = pnorm(abs(unlist(z)), lower.tail=FALSE)*2
  df[paste(traits[i],c('deltaZ'), sep = '_')] = z
  df[paste(traits[i],c('deltaP'), sep = '_')] = p
}

# get top p-value and top |rho|
df$top_pvalue =  df[,grep('pvalue',names(df))] %>% apply(1, FUN = min)
df$top_fdr = df[,grep('fdr',names(df))] %>% apply(1, FUN = min)
df$top_rhoAbs = df[,grep('rho',names(df))] %>% abs() %>% apply(1, FUN = max)
df$top_deltaP = df[,grep('deltaP',names(df))] %>% apply(1, FUN = min)

# create output
message('Writing txt file.')
df$`NA` = ""
cols = c('varName', 'description', 'top_deltaP', 'top_rhoAbs','top_pvalue','top_fdr')
for (i in 1:length(traits)) {
  for (j in 1:2) {
    cols = c(cols,'NA',paste(letters[j],traits[i],c('beta', 'se', 'rho', 'pvalue', 'fdr'), sep = '_'))
  }
  cols = c(cols,'NA',paste0(traits[i],c('_deltaP')))
}
cols = c(cols, c('NA', 'a_n', 'b_n', 'a_varType', 'a_resType', 'custom_category', 'Path'))
output = df[,cols]
output = output[order(output$top_fdr,output$top_pvalue),]
write.table(output, paste0(outputFile,'.txt'), sep = '\t', quote = F, row.names = F)
system(paste0('chmod 770 ', outputFile, '*'))
message('-- Completed: Create sex-stratified suppl. table of PheWAS results ---')
