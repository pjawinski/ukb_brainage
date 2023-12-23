#!/usr/bin/env Rscript

# ==================================================
# === Combine results of ANNOVAR enrichment test ===
# ==================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop(paste0('expected 3 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
traits=args[1] # traits="gap_gm,gap_wm,gap_gwm"
inFiles=args[2] # inFiles="results/gap_gm/fuma/FUMA_job174874/annov.stats.txt,results/gap_wm/fuma/FUMA_job174875/annov.stats.txt,results/gap_gwm/fuma/FUMA_job174876/annov.stats.txt"
outFile=args[3] # outFile="results/combined/annovar.txt"

message(paste0('\n--- Combine results of ANNOVAR enrichment test ---',
               '\ntraits: ', traits,
               '\ninFiles: ', inFiles,
               '\noutFile: ', outFile,'\n'))

# attach packages to current R session
for (pkg in c('dplyr', 'stringr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# transform variables
traits = str_split(traits, ',')[[1]]
inFiles = str_split(inFiles, ',')[[1]]

# open and merge datasets
message(' - loading and combinding results.')
for (i in 1:length(traits)) {
  
  # open datasets and keep relevant variables
  tmp = read.delim(inFiles[i], sep = '\t', header = T)
  tmp$fisher.P[tmp$fisher.P==0] = 1e-307
  tmp$fdr = p.adjust(tmp$fisher.P,method = 'BH')

  if (i ==1) {
    tmp = tmp[,c('annot', 'ref.count', 'ref.prop', 'count', 'prop', 'enrichment', 'fisher.P', 'fdr')]
  } else { 
    tmp = tmp[,c('annot', 'count', 'prop', 'enrichment', 'fisher.P', 'fdr')]
  }
  

  # create unique variable names
  names(tmp)[which(names(tmp) %in% c('count', 'prop', 'enrichment', 'fisher.P', 'fdr'))] = 
    paste0(traits[i],c('_count', '_prop', '_enrichment', '_fisher.P', '_fdr'))
  
  # merge datasets
  if (i ==1) { df = tmp } else { df = left_join(df,tmp, by = 'annot') }
}
  
# get top p-value and top |rho|
df$top_pvalue =  df[,grep('fisher.P',names(df))] %>% apply(1, FUN = min)
df$top_fdr =  df[,grep('fdr',names(df))] %>% apply(1, FUN = min)

# create output
message(' - creating output data frame.')
df$`NA` = ""
cols = c('annot', 'ref.count', 'ref.prop', 'top_pvalue','top_fdr')
for (i in 1:length(traits)) {
  cols = c(cols,'NA',paste0(traits[i],c('_count', '_prop', '_enrichment', '_fisher.P','_fdr')))
}
output = df[,cols]
output = output[order(output$top_fdr,output$top_pvalue),]

# write results
message(sprintf(' - writing %s',outFile))
write.table(output, file = outFile, row.names = F, quote = F, sep = '\t', na = "NA")
system(paste0('chmod 770 ', outFile))
message('-- Create table of ANNOVAR enrichment test results ---')
