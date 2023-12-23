#!/usr/bin/env Rscript

# ========================================================
# === Combine results of GOfuncR gsea enrichment tests ===
# ========================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop(paste0('expected 3 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
traits=args[1] # traits="gap_gm,gap_wm,gap_gwm"
gseaFiles=args[2] # gseaFiles="results/gap_gm/gofuncr/gofuncr.gsea.results.txt,results/gap_wm/gofuncr/gofuncr.gsea.results.txt,results/gap_gwm/gofuncr/gofuncr.gsea.results.txt"
outFile=args[3] # outFile="results/combined/gofuncr.gsea.txt"

message(paste0('\n--- Combine results of GOfuncR gene set enrichment tests ---',
               '\ntraits: ', traits,
               '\ngseaFiles: ', gseaFiles,
               '\noutFile: ', outFile,'\n'))

# attach packages to current R session
for (pkg in c('dplyr', 'stringr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# transform variables
traits = str_split(traits, ',')[[1]]
gseaFiles = str_split(gseaFiles, ',')[[1]]

# open and merge datasets
message(' - loading and combining results.')
for (i in 1:length(traits)) {
  
  # open datasets and keep relevant variables
  tmp = data.table::fread(gseaFiles[i], sep = '\t', header = T)

  # keep relevant columns
  if (i == 1) {
    tmp = tmp[,c('node_id','node_name','ontology','size','locus_count','rank','meanRank','raw_p_low_rank','FWER_low_rank','refined_p_low_rank','refined_FWER_low_rank')]
  } else { 
    tmp = tmp[,c('node_id','rank','meanRank','raw_p_low_rank','FWER_low_rank','refined_p_low_rank','refined_FWER_low_rank')]
  }
  
  # create unique variable names
  names(tmp)[which(names(tmp) %in% c('rank','meanRank','raw_p_low_rank','FWER_low_rank','refined_p_low_rank','refined_FWER_low_rank'))] = 
    paste0(traits[i],c('_rank','_meanRank','_raw_p_low_rank','_FWER_low_rank','_refinedP_low_rank','_refinedFWER_low_rank'))
  
  # merge datasets
  if (i == 1) { df = tmp } else { df = left_join(df,tmp, by = 'node_id') }
}
  
# get top p-value and top |rho|
df = data.frame(df)
df$top_pvalue =  df[,grep('_raw_p_low_rank',names(df))] %>% apply(1, FUN = min, na.rm = T)
df$top_fwer =  df[,grep('_FWER_low_rank',names(df))] %>% apply(1, FUN = min, na.rm = T)
df$top_refinedP =  suppressWarnings(df[,grep('_refinedP_low_rank',names(df))] %>% apply(1, FUN = min, na.rm = T))
df$top_refinedFWER =  suppressWarnings(df[,grep('_refinedFWER_low_rank',names(df))] %>% apply(1, FUN = min, na.rm = T))
df$top_refinedP[df$top_refinedP == Inf] = NA
df$top_refinedFWER[df$top_refinedFWER == Inf] = NA
df = df[order(df$top_refinedFWER,df$top_fwer,df$top_pvalue),]

# convert p = 0 to p < 0.001
for (i in c(grep('_FWER_low_rank',names(df)),grep('_refinedFWER_low_rank',names(df)),grep('top_fwer',names(df)),grep('top_refinedFWER',names(df)))) {
  df[df[,i]==0 & !is.na(df[,i]),i] = '< 0.001'
}

# create output
message(' - creating output data frame.')
df$`NA` = ""
cols = c('node_id','node_name','ontology','size','locus_count','top_pvalue','top_fwer','top_refinedFWER')
for (i in 1:length(traits)) {
    cols = c(cols,'NA',paste0(traits[i],c('_rank','_meanRank','_raw_p_low_rank','_FWER_low_rank','_refinedFWER_low_rank')))
}
output = df[,cols]

# write results
message(sprintf(' - writing %s',outFile))
write.table(output, file = outFile, row.names = F, quote = F, sep = '\t', na = "-")
system(paste0('chmod 770 ', outFile))
message('-- Completed: Combine results of GOfuncR gene set enrichment tests ---')
