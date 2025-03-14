#!/usr/bin/env Rscript

# ================================
# === rename/rearrange columns ===
# ================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=8) {
  stop(paste0('expected 8 arguments, but ', length(args), ' arguments provided.'), call.=FALSE)
}

# get arguments from command line
inFile = args[1] # inFile="results/pleiofdr/combined/pleio.crosstrait.eqtl.singleTissue.txt"
outFile = args[2] # outFile="results/pleiofdr/combined/pleio.crosstrait.eqtl.singleTissue.suppl.txt"
compress =  args[3] # compress="none"
skipLines = as.numeric(args[4]) # skipLines=0
sep = args[5] # sep="auto"
colsIn = args[6] # colsIn="locusnum,leadsnp,id,variant_id,gene_id,hgnc_symbol,tissue,tss_distance,ma_samples,ma_count,maf,slope,slope_se,pval_nominal,pval_nominal_threshold,min_pval_nominal,pval_beta"
colsRename = args[7] # colsRename="locusnum,leadsnp,id,variant_id,gene_id,hgnc_symbol,tissue,tss_distance,ma_samples,ma_count,maf,slope,slope_se,pval_nominal,pval_nominal_threshold,min_pval_nominal,pval_beta"
colsOut = args[8] # colsOut="locusnum,leadsnp,id,variant_id,gene_id,hgnc_symbol,NA,tissue,tss_distance,ma_samples,ma_count,NA,maf,slope,slope_se,pval_nominal,pval_nominal_threshold,min_pval_nominal,pval_beta"

message('\n--- Rename/rearrange columns | Settings ---',
               '\ninFile: ', inFile,
               '\noutFile: ', outFile,
               '\nskipLines: ', skipLines,
               '\nsep: ', sep,
               '\ncolsIn: ', colsIn,
               '\ncolsRename: ', colsRename,
               '\ncolsOut: ', colsOut,'\n')

# transform variables
colsIn = stringr::str_split(colsIn, ',')[[1]]
colsRename = stringr::str_split(colsRename, ',')[[1]]
colsOut = stringr::str_split(colsOut, ',')[[1]]
if (sep == "tab") { sep = "\t" }

# import data
message(paste0(' - importing data.'))
if (stringr::str_sub(inFile,-3,-1) == '.gz') {
  df = data.frame(data.table::fread(cmd=sprintf("gzip -dc %s", inFile), select = colsIn, fill = T, tmpdir = getwd(), header=T, skip = skipLines, sep = sep, stringsAsFactors=FALSE))
} else {
  df = data.frame(data.table::fread(file = inFile, select = colsIn, fill = T, tmpdir = getwd(), header=T, skip = skipLines, sep = sep, stringsAsFactors=FALSE))
}

library(dplyr)



# rename columns
message(paste0(' - renaming columns.'))
names(df) = colsRename

# add Empty columns
message(paste0(' - preparing output dataframe.'))
df$`NA` = NA
df = df[,colsOut]

# write output
message(sprintf(' - writing %s', outFile))
data.table::fwrite(df, file = outFile, sep = '\t', compress = compress)
system(sprintf('chmod -R 770 %s', outFile))
message('--- Completed: Rename/rearrange columns ---')

