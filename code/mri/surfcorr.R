#!/usr/bin/env Rscript

# ============================================================
# === Calculate aparc surface correlations for enigma plot ===
# ============================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=6) {
  stop(paste0('expected 6 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
traitName = args[1] # traitName = 'gap_gm'
traitFile = args[2] # traitFile = 'data/traits/gap_gm.txt'
covsFile = args[3] # covsFile = 'data/traits/covs.txt'
covNames = args[4] # covNames = 'sex,age,age2,ac1,ac2,TIV'
freesurferFile = args[5] # freesurferFile = 'results/mri/freesurfer.tsv.gz'
targetDir = args[6] # targetDir = 'results/gap_gm/discovery/surfplot'

message(paste0('\n--- Calculate aparc surface correlations ---',
	           '\ntraitName: ', traitName,
               '\ntraitFile: ', traitFile,
               '\ncovsFile: ', covsFile,
               '\ncovNames: ', covNames,
               '\nfreesurferFile: ', freesurferFile,
               '\ntargetDir: ', targetDir,'\n'))

# attach packages to current R session
for (pkg in c('data.table','dplyr','stringr','ppcor','plotly','psych','reshape2')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# load data
freesurfer = data.frame(fread(cmd = paste0('gzip -dc ', getwd(), '/', freesurferFile), sep = '\t', header = T))
trait = read.delim(traitFile, sep = '\t', header = T)
covs = read.delim(covsFile, sep = '\t', header = T)

# transform variables
covNames = str_split(covNames, ',')[[1]]

# merge trait, freesurfer, and covs information
df = trait %>% inner_join(freesurfer, by = 'IID') %>% inner_join(covs[,c('IID',covNames)], by = 'IID')

# define function to compute correlation matrix
makecorr = function(df, x = NULL, xlabels = x, covNames = NULL, colvars = NULL, type = 'pearson', cluster = F) {  
  
  # prepare data.frame
  if (!is.null(covNames)) { covs = df[,covNames] }
  df = df[,x]
  names(df) = xlabels
  
  # compute rho
  message('Calculating correlations.')
  if (is.null(covNames)) { 
  	rho = cor(x = df, use='pairwise.complete.obs', method=type)
  } else {
  	rho = partial.r(data = cbind(df,covs), x = names(df), y = covNames, use = 'pairwise', method = type)
  }

  # get n
  message('Getting pairwise number of observations.')
  n = data.frame(matrix(data = NA, nrow = length(df), ncol = length(df)))
  names(n) = names(df)
  pb = txtProgressBar(min = 1, max = ncol(n)*(ncol(n)-1)/2, style = 3, title = 'Calculating Correlations.')
  k = 0
  for (i in 1:(length(df)-1)) {
    for (j in (i+1):length(df)) {
    	k = k+1
	  	setTxtProgressBar(pb, k)
	    n[i,j] = n[j,i] = sum(!is.na(df[,i]) & !is.na(df[,j]))
    }
  }
  close(pb)

  # calculate p
  message('Calculating p-values.')
  corrs = corr.p(as.matrix(rho),as.matrix(n)-length(covNames),ci = F,adjust="none")

  # convert to data.frame
  df.rho = data.frame(corrs$r)
  df.pval = data.frame(corrs$p)
  df.n = data.frame(corrs$n)+length(covNames)
  for (i in 1:(length(df))) {
    df.rho[i,i] = 1 
    df.pval[i,i] = 0
    df.n[i,i] = sum(!is.na(df[,i]))
  }
 
  # add id variable
  df.rho$id = df.pval$id = df.n$id = names(df)
  df.rho = df.rho[,c(ncol(df.rho),1:(ncol(df.rho)-1))]
  df.pval = df.pval[,c(length(df.pval),1:(length(df.pval)-1))]
  df.n = df.n[,c(length(df.n),1:(length(df.n)-1))]

  # only keep certain columns
  if (!is.null(colvars)) { 
    df.rho = df.rho[!(x %in% colvars), c(TRUE, x %in% colvars)]
    df.pval = df.pval[!(x %in% colvars), c(TRUE, x %in% colvars)]
    df.n = df.n[!(x %in% colvars), c(TRUE, x %in% colvars)]
    df = df[, !(x %in% colvars)]
  }
  row.names(df.rho) = row.names(df.pval) = row.names(df.n) = NULL

  # return results (rho + pval + n + plot)
  results = list("rho" = df.rho, "pval" = df.pval, "n" = df.n)
  results
  
}

# calculate correlations
corrs = makecorr(df = df, x = c(traitName,names(freesurfer)[4:ncol(freesurfer)]), covNames = covNames, colvars = traitName, type = 'pearson')

# save results
message(sprintf('Saving %s/surfcorr.txt',targetDir))
output = corrs$rho %>% left_join(corrs$pval, by = 'id') %>% left_join(corrs$n, by = 'id') 
names(output) = c('id','rho','pval','n')
system(paste0('mkdir -p ', targetDir))
write.table(output, sprintf('%s/surfcorr.txt',targetDir), sep = '\t', quote = F, col.names = T, row.names = F)
message('\n--- Completed: calculate aparc surface correlations ---')
