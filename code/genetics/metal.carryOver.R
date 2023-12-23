#!/usr/bin/env Rscript

# =======================================================
# === renaming and reordering columns of metal output ===
# =======================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop(paste0('expected 2 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# rename and change order of metal output columns
inputFile = args[1] # inputFile="results/pwr_cz_alpha/metal/metal.ivweight.1.txt.tmp"
outputFile = args[2] # outputFile="results/pwr_cz_alpha/metal/metal.ivweight.gz"

# load data
message(' - loading data.')
df = data.frame(data.table::fread(inputFile, header = T, tmpdir=getwd()))

# rename columns
message(' - renaming columns.')
names(df)[names(df)=='MarkerName'] = 'ID'
names(df)[names(df)=='Allele1'] = 'A1'
names(df)[names(df)=='Allele2'] = 'A2'
names(df)[names(df)=='Freq1'] = 'A1_FREQ'
names(df)[names(df)=='FreqSE'] = 'A1_FREQ_SE'
names(df)[names(df)=='MinFreq'] = 'A1_FREQ_MIN'
names(df)[names(df)=='MaxFreq'] = 'A1_FREQ_MAX'
names(df)[names(df)=='Effect'] = 'BETA'
names(df)[names(df)=='StdErr'] = 'SE'
names(df)[names(df)=='Zscore'] = 'Z'
names(df)[names(df)=='P.value'] = 'P'
names(df)[names(df)=='Direction'] = 'DIRECTION'
names(df)[names(df)=='HetISq'] = 'HET_ISQ'
names(df)[names(df)=='HetCHiSq'] = 'HET_CHISQ'
names(df)[names(df)=='HetDf'] = 'HET_DF'
names(df)[names(df)=='HetPVal'] = 'HET_P'
names(df)[names(df)=='Ntotal'] = 'N'

# change column order
message(' - changing column order.')
cols_ivweight = c('ID','CHR','BP','A1','A2','A1_FREQ','A1_FREQ_SE','A1_FREQ_MIN','A1_FREQ_MAX','BETA','SE','P','N','DIRECTION','HET_ISQ','HetChiSq','HET_DF','HET_P')
cols_nweight = c('ID','CHR','BP','A1','A2','A1_FREQ','A1_FREQ_SE','A1_FREQ_MIN','A1_FREQ_MAX','Z','P','N','DIRECTION','HET_ISQ','HetChiSq','HET_DF','HET_P')
if (sum(!(cols_ivweight %in% names(df))) == 0) {
  df = df[,c(cols_ivweight,names(df)[which(!(names(df) %in% cols_ivweight))])]
} else if (sum(!(cols_nweight %in% names(df))) == 0) {
  df = df[,c(cols_nweight,names(df)[which(!(names(df) %in% cols_nweight))])]
}

# change alleles to upper case
df$A1 = toupper(df$A1)
df$A2 = toupper(df$A2)

# drop Weight
if ('Weight' %in% names(df)) {
  df$N = df$Weight
  df = df[,-which(names(df)=='Weight')]
}

# write output
message(' - writing file.')
data.table::fwrite(df, file = outputFile, sep = ' ', compress = 'gzip')
system(paste0('chmod -R 770 ', outputFile))

