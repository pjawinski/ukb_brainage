#!/usr/bin/env Rscript

# ==========================================
# === get variables for genetic analysis ===
# ==========================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 6) {
  stop(paste0('expected 6 arguments, but ', length(args), ' arguments provided.'), call.=FALSE)
}

# get arguments from command line
subsFile = args[1] # subsFile='results/mri/iid.discovery.txt'
varsFile = args[2] # varsFile='results/genetics/r2022.vars.pca.txt'
idCol = args[3] # idCol='IID'
vars = args[4] # vars='FID,IID,sex,t1.age,t1.age2,t1.ac1,t1.ac2,array,t1.TIV,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20'
varsRename = args[5] # varsRename='FID,IID,sex,age,age2,ac1,ac2,array,TIV,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20'
outFile = args[6] # outFile="data/traits/covs.txt"
message(paste0('\n--- Get subjects and variables of interest for genetic analyses ---',
               '\nsubsFile: ', subsFile,
               '\nvarsFile: ', varsFile,
               '\nidCol: ', idCol,
               '\nvars: ', vars,
               '\nvarsRename: ', varsRename,
               '\noutFile: ', outFile, '\n'))

# load packages
for (pkg in c('dplyr', 'stringr')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# transform variables
vars = str_split(vars, ',')[[1]]
varsRename = str_split(varsRename, ',')[[1]]

# load data
message('Loading data.')
df = read.delim(varsFile, header = T, colClasses = 'character')
subs = read.delim(subsFile, header = T, colClasses = 'character')

# extract subjects of interest
message('Extracting subjects of interest.')
subs = data.frame(subs[,idCol])
names(subs) = idCol
df = left_join(subs, df, by = idCol)

# extract variables of interest and rename them
message('Extracting variables of interest.')
df = df[,vars]
names(df) = varsRename

# remove columns only containing a single value
message('Removing variables containing no more than one unique value.')
idx = sapply(df, function(x) length(unique(x)) >1)
df = df[,idx]

# save output file
message('Saving output file.')
write.table(df, outFile, sep = '\t', row.names = F, col.names = T, quote = F)

message('--- Completed: Get subjects and variables of interest for genetic analyses ---\n')
