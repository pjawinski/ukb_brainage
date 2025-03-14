#!/usr/bin/env Rscript

# ==================================================
# === prepare replication covsFile and traitFile ===
# ==================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=6) {
  stop(paste0('expected 6 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
trait = args[1] # trait="t1.gap_gm"
targetDir = args[2] # targetDir="results/gap_gm/replicate/EUR"
data = args[3] # data="data/replicate/mri_repl_pan.txt"
covs = args[4] # covs="sex,t1.age,t1.age2,t1.ac1,t1.ac2,t1.ac3,t1.TIV,array,PanC1,PanC2,PanC3,PanC4,PanC5,PanC6,PanC7,PanC8,PanC9,PanC10,PanC11,PanC12,PanC13,PanC14,PanC15,PanC16,PanC17,PanC18,PanC19,PanC20"
ancestry = args[5] # ancestry="EUR"
ancestriesCol = args[6] # ancestriesCol="pan"

message(paste0('\n--- Settings for selection of replication cohort ---',
               '\ntrait: ', trait,
               '\ntargetDir: ', targetDir,
               '\ndata: ', data,
               '\ncovs: ', covs,
               '\nancestry: ', ancestry,
               '\nancestriesCol: ', ancestriesCol,'\n'))

# attach packages to current R session
for (pkg in c('dplyr','stringr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# transform variables
covs = str_split(covs, ',')[[1]]

# load data and select ancestry + covs
df = read.delim(data, sep = '\t', header = T)
df = df[df[[ancestriesCol]]==ancestry & !is.na(df[[ancestriesCol]]),c('FID','IID',covs,trait)]

# keep only cols with more than one value
idx = which(sapply(df,function(x) unique(x) %>% length) %>% as.vector() > 1)
df = df[,idx]

# export traitFile
df.trait = df[,c('FID','IID',trait)] 
# names(df.trait) = names(df.trait) %>% sub(pattern='t1.',replacement='') %>% sub(pattern='t2.',replacement='')
write.table(df.trait, sprintf('%s/trait.txt',targetDir), sep = '\t', row.names = F, quote = F)

# export covsFile
df.covs = df[,-which(names(df) %in% trait)]
# names(df.covs) = names(df.covs) %>% sub(pattern='t1.',replacement='') %>% sub(pattern='t2.',replacement='')
write.table(df.covs, sprintf('%s/covs.txt',targetDir),, sep = '\t', row.names = F, quote = F)

