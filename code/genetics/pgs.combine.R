#!/usr/bin/env Rscript

# =======================================================
# === Combine pgs prediction accuracies across traits ===
# =======================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop(paste0('expected 3 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
traits=args[1] # traits="gap_gm,gap_wm,gap_gwm"
pgsFiles=args[2] # pgsFiles="results/gap_gm/replicate/EUR/pgs.base32k.target20k.assoc.txt,results/gap_wm/replicate/EUR/pgs.base32k.target20k.assoc.txt,results/gap_gwm/replicate/EUR/pgs.base32k.target20k.assoc.txt"
outFile=args[3] # outFile="results/combined/pgs.base32k.target20k.assoc.txt"

message(paste0('\n--- Combine pgs prediction accuracies across traits ---',
               '\ntraits: ', traits,
               '\npgsFiles: ', pgsFiles,
               '\noutFile: ', outFile,'\n'))

# attach packages to current R session
for (pkg in c('dplyr','stringr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))}

# transform variables
traits = str_split(traits, ',')[[1]]
pgsFiles = str_split(pgsFiles, ',')[[1]]

# open and merge datasets
for (i in 1:length(traits)) {
  
  # open datasets and keep relevant variables
  message(sprintf('[%d/%d] Loading %s',i,length(traits),pgsFiles[i]))
  tmp = read.delim(pgsFiles[i], sep = '\t', header = T)
  tmp = tmp[,c('pgsMethod', 'R2', 'estimate', 'statistic', 'df', 'p.value')]
  
  # create unique variable names
  names(tmp)[which(names(tmp) %in% c('estimate', 'R2', 'statistic', 'df', 'p.value'))] = 
    paste0(traits[i],c('_R2', '_estimate', '_statistic', '_df', '_p.value'))
  
  # merge datasets
  if (i ==1) { df = tmp } else { df = left_join(df,tmp, by = 'pgsMethod') }
}

# convert pgsMethod label
df$pgsMethod = df$pgsMethod %>%
  sub(pattern = 'Pt_', replacement = 'C+T(') %>%
  sub(pattern = 'e[.]', replacement = 'e-')
idx = startsWith(df$pgsMethod,"C+T")
df$pgsMethod[idx] = paste0(df$pgsMethod[idx],")")

# create output
message(sprintf('Writing %s',outFile))
df$`NA` = ""
cols = c('pgsMethod')
for (i in 1:length(traits)) {
  cols = c(cols,'NA',paste0(traits[i],c('_R2', '_estimate', '_statistic', '_df', '_p.value')))
}
output = df[,cols]
write.table(output, outFile, sep = '\t', quote = F, row.names = F)
system(sprintf('chmod 770 %s', outFile))
message('--- Completed: Combine pgs prediction accuracies across traits ---')
