#!/usr/bin/env Rscript

# =======================================================
# === Combine prs prediction accuracies across traits ===
# =======================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop(paste0('expected 3 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
traits=args[1] # traits="gap_gm,gap_wm,gap_gwm"
prsFiles=args[2] # prsFiles="results/gap_gm/replicate/EUR/prs.assoc.txt,results/gap_wm/replicate/EUR/prs.assoc.txt,results/gap_gwm/replicate/EUR/prs.assoc.txt"
outFile=args[3] # outFile="results/combined/replicate.prs.txt"

message(paste0('\n--- Combine prs prediction accuracies across traits ---',
               '\ntraits: ', traits,
               '\nprsFiles: ', prsFiles,
               '\noutFile: ', outFile,'\n'))

# attach packages to current R session
for (pkg in c('dplyr','stringr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))}

# transform variables
traits = str_split(traits, ',')[[1]]
prsFiles = str_split(prsFiles, ',')[[1]]

# open and merge datasets
for (i in 1:length(traits)) {
  
  # open datasets and keep relevant variables
  message(sprintf('[%d/%d] Loading %s',i,length(traits),prsFiles[i]))
  tmp = read.delim(prsFiles[i], sep = '\t', header = T)
  tmp = tmp[,c('pthresh', 'R2', 'estimate', 'statistic', 'df', 'p.value')]
  
  # create unique variable names
  names(tmp)[which(names(tmp) %in% c('estimate', 'R2', 'statistic', 'df', 'p.value'))] = 
    paste0(traits[i],c('_R2', '_estimate', '_statistic', '_df', '_p.value'))
  
  # merge datasets
  if (i ==1) { df = tmp } else { df = left_join(df,tmp, by = 'pthresh') }
}

# convert to numbers
df$pthresh = df$pthresh %>%
  sub(pattern = 'Pt_', replacement = '') %>%
  sub(pattern = 'e.', replacement = 'e-') %>%
  as.numeric()

# create output
message(sprintf('Writing %s',outFile))
df$`NA` = ""
cols = c('pthresh')
for (i in 1:length(traits)) {
  cols = c(cols,'NA',paste0(traits[i],c('_R2', '_estimate', '_statistic', '_df', '_p.value')))
}
output = df[,cols]
write.table(output, outFile, sep = '\t', quote = F, row.names = F)
system(sprintf('chmod 770 %s', outFile))
message('--- Completed: Combine prs prediction accuracies across traits ---')
