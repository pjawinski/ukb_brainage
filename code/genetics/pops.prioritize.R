#!/usr/bin/env Rscript

# ========================
# === Prioritize genes ===
# ========================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=6) {
  stop(paste0('expected 6 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
conditionalFile=args[1] # conditionalFile="results/gap_gm/conditional/conditional.cleaned.txt"
conditionalPthresh=as.numeric(args[2]) # conditionalPthresh=5E-8
popsFile=args[3] # popsFile="results/gap_gm/pops/pops.preds.symbols"
annotationFile=args[4] # annotationFile="/fast/software/pops/data/utils/gene_annot_jun10.txt"
outputFile=args[5] # outputFile="results/gap_gm/pops/pops.prio.txt"
window=as.numeric(args[6])*1000 # window=500*1000

message(paste0('\n--- Prioritize genes from PoPS | Settings ---',
               '\nconditionalFile: ', conditionalFile,
               '\nconditionalPthresh: ', conditionalPthresh,
               '\npopsFile: ', popsFile,
               '\nannotationFile: ', annotationFile,
               '\noutputFile: ', outputFile, '\n'))

# attach packages to current R session
for (pkg in c('dplyr','stringr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# load data
message(' - loading data.')
cond = read.delim(conditionalFile, sep = '\t', header = T, quote = "")
pops = read.delim(popsFile, sep = '\t', header = T, quote = "")
anno = read.delim(annotationFile, sep = '\t', header = T, quote = "")

# filter conditional
message(sprintf(' - filtering conditional signals at pthresh %.2e.',conditionalPthresh))
cond = cond[cond$pJ < conditionalPthresh,]

# merge pops and anno
anno = anno[,c('ENSGID','CHR','START','END')]
genes = left_join(pops,anno,'ENSGID')

# get genes
message(sprintf(' - getting relevant genes within %d kbp', window/1000))
df = cond[,c('SNP','Chr','bp')]
df$prioritized = NA
for (i in 1:nrow(df)) {
  chr = df$Chr[i] 
  bp = df$bp[i]
  tmp.genes = genes[genes$CHR==chr & (abs(genes$START-bp) < window | abs(genes$START-bp) < window) & genes$PoPS_Score > 0,]
  if (nrow(tmp.genes) > 0) {
    prio = data.frame(genes = tmp.genes$SYMBOL, scores = tmp.genes$PoPS_Score)
    prio = prio[order(-prio$scores),]
    prio = prio[1:(min(3,nrow(prio))),]
    df$prioritized[i] = paste0(sprintf('%s (%0.2f)', prio$genes,prio$scores), collapse = ' | ')
    rm(chr,bp,prio,tmp.genes)
  }
}

# write results
message(sprintf(' - writing %s',outputFile))
write.table(df, file = outputFile, row.names = F, quote = F, sep = '\t', na = "NA")
system(sprintf('chmod 750 %s',outputFile))
message('--- Completed: Prioritize genes from PoPS')
