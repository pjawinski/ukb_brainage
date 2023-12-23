#!/usr/bin/env Rscript

# ==================================================
# === get exact p-value of heritability estimate ===
# ==================================================

# get arguments from command line
args = commandArgs(trailingOnly=TRUE)

# set arguments
targetDir = args[1]

# calculate exact p-value
message(paste0('\n--- Calculating p-value of heritability estimate ---'))
df = read.delim(paste0(targetDir, '/greml.hsq'))
LRT = df$Variance[df$Source=='LRT']
Pval = 0.5 * pchisq(LRT, df=1, lower.tail=FALSE)
if (Pval == 0) Pval = "< 5E-324"
df$Variance[df$Source=='Pval'] = Pval

# write output
write.table(df, file = paste0(targetDir,'/greml.pval.hsq'), sep = '\t', quote = F, row.names = F, na = "")
message(paste0('Done.'))
