#!/usr/bin/env Rscript

# ==================================================
# === Get credible sets using SusieR finemapping ===
# ==================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=6) {
  stop(paste0('expected 6 arguments, but ', length(args), ' arguments provided.'), call.=FALSE)
}

# get arguments from command line
ldfile = args[1] # ldfile="/slow/projects/ukb_brainage/results/gap_wm/gwama/eur/susieR/1_212557599/ldstore.ld"
zfile = args[2] # zfile="/slow/projects/ukb_brainage/results/gap_wm/gwama/eur/susieR/1_212557599/ldstore.z"
leadsnp = args[3] # leadsnp="rs3767867"
n = as.numeric(args[4]) # n=54890
outputDF = args[5] # outputDF="/slow/projects/ukb_brainage/results/gap_wm/gwama/eur/susieR/1_212557599/susier.1_212557599.DF.txt"
outputLS = args[6] # outputLS="/slow/projects/ukb_brainage/results/gap_wm/gwama/eur/susieR/1_212557599/susier.1_212557599.LS.txt"

message('\n--- Get credible set of variants using SusieR ---',
 '\nldfile: ', ldfile,
 '\nzfile: ', zfile,
 '\nleadsnp: ', leadsnp,
 '\nn: ', n,
 '\noutputDF: ', outputDF,
 '\noutputLS: ', outputLS,'\n')

# attach required packages and create target directory
for (pkg in c('bigsnpr','data.table','dplyr','susieR')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# import data
message('Importing data.')
ld = data.table::fread(ldfile, header = F, sep = ' ')
ld = unname(as.matrix(ld))
z = read.delim(zfile, header = T, sep = ' ')

# identify variants with missing r2 values
if (any(is.na(ld))) {
  message('Removing variants with missing r2 values.')
  row_na_counts = rowSums(is.na(ld))
  row_na_ordered = order(-row_na_counts)
  k = 1
  idx = row_na_ordered[k]
  while (any(is.na(ld[-idx,-idx]))) {
    k = k + 1
    idx = c(idx,row_na_ordered[k])
  }
  outputR2missing = ifelse(grepl("\\.txt$", outputLS), 
                     gsub(".txt$", ".removed.txt", outputLS), 
                     paste0(outputLS, ".removed.txt"))
  message(sprintf('Writing %s', outputR2missing))
  data.table::fwrite(z[idx,], file = outputR2missing, compress = 'none', sep = '\t', quote = F, row.names = F, col.names = T)
  ld = ld[-idx,-idx]
  z = z[-idx,]
}

# add R2 with index snp
z$r2 = ld[,grep(leadsnp,z$rsid)]^2

# get credible set of variants
message('Getting credible set of variants.')
fit = susie_rss(bhat = z$beta, shat = z$se, R = ld, n = n, L = 10, refine = T, n_purity = 1000)
cs = susie_get_cs(fit, Xcorr = ld, min_abs_corr = 0.0, n_purity = 1000)

for (i in 1:max(length(fit$sets$cs),1)) {
  idx = cs$cs[[i]]
  purity = cs$purity[i,1:3] %>% unlist() %>% as.vector()
  tmp = z[idx,]
  tmp$index_rsid = leadsnp
  tmp$cs = i
  tmp$cs.size = length(idx)
  tmp$coverage = cs$coverage[i]
  tmp$purity.min = purity[1]
  tmp$purity.median = purity[3]
  tmp$pip = fit$pip[idx]
  tmp = tmp[order(tmp$pip, decreasing = T),]
  tmp$cumpip = cumsum(tmp$pip)
  if (i==1) { df = tmp } else { df = rbind(df, tmp) }
}

# add snpTotal, snpCausal, regionalh2, minBP, maxBP
df$snpTotal = nrow(ld)
df$snpCausal = max(df$cs)
df$minBP = min(z$position)  
df$maxBP = max(z$position)
df$rangeBP = max(z$position) - min(z$position) 

# get strings of concatenated credible set sizes, and credible snps
cssize = paste(df$cs.size[!duplicated(df$cs)],collapse = ' | ')
minpurity = df$purity.min[!duplicated(df$cs)] %>% sprintf("%.2f", .) %>% paste(collapse = " | ")
cssnps = df %>%
  group_by(cs) %>%
  summarise(y = paste0(rsid, collapse = ",")) %>%
  summarise(y = paste(y, collapse = " | ")) %>%
  pull(y)

# sort columns
DF = df[,c('index_rsid','chromosome','minBP','maxBP','rangeBP','snpTotal','snpCausal','cs','cs.size','coverage','purity.min','purity.median','rsid','position','allele1','allele2','maf','beta','se','r2','pip','cumpip')]
LS = cbind(df[1,c('index_rsid','chromosome','minBP','maxBP','rangeBP','snpTotal','snpCausal')], data.frame(minpurity = minpurity, cssizes = cssize, cssnps = cssnps))

# output
message(sprintf('Writing %s', outputDF))
data.table::fwrite(DF, file = outputDF, compress = 'none', sep = '\t', quote = F, row.names = F, col.names = T)
message(sprintf('Writing %s', outputLS))
data.table::fwrite(LS, file = outputLS, compress = 'none', sep = '\t', quote = F, row.names = F, col.names = T)
system(sprintf('chmod -R 770 %s', outputDF))
system(sprintf('chmod -R 770 %s', outputLS))
message('--- Completed: Get credible set of variants using SusieR ---')
