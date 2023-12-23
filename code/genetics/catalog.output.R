#!/usr/bin/env Rscript

# ======================================
# === Summarise GWAS catalog results ===
# ======================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Input file must be supplied", call.=FALSE)
}

filepath=args[1]

# select columns
df = read.delim(filepath, sep = '\t', header = T)
df = df[,c('LOCUS_CNT','SNP_CNT','CHR','BP','ID','CYTOBAND','LEAD_SNP','LEAD_SNP_P','LEAD_SNP_KB','LEAD_SNP_RSQ','LEAD_SNP_ALLELES','REGION','NEAREST_GENE','NEAREST_GENE_DESCRIPTION','NEAREST_GENE_BIOTYPE','DISTANCE','TYPED','INFO','A1','A2','A1_FREQ','BETA','SE','Z','P','ETA2','N','CATALOG_DISEASE.TRAIT','CATALOG_PUBMEDID','CATALOG_FIRST_AUTHOR','CATALOG_DATE','CATALOG_RISK_ALLELE','CATALOG_P.VALUE')]

# remove rows without catalog match
df = df[!is.na(df$CATALOG_DISEASE.TRAIT),]

# save table
write.table(df, file = filepath, row.names = F, quote = F, sep = '\t')

# --------------------------
# --- get some meta info ---
# --------------------------

# how many loci show associations
loci = unique(df$LOCUS_CNT[!is.na(df$CATALOG_DISEASE.TRAIT)])
nloci = length(loci)

# show loci with associated traits
meta_loci = data.frame(matrix(NA, nrow = nloci, ncol = 4))
names(meta_loci) = c('LOCUS_CNT', 'CYTOBAND', 'LEAD_SNP', 'CATALOG_TRAITS')
k = 0
for (locus in loci) { 
 k = k + 1
 cytoband = unique(df$CYTOBAND[df$LOCUS_CNT == locus])
 lead_snp = unique(df$LEAD_SNP[df$LOCUS_CNT == locus])
 catalog_traits = unique(df$CATALOG_DISEASE.TRAIT[df$LOCUS_CNT == locus & !is.na(df$CATALOG_DISEASE.TRAIT)])
 catalog_traits = paste(catalog_traits[order(catalog_traits)], collapse = ' | ')
 meta_loci[k,] = cbind(locus,cytoband,lead_snp,catalog_traits)
}

# show traits with associated loci
traits = unique(df$CATALOG_DISEASE.TRAIT[!is.na(df$CATALOG_DISEASE.TRAIT)])
traits = traits[order(traits)]
ntraits = length(traits)

meta_traits = data.frame(matrix(NA, nrow = ntraits, ncol = 2))
names(meta_traits) = c('CATALOG_TRAITS', 'LOCI')
k = 0
for (trait in traits) { 
  k = k + 1
  loci = unique(df$CYTOBAND[df$CATALOG_DISEASE.TRAIT == trait & !is.na(df$CATALOG_DISEASE.TRAIT)])
  loci = paste(loci[order(loci)], collapse = ' | ')
  meta_traits[k,] = cbind(trait,loci)
}

# save meta info
targetDir = sub("/[^/]+$", "", filepath)
write.table(meta_loci, file = paste0(targetDir,'/catalog.by.locus.txt'), row.names = F, quote = F, sep = '\t')
write.table(meta_traits, file = paste0(targetDir,'/catalog.by.trait.txt'), row.names = F, quote = F, sep = '\t')
