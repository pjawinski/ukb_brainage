#!/usr/bin/env Rscript

# =====================================================
# === Get per-locus-summary of GWAS catalog results ===
# =====================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop(paste0('expected 5 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
inFile=args[1] # inFile="results/pleiofdr/combined/pleio.crosstrait.catalog.merged.txt"
outFile=args[2] # outFile="results/pleiofdr/combined/pleio.crosstrait.eqtl.summary.txt"
locusCol=args[3] # locusCol="discovnum"
leadsnpCol=args[4] # leadsnpCol="leadsnp"
catalogTraitCol=args[5] # catalogTraitCol="catalog_trait"
message(paste0('\n--- Get per-locus-summary of GWAS catalog results | Settings ---',
               '\ninFile: ', inFile,
               '\noutFile: ', outFile,
               '\nlocusCol: ', locusCol,
               '\nleadsnpCol: ', leadsnpCol,
               '\ncatalogTraitCol: ', catalogTraitCol, '\n'))

# load data
message(' - loading data.')
df = read.delim(inFile, sep = '\t', header = T)

# remove rows without catalog match
df = df[!is.na(df[[catalogTraitCol]]),]

# how many loci show associations
loci = unique(df[[locusCol]])
nloci = length(loci)

# show loci with associated traits
meta_loci = data.frame(matrix(NA, nrow = nloci, ncol = 3))
names(meta_loci) = c('locusnum', 'leadsnp', 'catalog_traits')
k = 0
for (locus in loci) { 
 k = k + 1
 leadsnp = unique(df[df[[locusCol]] == locus,leadsnpCol])
 catalog_traits = unique(df[df[[locusCol]] == locus,catalogTraitCol])
 catalog_traits = paste(catalog_traits[order(catalog_traits)], collapse = ' | ')
 meta_loci[k,] = cbind(locus,leadsnp,catalog_traits)
}

# show traits with associated loci
traits = unique(df[[catalogTraitCol]])
traits = traits[order(traits)]
ntraits = length(traits)

meta_traits = data.frame(matrix(NA, nrow = ntraits, ncol = 2))
names(meta_traits) = c('catalog_trait', 'loci')
k = 0
for (trait in traits) { 
  k = k + 1
  loci = unique(df[df[[catalogTraitCol]]==trait,locusCol])
  loci = paste(loci[order(loci)], collapse = ' | ')
  meta_traits[k,] = cbind(trait,loci)
}

# save meta info
message(sprintf(' - writing %s.by.locus.txt',outFile))
write.table(meta_loci, file = sprintf('%s.by.locus.txt',outFile), row.names = F, quote = F, sep = '\t')
message(sprintf(' - writing %s.by.trait.txt',outFile))
write.table(meta_traits, file = sprintf('%s.by.trait.txt',outFile), row.names = F, quote = F, sep = '\t')
system(sprintf('chmod 770 %s.by.locus.txt',outFile))
system(sprintf('chmod 770 %s.by.trait.txt',outFile))
message('\n--- Completed: Get per-locus-summary of GWAS catalog results ---')
