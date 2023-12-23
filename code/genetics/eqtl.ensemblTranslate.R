#!/usr/bin/env Rscript

# ============================================
# === translate ensembl to hgnc and entrez ===
# ============================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Input file must be supplied", call.=FALSE)
}

# set arguments
inputFile=args[1] # inputFile="results/gap_gm/eqtl/eqtl.multiTissue.txt" # inputFile="results/gap_gm/eqtl/eqtl.singleTissue.txt"
outputFile=args[2] # outputFile="results/gap_gm/eqtl/eqtl.multiTissue.hgnc.txt" # outputFile="results/gap_gm/eqtl/eqtl.singleTissue.transl.txt"
columnName=args[3] # columnName="gene_id"
message(paste0('\n--- Translate Ensembl to HGNC IDs | Settings ---',
               '\ninputFile: ', inputFile,
               '\nouputFile: ', outputFile,
               '\ncolumnName: ', columnName,'\n'))

# attach packages to current R session
for (pkg in c('biomaRt', 'dplyr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# load data and add ensembl gene variable without gene version
message('Loading data.')
df = read.delim(inputFile, sep = '\t', header = T)
df$ensembl_gene_id = gsub("\\.[0-9]*$", "", df[[columnName]])

# get hgnc id and hgnc symbol
attempt = 1
ensembl = NULL
while(is.null(ensembl) & attempt<= 10) {
  message(paste0('Connecting to Ensembl mysql database (attempt ', attempt, '/10).'))
  attempt = attempt + 1
  try(ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"))
  if (!is.null(ensembl)) { message('Connection established.') }
} 

if (is.null(ensembl)) {stop("Connection failed. Try again later.", call.=FALSE)}
message('Converting Ensembl gene IDs to HGNC IDs and symbols.')
transl.part1 = getBM(filters = "ensembl_gene_id",
               attributes = c("ensembl_gene_id","hgnc_id", 'hgnc_symbol'),
               values = unique(df$ensembl_gene_id), 
               mart = ensembl)
transl.part2 = getBM(filters = "ensembl_gene_id",
                     attributes = c("ensembl_gene_id","entrezgene_id", 'entrezgene_description'),
                     values = unique(df$ensembl_gene_id), 
                     mart = ensembl)

# only keep first match (sorted by hgnc_id / entrez_id)
transl.part1 = transl.part1 %>% 
               mutate(hgnc.num = as.numeric(gsub(".*:","",hgnc_id))) %>%
               arrange(hgnc.num) %>%
               distinct(ensembl_gene_id, .keep_all = TRUE) %>%
               select(-c(hgnc.num))
transl.part2 = transl.part2 %>% 
               arrange(entrezgene_id) %>%
               distinct(ensembl_gene_id, .keep_all = TRUE)

# add hgnc id and symbol to data frame
message('Writing results.')
output = df %>% left_join(transl.part1, by = 'ensembl_gene_id') %>%
                left_join(transl.part2, by = 'ensembl_gene_id')
write.table(output, file = outputFile, row.names = F, quote = F, sep = '\t')

# # --------------------------
# # --- summarize ---
# # --------------------------
# 
# # how many loci show associations
# loci = unique(df$LOCUS_CNT[!is.na(df$CATALOG_DISEASE.TRAIT)])
# nloci = length(loci)
# 
# # show loci with associated traits
# meta_loci = data.frame(matrix(NA, nrow = nloci, ncol = 4))
# names(meta_loci) = c('LOCUS_CNT', 'CYTOBAND', 'LEAD_SNP', 'CATALOG_TRAITS')
# k = 0
# for (locus in loci) { 
#  k = k + 1
#  cytoband = unique(df$CYTOBAND[df$LOCUS_CNT == locus])
#  lead_snp = unique(df$LEAD_SNP[df$LOCUS_CNT == locus])
#  catalog_traits = unique(df$CATALOG_DISEASE.TRAIT[df$LOCUS_CNT == locus & !is.na(df$CATALOG_DISEASE.TRAIT)])
#  catalog_traits = paste(catalog_traits[order(catalog_traits)], collapse = ' | ')
#  meta_loci[k,] = cbind(locus,cytoband,lead_snp,catalog_traits)
# }
# 
# # show traits with associated loci
# traits = unique(df$CATALOG_DISEASE.TRAIT[!is.na(df$CATALOG_DISEASE.TRAIT)])
# traits = traits[order(traits)]
# ntraits = length(traits)
# 
# meta_traits = data.frame(matrix(NA, nrow = ntraits, ncol = 2))
# names(meta_traits) = c('CATALOG_TRAITS', 'LOCI')
# k = 0
# for (trait in traits) { 
#   k = k + 1
#   loci = unique(df$CYTOBAND[df$CATALOG_DISEASE.TRAIT == trait & !is.na(df$CATALOG_DISEASE.TRAIT)])
#   loci = paste(loci[order(loci)], collapse = ' | ')
#   meta_traits[k,] = cbind(trait,loci)
# }
# 
# # save meta info
# targetDir = sub("/[^/]+$", "", filepath)
# write.table(meta_loci, file = paste0(targetDir,'/catalog.by.locus.txt'), row.names = F, quote = F, sep = '\t')
# write.table(meta_traits, file = paste0(targetDir,'/catalog.by.trait.txt'), row.names = F, quote = F, sep = '\t')
