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
inputFile=args[1] # inputFile="results/pleiofdr/combined/pleio.crosstrait.eqtl.singleTissue.txt"
outputFile=args[2] # outputFile="results/pleiofdr/combined/pleio.crosstrait.eqtl.singleTissue.txt.tmp"
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
output$hgnc_id[is.na(output$hgnc_id)]=""
output$hgnc_symbol[is.na(output$hgnc_symbol)]=""
output$entrezgene_id[is.na(output$entrezgene_id)]=""
output$entrezgene_description[is.na(output$entrezgene_description)]=""
write.table(output, file = outputFile, row.names = F, quote = F, sep = '\t')
