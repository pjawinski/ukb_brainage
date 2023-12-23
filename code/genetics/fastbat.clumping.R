#!/usr/bin/env Rscript

# ===============================
# === clump genes by position ===
# ===============================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop(paste0('expected 2 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# get arguments from command line
targetFile = args[1] # targetFile="results/gap_gwm/fastbat/fastbat.0kb.gene.fastbat.annovar"
clumpingWindow = as.numeric(args[2])*1000 # clumpingWindow = 10000*1000

# Start clumping
message(paste0('\n--- Clumping gene-based results ---',
               '\ntargetFile: ', targetFile,'\n'))
message(paste0('Loading ', targetFile))

# import sumstats
df = read.delim(targetFile, sep='\t', header=T, stringsAsFactors=FALSE)
df$FDR = p.adjust(df$Pvalue, method = 'BH', n = length(df$Pvalue))

# do clumping
message('Starting clumping procedure.')
df$INDEX_CLUMP = 999
df$INDEX_GENE = "None"
df$INDEX_GENE_PVAL = 999
clump = 1

for (i in names(table(df$Chr))) {
  message(paste0('Clumping genes on Chr', i, '...'))
  df_sub = df[df$Chr==i,]
  df_sub = df_sub[order(df_sub$Pvalue),]
  
  for (j in 1:nrow(df_sub)) {
    if (df_sub$INDEX_CLUMP[j]==999) {
      gene1 = df_sub[j,]
      df_sub$INDEX_CLUMP[j] = clump
      df_sub$INDEX_GENE[j] = gene1$Gene
      df_sub$INDEX_GENE_PVAL[j] = gene1$Pvalue
      
      for (k in (j+1):nrow(df_sub)) {
        if (k > nrow(df_sub)) {
          break }
        
        if (df_sub$INDEX_CLUMP[k]==999) {
          gene2 = df_sub[k,]
          distance = min(c(abs(gene1$Start - gene2$Start),abs(gene1$Start - gene2$End), abs(gene1$End - gene2$Start),abs(gene1$End - gene2$End)))
          if (distance < clumpingWindow) {
            df_sub$INDEX_CLUMP[k] = clump
            df_sub$INDEX_GENE[k] = gene1$Gene
            df_sub$INDEX_GENE_PVAL[k] = gene1$Pvalue
          }
        }
      }
      clump = clump + 1  
    }
  }
  index = as.vector(sapply(df_sub$Gene, function(x) which(x == df$Gene)))
  df$INDEX_CLUMP[index] = df_sub$INDEX_CLUMP
  df$INDEX_GENE[index] = df_sub$INDEX_GENE
  df$INDEX_GENE_PVAL[index] = df_sub$INDEX_GENE_PVAL
}

# add locus count and gene_count
df = df[order(df$INDEX_GENE_PVAL, df$Pvalue, df$Gene),]
df$LOCUS_COUNT = -1
df$GENE_COUNT = -1

index_gene = "None"
locus_count = 0
gene_count = 0

for (i in 1:length(df[,1])){
  gene_count = gene_count + 1
  if (index_gene != df$INDEX_GENE[i]) {
    index_gene = df$INDEX_GENE[i]
    locus_count = locus_count + 1
    gene_count = 1
  }
  df$LOCUS_COUNT[i] = locus_count
  df$GENE_COUNT[i] = gene_count
}

# sort 
df = df[,c('LOCUS_COUNT','GENE_COUNT','Gene', 'GENE_DESCRIPTION','GENE_BIOTYPE', 'ENTREZ_ID', 'HGNC_ID', 'CYTOBAND', 'Chr', 'Start','End','No.SNPs','SNP_start','SNP_end','Chisq.Obs.','Pvalue','FDR','TopSNP','TopSNP.Pvalue','INDEX_GENE','INDEX_GENE_PVAL')]

# write file
write.table(df, file=paste0(targetFile, '.clumped'), quote = FALSE, row.names = FALSE, sep = '\t')

# output some summary statistics
output = data.frame(geneCount_total = nrow(df),
                    geneCount_FDR = sum(df$FDR < 0.05),
                    geneCount_Bonf = sum(df$Pvalue * nrow(df) < 0.05),
                    locusCount_total = max(df$LOCUS_COUNT),
                    locusCount_FDR = sum(df$GENE_COUNT == 1 & df$FDR < 0.05),
                    locusCount_Bonf = sum(df$GENE_COUNT == 1 & df$Pvalue * nrow(df) < 0.05))
write.table(output, file=paste0(targetFile, '.clumped.summary'), quote = FALSE, row.names = FALSE, sep = '\t')
message(paste0('--- Clumping gene-based results finished ---\n'))

