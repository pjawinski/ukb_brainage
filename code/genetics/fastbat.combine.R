#!/usr/bin/env Rscript

# ==================================================
# === Combine GCTA-fastBAT results across traits ===
# ==================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop(paste0('expected 4 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
traits = args[1] # traits="gap_gm,gap_wm,gap_gwm"
inputFiles = args[2] # inputFiles="results/gap_gm/fastbat/fastbat.0kb.gene.fastbat.annovar.clumped,results/gap_wm/fastbat/fastbat.0kb.gene.fastbat.annovar.clumped,results/gap_gwm/fastbat/fastbat.0kb.gene.fastbat.annovar.clumped"
clumpingWindow = as.numeric(args[3]) * 1000 # clumpingWindow = as.numeric(args[3]) * 1000 # clumpingWindow = 10000 * 1000
outputFile = args[4] # outputFile="results/combined/fastbat.txt"

message(paste0('\n--- Combine GCTA-fastBAT results across traits ---',
               '\ntraits: ', traits,
               '\ninputFiles: ', inputFiles,
               '\nclumpingWindow: ', clumpingWindow,
               '\noutputFile: ', outputFile,'\n'))

# attach required packages
for (pkg in c('dplyr','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# transform variables
traits = str_split(traits, ',')[[1]]
inputFiles = str_split(inputFiles, ',')[[1]]

# open and merge datasets
for (i in 1:length(traits)) {
  
  # open datasets and keep relevant variables
  message(sprintf(' - reading input file %s',inputFiles[i]))
  tmp = read.delim(inputFiles[i], sep = '\t', header = T)
  if (i ==1) {
    tmp = tmp[,c('Gene','GENE_DESCRIPTION','GENE_BIOTYPE','ENTREZ_ID','HGNC_ID','CYTOBAND','Chr','Start','End','No.SNPs',
               'Chisq.Obs.', 'Pvalue', 'FDR','LOCUS_COUNT', 'GENE_COUNT')]
  } else { 
    tmp = tmp[,c('Gene','Chisq.Obs.','Pvalue','FDR','LOCUS_COUNT','GENE_COUNT')]
  }
  
  # create unique variable names
  names(tmp)[which(names(tmp) %in% c('Chisq.Obs.','Pvalue','FDR','LOCUS_COUNT','GENE_COUNT'))] = 
    paste0(traits[i],c('_Chisq.Obs.','_Pvalue','_FDR','_LOCUS_COUNT','_GENE_COUNT'))
  
  # merge datasets
  if (i ==1) { df = tmp } else { df = left_join(df,tmp, by = 'Gene') }
}
  
# get top p-value and top |rho|
df$topPvalue =  df[,grep('Pvalue',names(df))] %>% apply(1, FUN = min)
df$topFDR = df[,grep('FDR',names(df))] %>% apply(1, FUN = min)

# define function for clumping
clumpGenes <- function(df, pval = 'top_pvalue', gene = 'Gene', chr = 'Chr', start = 'Start', end = 'End', clumpingWindow = 3000000) {
   
  # create variables
  df$INDEX_CLUMP = 999
  df$INDEX_GENE = "None"
  df$INDEX_GENE_PVAL = 999
  clump = 1

  # loop across chromosomes
  for (i in names(table(df[,chr]))) {
    message(sprintf(' - clumping genes on chr%s...', i))
    df_sub = df[df[,chr]==i,]
    df_sub = df_sub[order(df_sub[,pval]),]
    
    # loop across genes
    for (j in 1:nrow(df_sub)) {
      if (df_sub$INDEX_CLUMP[j]==999) {
        gene1 = df_sub[j,]
        df_sub$INDEX_CLUMP[j] = clump
        df_sub$INDEX_GENE[j] = gene1[,gene]
        df_sub$INDEX_GENE_PVAL[j] = gene1[,pval]
        
        # get distance with other genes amd clump if necessary
        for (k in (j+1):nrow(df_sub)) {
          if (k > nrow(df_sub)) {
            break
          }
          if (df_sub$INDEX_CLUMP[k]==999) {
            gene2 = df_sub[k,]
            distance = min(abs(c(gene1[,start]-gene2[,start], gene1[,start]-gene2[,end], gene1[,end]-gene2[,start], gene1[,end]-gene2[,start])))
            if (distance < clumpingWindow) {
              df_sub$INDEX_CLUMP[k] = clump
              df_sub$INDEX_GENE[k] = gene1[,gene]
              df_sub$INDEX_GENE_PVAL[k] = gene1[,pval]
            }
          }
        }
        clump = clump + 1  
      }
    }
    idx = as.vector(sapply(df_sub[,gene], function(x) which(x == df[,gene])))
    df$INDEX_CLUMP[idx] = df_sub$INDEX_CLUMP
    df$INDEX_GENE[idx] = df_sub$INDEX_GENE
    df$INDEX_GENE_PVAL[idx] = df_sub$INDEX_GENE_PVAL
  }

  # add locus count and gene_count
  df = df[order(df$INDEX_GENE_PVAL, df[,pval], df[,gene]),]
  df$LOCUS_COUNT = -1
  df$GENE_COUNT = -1

  index_gene = "None"
  locus_count = 0
  gene_count = 0
  for (i in 1:nrow(df)){
    gene_count = gene_count + 1
    if (index_gene != df$INDEX_GENE[i]) {
      index_gene = df$INDEX_GENE[i]
      locus_count = locus_count + 1
      gene_count = 1
    }
    df$LOCUS_COUNT[i] = locus_count
    df$GENE_COUNT[i] = gene_count
  } 
  return(df)
}

# do 2nd-level clumping based on top-pvalue
message(' - starting clumping procedure.')
clumped = clumpGenes(df, pval = 'topPvalue', gene = 'Gene', chr = 'Chr', start = 'Start', end = 'End', clumpingWindow = clumpingWindow)

# add cross-trait clumps
df = left_join(df,clumped[,c('Gene','LOCUS_COUNT','GENE_COUNT')],by = 'Gene')

# create output
df$`NA` = ""
cols = c('Gene','GENE_DESCRIPTION','topPvalue','topFDR')
for (i in 1:length(traits)) {
  cols = c(cols,'NA',paste0(traits[i],c('_Chisq.Obs.','_Pvalue','_FDR','_LOCUS_COUNT','_GENE_COUNT')))
}
cols = c(cols, c('NA','ENTREZ_ID','HGNC_ID','CYTOBAND','Chr','Start','End','No.SNPs','LOCUS_COUNT','GENE_COUNT'))
output = df[,cols]
output = output[order(output$topFDR,output$topPvalue),]

# write file
message(sprintf(' - writing %s',outputFile))
write.table(output, outputFile, sep = '\t', quote = F, row.names = F)
message('--- Completed: Combine GCTA-fastBAT results across traits ---')


