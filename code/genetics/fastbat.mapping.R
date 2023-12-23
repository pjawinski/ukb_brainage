#!/usr/bin/env Rscript

# ===============================================
# === Map gene-level to snp-level discoveries ===
# ===============================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop(paste0('expected 5 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# get arguments from command line
conditionalFile = args[1] # conditionalFile="results/gap_gwm/conditional/conditional.cleaned.tophits.annovar.txt"
fastbatFile = args[2] # fastbatFile = "results/gap_gwm/fastbat/fastbat.0kb.gene.fastbat.annovar.clumped"
pthreshMapping = as.numeric(args[3]) # pthreshMapping = 5E-8
windowSize = as.numeric(args[4]) * 1000 # windowSize = 10000 * 1000
outputPrefix = args[5] # outputPrefix = "results/gap_gwm/fastbat/fastbat.0kb.mapping"

# Start clumping
message(paste0('\n--- Clumping gene-based results ---',
               '\nconditionalFile: ', conditionalFile,
               '\nfastbatFile: ', fastbatFile,
               '\npthreshMapping: ', pthreshMapping,
               '\nwindowSize: ', windowSize,
               '\noutputPrefix: ', outputPrefix,'\n'))

# import data
message('Loading files.')
conditional.raw = read.delim(conditionalFile, sep='\t', header=T)
fastbat.raw = read.delim(fastbatFile, sep='\t', header=T)

# prepare data.frames (only keep relevant results)
conditional = conditional.raw
conditional = conditional[conditional$SNP_CNT == 1 & conditional$P < pthreshMapping & conditional$LEAD_SNP_pJ < pthreshMapping,]
conditional$CHR = as.character(factor(conditional$CHR, c(as.character(1:22),'X','Y','XY','MT'), as.character(1:26)))
fastbat = fastbat.raw
fastbat = fastbat[order(fastbat$Pvalue),]
fastbat$Bonf = as.numeric(fastbat$Pvalue < 0.05/nrow(fastbat))
fastbat = fastbat[fastbat$FDR < 0.05,]

# map gene-level to snp-level discoveries
message('Mapping gene-level to snp-level discoveries.')
mapping = data.frame(matrix(data = NA, nrow = nrow(fastbat), ncol = 5))
names(mapping) = c('Gene', 'IndexGene', 'Bonf', 'nearestSNP', 'nearestSNPdistance')

for (i in 1:nrow(fastbat)) {
  mapping$Gene[i] = fastbat$Gene[i]
  mapping$IndexGene[i] = fastbat$INDEX_GENE[i]
  mapping$Bonf[i] = fastbat$Bonf[i]
  
  # get distance between gene-level discovery and nearest SNP-level discovery
  tmp = conditional[conditional$CHR == fastbat$Chr[i],]
  if (nrow(tmp) > 0) {
    for (j in 1:nrow(tmp)) {
      if (tmp$BP[j] > fastbat$Start[i] & tmp$BP[j] < fastbat$End[i]) {
        mapping$nearestSNP[i] = tmp$ID[j]
        mapping$nearestSNPdistance[i] = 0
        next
      } else {
        distance = min(abs(tmp$BP[j] - fastbat$Start[i]), abs(tmp$BP[j] - fastbat$End[i]))
        if (distance < windowSize & (distance < mapping$nearestSNPdistance[i] | is.na(mapping$nearestSNPdistance[i]))) {
          mapping$nearestSNP[i] = tmp$ID[j]
          mapping$nearestSNPdistance[i] = distance
        }
      }
    }
  }
}
mapping$novelDiscovery = as.numeric(is.na(mapping$nearestSNPdistance))

# get summary of mapped and novel discoveries
mapping.summary = data.frame(nGenes = nrow(mapping),
                             nGenes.mapped = sum(!is.na(mapping$nearestSNP)),
                             nGenes.novel = sum(is.na(mapping$nearestSNP)),
                             nGenesBonf = sum(mapping$Bonf == 1),
                             nGenesBonf.mapped = sum(mapping$Bonf == 1 & !is.na(mapping$nearestSNP)),
                             nGenesBonf.novel = sum(mapping$Bonf == 1 & is.na(mapping$nearestSNP)),
                             nIndexGenes = sum(mapping$Gene == mapping$IndexGene),
                             nIndexGenes.mapped = sum(mapping$Gene == mapping$IndexGene & !is.na(mapping$nearestSNP)),
                             nIndexGenes.novel = sum(mapping$Gene == mapping$IndexGene & is.na(mapping$nearestSNP)),
                             nIndexGenesBonf = sum(mapping$Gene == mapping$IndexGene & mapping$Bonf == 1),
                             nIndexGenesBonf.mapped = sum(mapping$Gene == mapping$IndexGene & mapping$Bonf == 1 & !is.na(mapping$nearestSNP)),
                             nIndexGenesBonf.novel = sum(mapping$Gene == mapping$IndexGene & mapping$Bonf == 1 & is.na(mapping$nearestSNP)))

# rename mapping variables
gene2snp = mapping
gene2snp.summary = mapping.summary

# ===============================================
# === Map snp-level to gene-level discoveries ===
# ===============================================

# (a) map genes annotated at snp-level to gene-level discoveries 
message('Mapping snp-level to gene-level discoveries.')
gene2gene = data.frame(geneCond = conditional$NEAREST_GENE,
                       geneCond_proteinCoding = as.numeric(conditional$NEAREST_GENE_BIOTYPE == 'protein_coding'),
                       geneCond_in_fastbat = as.numeric(conditional$NEAREST_GENE %in% gsub('.XY', '', fastbat.raw$Gene)),
                       geneCond_in_fastbatFDR = as.numeric(conditional$NEAREST_GENE %in% gsub('.XY', '', fastbat$Gene)),
                       geneCond_in_fastbatBonf = as.numeric(conditional$NEAREST_GENE %in% gsub('.XY', '', fastbat$Gene[fastbat$Bonf == 1])))
gene2gene.summary = data.frame(geneCond = nrow(conditional),
                               geneCond_proteinCoding = sum(conditional$NEAREST_GENE_BIOTYPE == 'protein_coding'),
                               geneCond_in_fastbat = sum(conditional$NEAREST_GENE %in%  gsub('.XY', '', fastbat.raw$Gene)),
                               geneCond_in_fastbatFDR = sum(conditional$NEAREST_GENE %in% gsub('.XY', '', fastbat$Gene)),
                               geneCond_in_fastbatBonf = sum(conditional$NEAREST_GENE %in% gsub('.XY', '', fastbat$Gene[fastbat$Bonf == 1])))

# (b) map snp-level to gene-level discoveries
fastbat.fdr = fastbat.raw[fastbat.raw$FDR < 0.05,]
fastbat.fdr = fastbat.fdr[order(fastbat.fdr$Pvalue),]
fastbat.bonf = fastbat.raw[fastbat.raw$Pvalue * nrow(fastbat.raw) < 0.05,]
fastbat.bonf = fastbat.bonf[order(fastbat.bonf$Pvalue),]

mapping = data.frame(matrix(data = NA, nrow = nrow(conditional), ncol = 6))
names(mapping) = c('snpCond', 'geneCond', 'nearestGeneFDR', 'nearestGeneFDR.distance', 'nearestGeneBonf', 'nearestGeneBonf.distance')

for (i in 1:nrow(conditional)) {
  mapping$snpCond[i] = conditional$ID[i]
  mapping$geneCond[i] = conditional$NEAREST_GENE[i]
  
  # get distance to FDR-significant genes
  tmp = fastbat.fdr[fastbat.fdr$Chr == conditional$CHR[i],]
  if (nrow(tmp) > 0) {
    for (j in 1:nrow(tmp)) {
      if (conditional$BP[i] > tmp$Start[j] & conditional$BP[i] < tmp$End[j]) {
        mapping$nearestGeneFDR[i] = tmp$Gene[j]
        mapping$nearestGeneFDR.distance[i] = 0
        next
      } else {
        distance = min(abs(conditional$BP[i] - tmp$Start[j]), abs(conditional$BP[i] - tmp$End[j]))
        if (distance < windowSize & (distance < mapping$nearestGeneFDR.distance[i] | is.na(mapping$nearestGeneFDR.distance[i]))) {
          mapping$nearestGeneFDR.distance[i] = distance
          mapping$nearestGeneFDR[i] = tmp$Gene[j]
        }
      }
    }
  }
  
  # get distance to Bonferroni-significant genes
  tmp = fastbat.bonf[fastbat.bonf$Chr == conditional$CHR[i],]
  if (nrow(tmp) > 0) {
    for (j in 1:nrow(tmp)) {
      if (conditional$BP[i] > tmp$Start[j] & conditional$BP[i] < tmp$End[j]) {
        mapping$nearestGeneBonf[i] = tmp$Gene[j]
        mapping$nearestGeneBonf.distance[i] = 0
        next
      } else {
        distance = min(abs(conditional$BP[i] - tmp$Start[j]), abs(conditional$BP[i] - tmp$End[j]))
        if (distance < windowSize & (distance < mapping$nearestGeneBonf.distance[i] | is.na(mapping$nearestGeneBonf.distance[i]))) {
          mapping$nearestGeneBonf[i] = tmp$Gene[j]
          mapping$nearestGeneBonf.distance[i] = distance
        }
      }
    }
  }
}

# get summary
mapping.summary = data.frame(snpCond = nrow(mapping),
                              geneCond = nrow(mapping),
                              nearestGeneFDR = sum(!is.na(mapping$nearestGeneFDR)),
                              nearestGeneFDR.maxDistance = max(mapping$nearestGeneFDR.distance, na.rm = T),
                              nearestGeneBonf = sum(!is.na(mapping$nearestGeneBonf)),
                              nearestGeneBonf.maxDistance = max(mapping$nearestGeneBonf.distance, na.rm = T))

# merge (a) and (b)
snp2gene = data.frame(snpCond = mapping$snpCond,
                      geneCond = mapping$geneCond,
                      geneCond.proteinCoding = gene2gene$geneCond_proteinCoding,
                      geneCond.inGeneBased = gene2gene$geneCond_in_fastbat,
                      geneCond.inGeneBasedFDR = gene2gene$geneCond_in_fastbatFDR,
                      geneCond.inGeneBasedBonf = gene2gene$geneCond_in_fastbatBonf,
                      nearestGeneFDR = mapping$nearestGeneFDR,
                      nearestGeneFDR.distance = mapping$nearestGeneFDR.distance,
                      nearestGeneBonf = mapping$nearestGeneBonf,
                      nearestGeneBonf.distance = mapping$nearestGeneBonf.distance)

snp2gene.summary = data.frame(snpCond = mapping.summary$snpCond,
                              geneCond = mapping.summary$geneCond,
                              geneCond.proteinCoding = gene2gene.summary$geneCond_proteinCoding,
                              geneCond.inGeneBased = gene2gene.summary$geneCond_in_fastbat,
                              geneCond.inGeneBasedFDR = gene2gene.summary$geneCond_in_fastbatFDR,
                              geneCond.inGeneBasedBonf = gene2gene.summary$geneCond_in_fastbatBonf,
                              nearestGeneFDR = mapping.summary$nearestGeneFDR,
                              nearestGeneFDR.maxDistance = mapping.summary$nearestGeneFDR.maxDistance,
                              nearestGeneBonf = mapping.summary$nearestGeneBonf,
                              nearestGeneBonf.maxDistance = mapping.summary$nearestGeneBonf.maxDistance)
  
# write files
message('Writing output files.')
write.table(snp2gene, file = paste0(outputPrefix, '.snp2gene'), quote = FALSE, row.names = FALSE, sep = '\t')
write.table(snp2gene.summary, file=paste0(outputPrefix, '.snp2gene.summary'), quote = FALSE, row.names = FALSE, sep = '\t')
write.table(gene2snp, file = paste0(outputPrefix, '.gene2snp'), quote = FALSE, row.names = FALSE, sep = '\t')
write.table(gene2snp.summary, file=paste0(outputPrefix, '.gene2snp.summary'), quote = FALSE, row.names = FALSE, sep = '\t')
message(paste0('--- Mapping gene-level to snp-level results completed. ---\n'))


