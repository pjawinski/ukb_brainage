#!/usr/bin/env Rscript

# ==================================================
# === Run gene set enrichment test using GOfuncR ===
# ==================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=6) {
  stop(paste0('expected 6 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
fastbatFile = args[1] # fastbatFile="results/combined/fastbat.txt"
fastbatpCol = args[2] # fastbatpCol="gap_gm_Pvalue"
fastbatGeneCol = args[3] # fastbatGeneCol="Gene"
fastbatClumpCol = args[4] # fastbatClumpCol="LOCUS_COUNT"
crossFWER = args[5] # crossFWER=TRUE
outFile = args[6] # outFile="results/gap_gm/gofuncr/gofuncr.gsea"

logInfo = paste0('\n--- Run gene set enrichment test using GOfuncR ---',
               '\nfastbatFile: ', fastbatFile,
               '\nfastbatpCol: ', fastbatpCol,
               '\nfastbatGeneCol: ', fastbatGeneCol,
               '\nfastbatClumpCol: ', fastbatClumpCol,
               '\ncrossFWER: ', crossFWER,
               '\noutFile: ', outFile,'\n')
message(logInfo)

# attach required packages
for (pkg in c('dplyr','stringr','GOfuncR')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# create input dataframe with candidate genes
message('Load input genes....')
fb = read.delim(fastbatFile, sep = '\t', header = T)
fb = fb[,c(fastbatGeneCol,fastbatpCol,fastbatClumpCol)]
fb = fb[order(fb[,fastbatpCol],fb[,fastbatGeneCol]),]
input_willi = data.frame(gene_id = fb[,fastbatGeneCol], gene_scores = fb[,fastbatpCol])

# run overrepresentation test
message('Run gene set enrichment test...')
set.seed(2340)
res_willi = go_enrich(input_willi, test='wilcoxon', n_randsets=10000, silent = T)
stats = res_willi[[1]]

# get size
message('Add GO size...')
GOgenes = get_anno_genes(stats$node_id,genes = res_willi[[2]]$gene_id)
refSize = data.frame(table(GOgenes$go_id))
names(refSize) = c('node_id','size')
gsea = left_join(stats,refSize,by = 'node_id')

# add mean rank and number of unique loci
message('Add mean rank...')
input_willi = cbind(input_willi, locus = fb[,fastbatClumpCol])
input_willi_mapped = input_willi[input_willi$gene_id %in% res_willi[[2]]$gene_id,]
input_willi_mapped$rank = rank(input_willi_mapped$gene_scores)
gsea$rank = gsea$locus_count = NA
for (i in 1:nrow(gsea)) {
  if (i%%100 == 0) { message(sprintf(' - %d out of %d completed.',i, nrow(gsea))) }
  idx = input_willi_mapped$gene_id %in% GOgenes$gene[GOgenes$go_id == gsea$node_id[i]]
  go = input_willi_mapped$rank[idx]
  gsea$rank[i] = mean(go)
  gsea$locus_count[i] = length(unique(input_willi_mapped$locus[idx]))
}
gsea$meanRank = mean(input_willi_mapped$rank)

# calculate FWER across the three ontologies
if (crossFWER == TRUE) {
  message('Calculate FWER across ontologies..')
  randomPvals = cbind(res_willi[[4]]$lower[res_willi[[4]]$ontology == 'biological_process'],
    res_willi[[4]]$lower[res_willi[[4]]$ontology == 'cellular_component'],
    res_willi[[4]]$lower[res_willi[[4]]$ontology == 'molecular_function']) %>%
    apply(1, FUN = min, na.rm = T)
  for (i in 1:nrow(gsea)) {
    pval = gsea$raw_p_low_rank[i]
    gsea$FWER_low_rank[i] = sum(randomPvals <= pval)/10000
  }
}
# perform refinement for categories with FWER < 0.05
message('Perform refinement using elim algorithm...')
if (crossFWER == TRUE) { res_willi[[1]]$FWER_low_rank = gsea$FWER_low_rank }
refined = refine(res_willi, fwer = 0.05, fwer_col = 6)

# add post-refinement FDR
message('Add post-refinement FDR...')
for (i in 1:nrow(refined)) {
  pval = refined$refined_p_low_rank[i]
  if (crossFWER != TRUE) {
    ontology = refined$ontology[i]
    randomPvals = res_willi[[4]][res_willi[[4]][,1] == ontology, 2]
  }
  refined$refined_FWER_low_rank[i] = sum(randomPvals <= pval)/10000
}
refined$signif = refined$refined_FWER_low_rank < 0.05
gsea = left_join(gsea,refined[,c('node_id','refined_p_low_rank','refined_FWER_low_rank','signif')],by = 'node_id')

# create log
tmpMessage = paste0(capture.output(res_willi[[3]]), collapse = '\n')
message(tmpMessage)
logInfo = c(logInfo,'\n',tmpMessage)

# write results
message(sprintf(' - writing %s.results.txt',outFile))
write.table(gsea, sprintf('%s.results.txt',outFile), sep = '\t', quote = F, row.names = F)
gsea = gsea[order(gsea$raw_p_low_rank),]

# write list of valid genes
message(sprintf(' - writing %s.mappedGenes.txt',outFile))
write.table(fb[fb[,fastbatGeneCol] %in% res_willi[[2]]$gene_id,], sprintf('%s.mappedGenes.txt',outFile), sep = '\t', quote = F, row.names = F)

# save log file
message(sprintf(' - writing %s.log', outFile)) 
sink(sprintf('%s.log', outFile))
sink(stdout(), type = "message")
message(logInfo)
sink()
system(sprintf('chmod -R 770 %s*', outFile))
message('--- Completed: Run gene set enrichment test using GOfuncR ---')
