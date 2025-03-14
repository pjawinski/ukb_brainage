#!/usr/bin/env Rscript

# ============================================================================
# === Calculate credible set overlap between SbayesRC, SusieR, and FINEMAP ===
# ============================================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=6) {
  stop(paste0('expected 6 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
sbayesfile=args[1] # sbayesfile="results/combined/gwama.eur.credible.sbayesrc.suppl.txt"
susieRfile=args[2] # susieRfile="results/combined/gwama.eur.credible.susieR.suppl.txt"
finemapfile=args[3] # finemapfile="results/combined/gwama.eur.credible.finemap.suppl.txt"
topsetsize=as.numeric(args[4]) # topsetsize=20
topsetoverlap=as.numeric(args[5]) # topsetoverlap=0.75
outputfile=args[6] # outputfile="results/combined/gwama.eur.credible.overlap"

logInfo = paste0('\n--- Calculating credible set overlap between SbayesRC, SusieR, and FINEMAP ---',
               '\nsbayesfile: ', sbayesfile,
               '\nsusieRfile: ', susieRfile,
               '\nfinemapfile: ', finemapfile,
               '\ntopsetsize: ', topsetsize,
               '\ntopsetoverlap: ', topsetoverlap,
               '\noutputfile: ', outputfile,'\n')
message(logInfo)

# attach packages to current R session
for (pkg in c('dplyr', 'stringr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# load data
message('Loading Files.')
sbayes = read.delim(sbayesfile, sep = '\t', header = T)
susieR = read.delim(susieRfile, sep = '\t', header = T)
finemap = read.delim(finemapfile, sep = '\t', header = T)
sbayes$key = paste0(sbayes$LOCUS_COUNT, ":", sbayes$trait)
susieR$key = paste0(susieR$LOCUS_COUNT, ":", susieR$trait)
finemap$key = paste0(finemap$LOCUS_COUNT, ":", finemap$trait)

# calculate overlap between susieR and finemap
message('Calculating overlap between SusieR and FINEMAP')
allkeys = intersect(unique(susieR$key),unique(finemap$key))
results = data.frame(matrix(nrow = length(allkeys), ncol = 9))
names(results) = c('locus','trait','susieR_count','finemap_count','total_count','overlap_count','prop_overlap','prop_susieR_in_finemap','prop_finemap_in_susieR')
k = 0
for (key in allkeys) {

  # subset data for the current locus:trait
  k = k + 1
  susieR.tmp = susieR[susieR$key == key,]
  finemap.tmp = finemap[finemap$key == key,]

  # loop over credible sets
  results.tmp = data.frame(matrix(nrow = 1, ncol = 9))
  names(results.tmp) = c('cs_susieR','cs_finemap','susieR_count','finemap_count','total_count','overlap_count','prop_overlap','prop_susieR_in_finemap','prop_finemap_in_susieR')
  l = 0
  for (cssusieR in unique(susieR.tmp$cs)) {
    for (csfinemap in unique(finemap.tmp$signalCount)) {
      l = l + 1
      susieRvariants = susieR.tmp$rsid[susieR.tmp$cs==cssusieR]
      finemapvariants = finemap.tmp$rsid[finemap.tmp$signalCount==csfinemap]
      shared_variants = intersect(susieRvariants,finemapvariants)
      total_variants = union(susieRvariants,finemapvariants)

      # Metrics
      overlap_count = length(shared_variants)
      jaccard_index = length(shared_variants) / length(total_variants)
      proportion_susieR_in_finemap = length(shared_variants) / length(susieRvariants)
      proportion_finemap_in_susieR = length(shared_variants) / length(finemapvariants)
     
      # save temporary results
      results.tmp[l,] = c(cssusieR,csfinemap,length(susieRvariants),length(finemapvariants),length(total_variants),overlap_count,jaccard_index,proportion_susieR_in_finemap,proportion_finemap_in_susieR)
    }
  }

  # keep best combination of credible variants
  results.tmp = results.tmp[order(results.tmp$cs_susieR,-results.tmp$prop_overlap),]
  results.tmp = results.tmp[!duplicated(results.tmp$cs_susieR),]
  results[k,] = c(susieR.tmp$LOCUS_COUNT[1],
    susieR.tmp$trait[1],
    paste0(results.tmp$susieR_count,collapse = " | "),
    paste0(results.tmp$finemap_count,collapse = " | "),
    paste0(results.tmp$total_count,collapse = " | "),
    paste0(results.tmp$overlap_count,collapse = " | "),
    paste0(results.tmp$prop_overlap,collapse = " | "),
    paste0(results.tmp$prop_susieR_in_finemap,collapse = " | "),
    paste0(results.tmp$prop_finemap_in_susieR,collapse = " | "))
}
values = strsplit(results$prop_overlap, " \\| ")
values = as.numeric(unlist(values))
median_values = median(values, na.rm = TRUE)
results_susieR_finemap = results
tmpInfo = sprintf('Median credible set overlap between SusieR and FINEMAP: %0.3f',median_values)
message(tmpInfo)
logInfo = paste(logInfo, tmpInfo, sep = '\n')


# calculate overlap between sbayes and susieR
allkeys = intersect(unique(sbayes$key),unique(susieR$key))
results = data.frame(matrix(nrow = length(allkeys), ncol = 9))
names(results) = c('locus','trait','sbayes_count','susieR_count','total_count','overlap_count','prop_overlap','prop_sbayes_in_susieR','prop_susieR_in_sbayes')
k = 0
for (key in allkeys) {

  # subset data for the current locus:trait
  k = k + 1
  sbayes.tmp = sbayes[sbayes$key == key,]
  susieR.tmp = susieR[susieR$key == key,]

  # loop over credible sets
  results.tmp = data.frame(matrix(nrow = 1, ncol = 9))
  names(results.tmp) = c('cs_sbayes','cs_susieR','sbayes_count','susieR_count','total_count','overlap_count','prop_overlap','prop_sbayes_in_susieR','prop_susieR_in_sbayes')
  l = 0
  for (cssbayes in unique(sbayes.tmp$csCount)) {
    for (cssusieR in unique(susieR.tmp$cs)) {
      l = l + 1
      sbayesvariants = sbayes.tmp$Name[sbayes.tmp$csCount==cssbayes]
      susieRvariants = susieR.tmp$rsid[susieR.tmp$cs==cssusieR]
      shared_variants = intersect(sbayesvariants,susieRvariants)
      total_variants = union(sbayesvariants,susieRvariants)

      # Metrics
      overlap_count = length(shared_variants)
      jaccard_index = length(shared_variants) / length(total_variants)
      proportion_sbayes_in_susieR = length(shared_variants) / length(sbayesvariants)
      proportion_susieR_in_sbayes = length(shared_variants) / length(susieRvariants)
     
      # save temporary results
      results.tmp[l,] = c(cssbayes,cssusieR,length(sbayesvariants),length(susieRvariants),length(total_variants),overlap_count,jaccard_index,proportion_sbayes_in_susieR,proportion_susieR_in_sbayes)
    }
  }

  # keep best combination of credible variants
  results.tmp = results.tmp[order(results.tmp$cs_sbayes,-results.tmp$prop_overlap),]
  results.tmp = results.tmp[!duplicated(results.tmp$cs_sbayes),]
  results[k,] = c(sbayes.tmp$LOCUS_COUNT[1],
    sbayes.tmp$trait[1],
    paste0(results.tmp$sbayes_count,collapse = " | "),
    paste0(results.tmp$susieR_count,collapse = " | "),
    paste0(results.tmp$total_count,collapse = " | "),
    paste0(results.tmp$overlap_count,collapse = " | "),
    paste0(results.tmp$prop_overlap,collapse = " | "),
    paste0(results.tmp$prop_sbayes_in_susieR,collapse = " | "),
    paste0(results.tmp$prop_susieR_in_sbayes,collapse = " | "))
}
values = strsplit(results$prop_sbayes_in_susieR, " \\| ")
values = as.numeric(unlist(values))
median_values = median(values, na.rm = TRUE)
results_sbayes_susieR = results
tmpInfo = sprintf('Median proportion of SBayesRC credible variants in SusieR credible sets: %0.3f',median_values)
message(tmpInfo)
logInfo = paste(logInfo, tmpInfo, sep = '\n')

# calculate overlap between sbayes and finemap
allkeys = intersect(unique(sbayes$key),unique(finemap$key))
results = data.frame(matrix(nrow = length(allkeys), ncol = 9))
names(results) = c('locus','trait','sbayes_count','finemap_count','total_count','overlap_count','prop_overlap','prop_sbayes_in_finemap','prop_finemap_in_sbayes')
k = 0
for (key in allkeys) {

  # subset data for the current locus:trait
  k = k + 1
  sbayes.tmp = sbayes[sbayes$key == key,]
  finemap.tmp = finemap[finemap$key == key,]

  # loop over credible sets
  results.tmp = data.frame(matrix(nrow = 1, ncol = 9))
  names(results.tmp) = c('cs_sbayes','cs_finemap','sbayes_count','finemap_count','total_count','overlap_count','prop_overlap','prop_sbayes_in_finemap','prop_finemap_in_sbayes')
  l = 0
  for (cssbayes in unique(sbayes.tmp$csCount)) {
    for (csfinemap in unique(finemap.tmp$signalCount)) {
      l = l + 1
      sbayesvariants = sbayes.tmp$Name[sbayes.tmp$csCount==cssbayes]
      finemapvariants = finemap.tmp$rsid[finemap.tmp$signalCount==csfinemap]
      shared_variants = intersect(sbayesvariants,finemapvariants)
      total_variants = union(sbayesvariants,finemapvariants)

      # Metrics
      overlap_count = length(shared_variants)
      jaccard_index = length(shared_variants) / length(total_variants)
      proportion_sbayes_in_finemap = length(shared_variants) / length(sbayesvariants)
      proportion_finemap_in_sbayes = length(shared_variants) / length(finemapvariants)
     
      # save temporary results
      results.tmp[l,] = c(cssbayes,csfinemap,length(sbayesvariants),length(finemapvariants),length(total_variants),overlap_count,jaccard_index,proportion_sbayes_in_finemap,proportion_finemap_in_sbayes)
    }
  }

  # keep best combination of credible variants
  results.tmp = results.tmp[order(results.tmp$cs_sbayes,-results.tmp$prop_overlap),]
  results.tmp = results.tmp[!duplicated(results.tmp$cs_sbayes),]
  results[k,] = c(sbayes.tmp$LOCUS_COUNT[1],
    sbayes.tmp$trait[1],
    paste0(results.tmp$sbayes_count,collapse = " | "),
    paste0(results.tmp$finemap_count,collapse = " | "),
    paste0(results.tmp$total_count,collapse = " | "),
    paste0(results.tmp$overlap_count,collapse = " | "),
    paste0(results.tmp$prop_overlap,collapse = " | "),
    paste0(results.tmp$prop_sbayes_in_finemap,collapse = " | "),
    paste0(results.tmp$prop_finemap_in_sbayes,collapse = " | "))
}
values = strsplit(results$prop_sbayes_in_finemap, " \\| ")
values = as.numeric(unlist(values))
median_values = median(values, na.rm = TRUE)
results_sbayes_finemap = results
tmpInfo = sprintf('Median proportion of SBayesRC credible variants and FINEMAP credible sets: %0.3f',median_values)
message(tmpInfo)
logInfo = paste(logInfo, tmpInfo, sep = '\n')

# get top credible sets (small size, high consistency actoss fine-mapping strategies)
check_size <- function(x) {
  values <- as.numeric(unlist(strsplit(x, " \\| ")))
  any(values <= topsetsize, na.rm = TRUE)
}
check_overlap <- function(x) {
  values <- as.numeric(unlist(strsplit(x, " \\| ")))
  any(values >= topsetoverlap, na.rm = TRUE)
}
topsets.part1 = results_susieR_finemap[sapply(results_susieR_finemap$susieR_count,check_size) & sapply(results_susieR_finemap$susieR_count,check_size),c('locus','trait','susieR_count','finemap_count','prop_overlap')]
names(topsets.part1) = c('locus','trait','susieR_count','finemap_count','susieR_finemap_overlap')
topsets.part2 = results_sbayes_susieR[sapply(results_sbayes_susieR$sbayes_count,check_size) & sapply(results_sbayes_susieR$susieR_count,check_size),c('locus','trait','sbayes_count','prop_overlap','prop_sbayes_in_susieR')]
names(topsets.part2) = c('locus','trait','sbayes_count','sbayes_susieR_overlap','prop_sbayes_in_susieR')
topsets.part3 = results_sbayes_finemap[sapply(results_sbayes_finemap$sbayes_count,check_size) & sapply(results_sbayes_finemap$finemap_count,check_size),c('locus','trait','prop_overlap','prop_sbayes_in_finemap')]
names(topsets.part3) = c('locus','trait','sbayes_finemap_overlap','prop_sbayes_in_finemap')
topsets = inner_join(topsets.part1, topsets.part2, by = c('locus','trait'))
topsets = inner_join(topsets, topsets.part3, by = c('locus','trait'))
topsets = topsets[sapply(topsets$susieR_finemap_overlap,check_overlap) & sapply(topsets$prop_sbayes_in_susieR,check_overlap) & sapply(topsets$prop_sbayes_in_finemap,check_overlap) ,]
topsets = topsets[,c('locus','trait','sbayes_count','susieR_count','finemap_count','susieR_finemap_overlap','prop_sbayes_in_susieR','prop_sbayes_in_finemap')]

# write results
message(sprintf('\nWriting %s.susieR.finemap.txt',outputfile))
  write.table(results_susieR_finemap, file = sprintf('%s.susieR.finemap.txt',outputfile), row.names = F, quote = F, sep = '\t', na = "NA")
  system(sprintf('chmod 770 %s.susieR.finemap.txt', outputfile))
message(sprintf('Writing %s.sbayes.susieR.txt',outputfile))
  write.table(results_sbayes_susieR, file = sprintf('%s.sbayes.susieR.txt',outputfile), row.names = F, quote = F, sep = '\t', na = "NA")
  system(sprintf('chmod 770 %s.sbayes.susieR.txt', outputfile))
message(sprintf('Writing %s.sbayes.finemap.txt',outputfile))
  write.table(results_sbayes_finemap, file = sprintf('%s.sbayes.finemap.txt',outputfile), row.names = F, quote = F, sep = '\t', na = "NA")
  system(sprintf('chmod 770 %s.sbayes.finemap.txt', outputfile))
message(sprintf('Writing %s.topsets.txt',outputfile))
  write.table(topsets, file = sprintf('%s.topsets.txt',outputfile), row.names = F, quote = F, sep = '\t', na = "NA")
  system(sprintf('chmod 770 %s.topsets.txt', outputfile))

# save log file
message(sprintf('Writing %s.log',outputfile))
sink(sprintf('%s.log',outputfile)) 
sink(stdout(), type = "message")
message(logInfo)
sink()
system(sprintf('chmod -R 770  %s.log', outputfile))
message('-- Completed: Calculating credible set overlap between SbayesRC, SusieR, and FINEMA ---')
