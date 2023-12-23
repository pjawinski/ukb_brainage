#!/usr/bin/env Rscript

# =====================================================================
# === Create table of replication results for snp-level discoveries ===
# =====================================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop(paste0('expected 5 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
snplevelFile=args[1] # snplevelFile="results/combined/snplevel.txt"
replicateFiles=args[2] # replicateFiles="results/gap_gm/replicate/replicateDiscoveries.txt,results/gap_wm/replicate/replicateDiscoveries.txt,results/gap_gwm/replicate/replicateDiscoveries.txt"
outFile=args[3] # outFile="results/combined/replicateDiscoveries"
traits=args[4] # traits="gap_gm,gap_wm,gap_gwm"
traitNames=args[5] # traitNames="grey matter,white matter,grey and white matter"

logInfo = paste0('\n--- Creating table of replication results for snp-level discoveries ---',
               '\nsnplevelFile: ', snplevelFile,
               '\nreplicateFiles: ', replicateFiles,
               '\noutFile: ', outFile,
               '\ntraits: ', traits,
               '\ntraitNames: ', traitNames, '\n')
message(logInfo)

# attach packages to current R session
for (pkg in c('dplyr', 'stringr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# transform variables
traits = str_split(traits, ',')[[1]]
traitNames = str_split(traitNames, ',')[[1]]
traitNames = data.frame(TRAIT = traits, TRAITNAME = traitNames)
replicateFiles = str_split(replicateFiles, ',')[[1]]

# load data
status = 'Loading Files.'; message(status); logInfo = sprintf('%s\n%s',logInfo,status)
df = read.delim(snplevelFile, sep = '\t', header = T)
df = left_join(df,traitNames,by = 'TRAIT')
df$key = sprintf('%s:%s',df$TRAIT,df$ID)
df = df[,c('key','LOCUS_COUNT','DISCOV_COUNT','TRAIT','TRAITNAME','CYTOBAND','LEAD_SNP_GWS','ID','CHR','BP','A1','A2','A1_FREQ','BETA','SE','Z','P','COJO_pJ','ETA2','N','NEAREST_GENE','DISTANCE')]
names(df)[2:length(names(df))] = sprintf('DISCOV_%s',names(df)[2:length(names(df))])

# load replicateFiles
for (i in 1:length(replicateFiles)) {
  tmp = read.delim(replicateFiles[i])
  key = sprintf('%s:%s',traits[i],tmp$DISCOV_LEAD_SNP)
  tmp = tmp[,-which(names(tmp) %in% c('DISCOV_A1_FREQ','DISCOV_BETA','DISCOV_SE','DISCOV_Z','DISCOV_P','DISCOV_ETA2','DISCOV_N','DISCOV_LEAD_SNP','DISCOV_LEAD_SNP_P','DISCOV_LEAD_SNP_pJ'))]
  cols2change=c('DISCOV_CHR','DISCOV_BP','DISCOV_A1','DISCOV_A2','DISCOV_CYTOBAND','DISCOV_REGION','DISCOV_NEAREST_GENE','DISCOV_NEAREST_GENE_DESCRIPTION','DISCOV_NEAREST_GENE_BIOTYPE','DISCOV_DISTANCE','DISCOV_ISLEAD')
  names(tmp) = str_replace_all(names(tmp), 'DISCOV_', '')
  names(tmp) = str_replace_all(names(tmp), 'REPLIC_', '')
  names(tmp) = sprintf('REPLIC_%s',names(tmp))
  tmp$key = key
  if (i==1) { repl = tmp } else { repl = rbind(repl, tmp) }
}

# merge files
status = 'Merging discovery and replication tables.'; message(status); logInfo = sprintf('%s\n%s',logInfo,status)
df = left_join(df,repl,by = 'key')

# add column that indicates whether association refers to a genome-wide significant discovery
df$DISCOV_GWS = df$DISCOV_COJO_pJ < 5E-8

# get some stats
status = 'Calculating statistics.'; message(status); logInfo = sprintf('%s\n%s',logInfo,status)
statistics = paste0('\nTRAIT\tpUpper\tpLower\tconsistentDirection\tconsistentDirectionP\tnominalSig\tnominalSigP\tnominalSigGenes')

  # trait-wise
  pUpper = c(1,1E-6,5E-8)
  pLower = c(0,5E-8,0)
  for (i in 1:length(traits)) {
    for (j in 1:length(pLower)) {
      tmp = df[df$DISCOV_TRAIT==traits[i] & df$DISCOV_COJO_pJ >= pLower[j] & df$DISCOV_COJO_pJ < pUpper[j],]
      ntotal = nrow(tmp)
      directionCount = sum(tmp$REPLIC_CONSIST_DIRECTION)
      nominalSigCount = sum(tmp$REPLIC_CONSIST_REPLICATED)
      directionP = pbinom(directionCount-1, ntotal, 0.5, lower.tail = FALSE)
      nominalSigP = pbinom(nominalSigCount-1, ntotal, 0.05, lower.tail = FALSE)
      location = rep('in',nrow(tmp))
      location[tmp$DISCOV_DISTANCE > 0] = 'near'
      idx = tmp$REPLIC_CONSIST_REPLICATED
      nominalSigGenes = paste0(sprintf('%s (%s %s %s)', tmp$DISCOV_CYTOBAND[idx], tmp$DISCOV_ID[idx], location[idx], tmp$DISCOV_NEAREST_GENE[idx]), collapse = ', ')
      statistics = sprintf('%s\n%s\t%.1e\t%.1e\t%d/%d\t%.1e\t%d/%d\t%.1e\t%s',statistics,traits[i],pUpper[j],pLower[j],directionCount,ntotal,directionP,nominalSigCount,ntotal,nominalSigP,nominalSigGenes)
    }
  }

  # across all traits
  for (j in 1:length(pLower)) {
    tmp = df[df$DISCOV_DISCOV_COUNT == 1 & df$DISCOV_COJO_pJ >= pLower[j] & df$DISCOV_COJO_pJ < pUpper[j],]
    ntotal = nrow(tmp)
    directionCount = sum(tmp$REPLIC_CONSIST_DIRECTION)
    nominalSigCount = sum(tmp$REPLIC_CONSIST_REPLICATED)
    directionP = pbinom(directionCount-1, ntotal, 0.5, lower.tail = FALSE)
    nominalSigP = pbinom(nominalSigCount-1, ntotal, 0.05, lower.tail = FALSE)
    location = rep('in',nrow(tmp))
    location[tmp$DISCOV_DISTANCE > 0] = 'near'
    idx = tmp$REPLIC_CONSIST_REPLICATED
    nominalSigGenes = paste0(sprintf('%s (%s %s %s)', tmp$DISCOV_CYTOBAND[idx], tmp$DISCOV_ID[idx], location[idx], tmp$DISCOV_NEAREST_GENE[idx]), collapse = ', ')
    statistics = sprintf('%s\n%s\t%.1e\t%.1e\t%d/%d\t%.1e\t%d/%d\t%.1e\t%s',statistics,'all',pUpper[j],pLower[j],directionCount,ntotal,directionP,nominalSigCount,ntotal,nominalSigP,nominalSigGenes)
  }
  message(statistics)
  logInfo = sprintf('%s\n%s',logInfo,statistics)

# remove R2 and phasing information if replication variant equals discovery variant
idx = df$DISCOV_ID == df$REPLIC_ID
df[idx,c('REPLIC_LEAD_SNP_RSQ','REPLIC_LEAD_SNP_ALLELES')] = '-'

# create output
df$`NA` = ""
cols=c('DISCOV_LOCUS_COUNT','DISCOV_TRAITNAME','DISCOV_CYTOBAND','DISCOV_GWS',
  'NA','DISCOV_ID','DISCOV_CHR','DISCOV_BP','DISCOV_A1','DISCOV_A2','DISCOV_A1_FREQ','DISCOV_BETA','DISCOV_SE','DISCOV_Z','DISCOV_P','DISCOV_COJO_pJ','DISCOV_ETA2','DISCOV_N',
  'NA','REPLIC_ID','REPLIC_CHR','REPLIC_BP','REPLIC_A1','REPLIC_A2','REPLIC_ISLEAD','REPLIC_LEAD_SNP_RSQ','REPLIC_LEAD_SNP_ALLELES','REPLIC_A1_FREQ','REPLIC_BETA','REPLIC_SE','REPLIC_Z','REPLIC_P','REPLIC_ETA2','REPLIC_N',
  'REPLIC_DIRECTION','REPLIC_HET_ISQ','REPLIC_HetChiSq','REPLIC_HET_DF','REPLIC_HET_P','REPLIC_CONSIST_DIRECTION','REPLIC_CONSIST_REPLICATED',
  'NA','REPLIC_REGION','REPLIC_NEAREST_GENE','REPLIC_NEAREST_GENE_DESCRIPTION','REPLIC_NEAREST_GENE_BIOTYPE','REPLIC_DISTANCE')
output = df[,cols]

# write results
status = sprintf('\nWriting %s.txt',outFile); message(status); logInfo = sprintf('%s\n%s',logInfo,status)
write.table(output, file = sprintf('%s.txt',outFile), row.names = F, quote = F, sep = '\t')
system(sprintf('chmod -R 770 %s.txt',outFile))

# save log file
status = sprintf('Writing %s.log',outFile); message(status); logInfo = sprintf('%s\n%s',logInfo,status)
sink(sprintf('%s.log',outFile)) 
sink(stdout(), type = "message")
message(logInfo)
sink()
system(sprintf('chmod -R 770 %s.log',outFile))
status = '-- Completed: Creating table of replication results for snp-level discoveries ---'; message(status); logInfo = sprintf('%s\n%s',logInfo,status)

