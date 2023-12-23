#!/usr/bin/env Rscript

# ======================================================
# === get replication results of discovered variants ===
# ======================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop(paste0('expected 4 arguments, but ', length(args), ' arguments provided.'), call.=FALSE)
}

# get arguments from command line
discoveryFile = args[1] # discoveryFile="results/gap_gm/conditional/conditional.cleaned.tophits.annovar.txt"
replicationFile = args[2] # replicationFile="results/gap_gm/replicate/metal.ivweight.qc.gz"
replicationNcovs = as.numeric(args[3]) # replicationNcovs = 8 # sex,age,age2,TIV,PC1-4 (8-18)
outFile = args[4] # outFile="results/gap_gm/replicate/replicateDiscoveries"
logInfo = paste0('\n--- Get replication results of discovered variants ---',
               '\ndiscoveryFile: ', discoveryFile,
               '\nreplicationFile: ', replicationFile,
               '\nreplicationNcovs: ', replicationNcovs,
               '\noutFile: ', outFile,'\n')
message(logInfo)

# load required pages
for (pkg in c('dplyr','data.table')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# load data
message(' - loading data.')
discov = read.delim(discoveryFile, sep = '\t', head = TRUE)
replic = data.frame(fread(cmd=sprintf('gzip -dc %s', replicationFile), tmpdir = getwd(), header=T, stringsAsFactors=FALSE))

# change variable names
idx = names(discov)!='ID'
names(discov)[idx] = paste0('DISCOV_',names(discov)[idx])
idx = names(replic)!='ID'
names(replic)[idx] = paste0('REPLIC_',names(replic)[idx])

# merge data frames and prioritize replication variants
# priority list: order results by DISCOV_LEAD_SNP_P, REPLIC_N, DISCOV_LEAD_SNP_RSQ, DISCOV_P)
message(' - merging discovery and replication dataset.')
df = inner_join(discov,replic,'ID')
df = df[order(df$DISCOV_LEAD_SNP_P, -df$REPLIC_N, -df$DISCOV_LEAD_SNP_RSQ, df$DISCOV_P),]
df = df[!duplicated(df$DISCOV_LOCUS_CNT),]

# test if alleles match
if (nrow(df) == sum((df$DISCOV_A1 == df$REPLIC_A1 | df$DISCOV_A1 == df$REPLIC_A2) & (df$DISCOV_A2 == df$REPLIC_A1 | df$DISCOV_A2 == df$REPLIC_A2))) {
	message(' - discovery and replication alleles are the same')
}

# invert results if A1 == A2
idx = df$REPLIC_A1 == df$DISCOV_A2
A1 = df$REPLIC_A1[idx]
A2 = df$REPLIC_A2[idx]
df$REPLIC_A1[idx] = A2
df$REPLIC_A2[idx] = A1
df$REPLIC_A1_FREQ[idx] = 1-df$REPLIC_A1_FREQ[idx]
A1_FREQ_MAX = 1-df$REPLIC_A1_FREQ_MIN[idx]
A1_FREQ_MIN = 1-df$REPLIC_A1_FREQ_MAX[idx]
df$REPLIC_A1_FREQ_MAX[idx] = A1_FREQ_MAX
df$REPLIC_A1_FREQ_MIN[idx] = A1_FREQ_MIN
df$REPLIC_BETA[idx] = -df$REPLIC_BETA[idx]
df$REPLIC_DIRECTION[idx] = df$REPLIC_DIRECTION[idx] %>% gsub('-','!', .) %>% gsub('[+]','-', .) %>% gsub('!','+', .) 

# add ISLEAD, Z and ETA
df$DISCOV_ISLEAD = FALSE
df$DISCOV_ISLEAD[df$DISCOV_LEAD_SNP == df$ID] = TRUE
df$REPLIC_Z = df$REPLIC_BETA/df$REPLIC_SE
df$REPLIC_ETA2 = df$REPLIC_Z^2/(df$REPLIC_Z^2+(df$REPLIC_N-2-replicationNcovs))

# make columns for replication consistency
df$CONSIST_DIRECTION = sign(df$DISCOV_Z) == sign(df$REPLIC_Z)
df$CONSIST_REPLICATED = df$CONSIST_DIRECTION & df$REPLIC_P/2 < 0.05

# create table
message(sprintf(' - saving  %s.txt.',outFile))
df$`NA` = ""
cols = c('DISCOV_CHR','DISCOV_BP','ID','DISCOV_A1','DISCOV_A2',
		 'DISCOV_CYTOBAND','DISCOV_REGION','DISCOV_NEAREST_GENE','DISCOV_NEAREST_GENE_DESCRIPTION','DISCOV_NEAREST_GENE_BIOTYPE','DISCOV_DISTANCE',
		'NA','DISCOV_A1_FREQ','DISCOV_BETA','DISCOV_SE','DISCOV_Z','DISCOV_P','DISCOV_ETA2','DISCOV_N',
		'NA','REPLIC_A1_FREQ','REPLIC_BETA','REPLIC_SE','REPLIC_Z','REPLIC_P','REPLIC_ETA2','REPLIC_N',
		'REPLIC_DIRECTION','REPLIC_HET_ISQ','REPLIC_HetChiSq','REPLIC_HET_DF','REPLIC_HET_P',
		'CONSIST_DIRECTION','CONSIST_REPLICATED',
		'NA','DISCOV_ISLEAD','DISCOV_LEAD_SNP','DISCOV_LEAD_SNP_P','DISCOV_LEAD_SNP_pJ','DISCOV_LEAD_SNP_RSQ','DISCOV_LEAD_SNP_ALLELES')
output = df[,cols]
write.table(output, file = sprintf('%s.txt',outFile), sep = '\t', col.names = T, row.names = F, quote = F)
system(sprintf('chmod -R 770 %s.txt', outFile))

# save log file
message(sprintf(' - saving  %s.log.',outFile))
sink(sprintf('%s.log',outFile))
sink(stdout(), type = "message")
message(logInfo)
sink()
system(sprintf('chmod -R 770 %s.log', outFile))
message('--- Completed: Get replication results of discovered variants ---')

