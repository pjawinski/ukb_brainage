#!/usr/bin/env Rscript

# ====================================================
# === Create table of SMR results for supplementum ===
# ====================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop(paste0('expected 5 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
snplevelFile=args[1] # snplevelFile="results/combined/snplevel.txt"
smrFiles=args[2] # smrFiles="results/gap_gm/smr/smr.eqtl.filtered.assigned.txt,results/gap_wm/smr/smr.eqtl.filtered.assigned.txt,results/gap_gwm/smr/smr.eqtl.filtered.assigned.txt"
outputFile=args[3] # outputFile="results/combined/smr.suppl.txt"
traits=args[4] # traits="gap_gm,gap_wm,gap_gwm"
traitNames=args[5] # traitNames="grey matter,white matter,grey and white matter"

message(paste0('\n--- Create table of SMR results for supplementum ---',
               '\nsnplevelFile: ', snplevelFile,
               '\nsmrFiles: ', smrFiles,
               '\noutputFile: ', outputFile,
               '\ntraits: ', traits,
               '\ntraitNames: ', traitNames, '\n'))

# attach packages to current R session
for (pkg in c('dplyr', 'stringr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# transform variables
traits = str_split(traits, ',')[[1]]
traitNames = str_split(traitNames, ',')[[1]]
smrFiles = str_split(smrFiles, ',')[[1]]

# load data
message('Loading Files.')
snplevel = read.delim(snplevelFile, sep = '\t', header = T)
snplevel$key = paste0(snplevel$TRAIT,"_", snplevel$ID)
snplevel = snplevel[,c('LOCUS_COUNT','DISCOV_COUNT','CHR', 'BP','key')]

# for i 
for (i in 1:length(smrFiles)) {
  tmp = read.table(smrFiles[i], sep = '\t', head = T)
  tmp$trait = traitNames[i]
  tmp$key = paste0(traits[i],"_", tmp$COND_LEAD_SNP)
  if (i == 1) { smr = tmp } else { smr = rbind(smr, tmp)}
}

# combine
df = inner_join(smr, snplevel, by = 'key')
df = df[order(df$LOCUS_COUNT, df$DISCOV_COUNT, df$p_SMR),]

# create supplementum data.frame
suppl = df[,c('LOCUS_COUNT','trait','COND_LEAD_SNP','CHR','BP',
              'probeID','Probe_bp','Gene','topSNP','topSNP_bp', 'COND_LEAD_SNP_RSQ', 'A1','A2','Freq',
              'b_GWAS','se_GWAS','p_GWAS','b_eQTL','se_eQTL','p_eQTL','b_SMR','se_SMR','p_SMR','FDR','p_HEIDI','nsnp_HEIDI')]

# write results
message('Writing output file.')
write.table(suppl, file = outputFile, row.names = F, quote = F, sep = '\t', na = "NA")
system(paste0('chmod 770 ', outputFile))
message('-- Creating suppl. table of SMR results completed. ---')
