#!/usr/bin/env Rscript

# ========================================================
# === Create table of FINEMAP results for supplementum ===
# ========================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop(paste0('expected 5 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
snplevelFile=args[1] # snplevelFile="results/combined/gwama.eur.snplevel/snplevel.gws.txt"
finemapFiles=args[2] # finemapFiles="results/gap_gm/gwama/eur/susieR/susieR.df.purity50.annot.txt,results/gap_wm/gwama/eur/susieR/susieR.df.purity50.annot.txt,results/gap_gwm/gwama/eur/susieR/susieR.df.purity50.annot.txt"
outputFile=args[3] # outputFile="results/combined/gwama.eur.credible.susieR.suppl.txt"
traits=args[4] # traits="gap_gm,gap_wm,gap_gwm"
traitNames=args[5] # traitNames="grey matter,white matter,grey and white matter"

message(paste0('\n--- Creating table of SusieR fine-mapping results for supplementum ---',
               '\nsnplevelFile: ', snplevelFile,
               '\nfinemapFiles: ', finemapFiles,
               '\noutputFile: ', outputFile,
               '\ntraits: ', traits,
               '\ntraitNames: ', traitNames, '\n'))

# attach packages to current R session
for (pkg in c('dplyr', 'stringr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# transform variables
traits = str_split(traits, ',')[[1]]
traitNames = str_split(traitNames, ',')[[1]]
finemapFiles = str_split(finemapFiles, ',')[[1]]

# load data
message('Loading Files.')
snplevel = read.delim(snplevelFile, sep = '\t', header = T)
snplevel$key = paste0(snplevel$TRAIT,"_", snplevel$ID)
snplevel = snplevel[,c('LOCUS_COUNT','DISCOV_COUNT','BP','key')]

# add traitNames and identifier
for (i in 1:length(finemapFiles)) {
  tmp = read.table(finemapFiles[i], sep = '\t', head = T, quote = "")
  tmp$trait = traitNames[i]
  tmp$key = paste0(traits[i],"_", tmp$index_rsid)
  if (i == 1) { finemap = tmp } else { finemap = rbind(finemap, tmp)}
}

# combine
message('Join credible sets with snplevel results.')
df = inner_join(finemap, snplevel, by = 'key')
df = df[order(df$LOCUS_COUNT, df$DISCOV_COUNT, df$cs, df$cumpip),]

# create supplementum data.frame
suppl = df[,c('LOCUS_COUNT','trait','index_rsid','chromosome','BP','minBP','maxBP','rangeBP','snpTotal','snpCausal',
              'cs','cs.size','coverage','purity.min','purity.median',
              'rsid','position','allele1','allele2','maf','beta','se','pip','cumpip','r2',
              'REGION','NEAREST_GENE','DISTANCE')]

# write results
message(sprintf('Writing %s.',outputFile))
write.table(suppl, file = outputFile, row.names = F, quote = F, sep = '\t', na = "NA")
system(paste0('chmod 770 ', outputFile))
message('-- Completed: Creating table of SusieR fine-mapping results for supplementum ---')
