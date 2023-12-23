#!/usr/bin/env Rscript

# ==================================================================================
# === Create table of nonsynonymous variants in discovered loci for supplementum ===
# ==================================================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop(paste0('expected 5 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
snplevelFile=args[1] # snplevelFile="results/combined/snplevel.txt"
nsynFiles=args[2] # nsynFiles="results/gap_gm/conditional/conditional.cleaned.tophits.annovar.nonsynonymous.txt,results/gap_wm/conditional/conditional.cleaned.tophits.annovar.nonsynonymous.txt,results/gap_gwm/conditional/conditional.cleaned.tophits.annovar.nonsynonymous.txt"
outputFile=args[3] # outputFile="results/combined/nonsynonymous.suppl.txt"
traits=args[4] # traits="gap_gm,gap_wm,gap_gwm"
traitNames=args[5] # traitNames="grey matter,white matter,grey and white matter"

message(paste0('\n--- Create table of exonic non-synonymous variants for supplementum ---',
               '\nsnplevelFile: ', snplevelFile,
               '\nnsynFiles: ', nsynFiles,
               '\noutputFile: ', outputFile,
               '\ntraits: ', traits,
               '\ntraitNames: ', traitNames, '\n'))

# attach packages to current R session
for (pkg in c('dplyr', 'stringr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# transform variables
traits = str_split(traits, ',')[[1]]
traitNames = str_split(traitNames, ',')[[1]]
nsynFiles = str_split(nsynFiles, ',')[[1]]

# load data
message('Loading Files.')
snplevel = read.delim(snplevelFile, sep = '\t', header = T)
snplevel$key = paste0(snplevel$TRAIT,"_", snplevel$ID)
snplevel = snplevel[,c('LOCUS_COUNT','DISCOV_COUNT','CHR', 'BP','key')]
names(snplevel) = c('LOCUS_COUNT','DISCOV_COUNT','LEAD_SNP_CHR', 'LEAD_SNP_BP','key')

# load nsyn Files
for (i in 1:length(nsynFiles)) {
  tmp = read.table(nsynFiles[i], sep = '\t', head = T)
  tmp$trait = traitNames[i]
  tmp$key = paste0(traits[i],"_", tmp$LEAD_SNP)
  if (i == 1) { nsyn = tmp } else { nsyn = rbind(nsyn, tmp)}
}

# combine
df = inner_join(nsyn, snplevel, by = 'key')
df = df[order(df$LOCUS_COUNT, df$DISCOV_COUNT, df$P),]

# create supplementum data.frame
suppl = df[,c('LOCUS_COUNT','trait','LEAD_SNP','LEAD_SNP_CHR','LEAD_SNP_BP','LEAD_SNP_RSQ',
              'ID','BP','TYPED','INFO','A1','A2','A1_FREQ','BETA','SE','Z','P','ETA2','N',
              'REGION','NEAREST_GENE','NEAREST_GENE_DESCRIPTION','NEAREST_GENE_BIOTYPE',
              'EXONIC_FUNCTION','TRANSCRIPT_CONSEQUENCE','CADD_PHRED','CADD_RAW','CADD_RANK','DANN','DANN_RANK','REVEL','REVEL_RANK')]

# write results
message('Writing output file.')
write.table(suppl, file = outputFile, row.names = F, quote = F, sep = '\t', na = "NA")
system(paste0('chmod 770 ', outputFile))
message('-- Creating suppl. table of exonic non-synonymous variants completed. ---')
