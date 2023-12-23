#!/usr/bin/env Rscript

# ===========================================================
# === Create table of colocalized GTEx multi-tissue eQTLs ===
# ===========================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop(paste0('expected 5 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
snplevelFile=args[1] # snplevelFile="results/combined/snplevel.txt"
eqtlFiles=args[2] # eqtlFiles="results/gap_gm/eqtl/eqtl.multiTissue.txt,results/gap_wm/eqtl/eqtl.multiTissue.txt,results/gap_gwm/eqtl/eqtl.multiTissue.txt"
outputFile=args[3] # outputFile="results/combined/gtex.multitissue.txt"
traits=args[4] # traits="gap_gm,gap_wm,gap_gwm"
traitNames=args[5] # traitNames="grey matter,white matter,grey and white matter"

message(paste0('\n--- Create table of colocalized GTEx multi-tissue eQTLs ---',
               '\nsnplevelFile: ', snplevelFile,
               '\neqtlFiles: ', eqtlFiles,
               '\noutputFile: ', outputFile,
               '\ntraits: ', traits,
               '\ntraitNames: ', traitNames, '\n'))

# attach packages to current R session
for (pkg in c('dplyr', 'stringr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# transform variables
traits = str_split(traits, ',')[[1]]
traitNames = str_split(traitNames, ',')[[1]]
eqtlFiles = str_split(eqtlFiles, ',')[[1]]

# load data
message('Loading Files.')
snplevel = read.delim(snplevelFile, sep = '\t', header = T)
snplevel$key = paste0(snplevel$TRAIT,"_", snplevel$ID)
snplevel = snplevel[,c('LOCUS_COUNT','DISCOV_COUNT','CHR', 'BP','key')]

# load eqtl files
for (i in 1:length(eqtlFiles)) {
  tmp = read.table(eqtlFiles[i], sep = '\t', head = T)
  tmp = tmp[,-which(names(tmp) == "LOCUS_COUNT")]
  #tmp = tmp[,-grep("pval_",names(tmp))]
  tmp$trait = traitNames[i]
  tmp$key = paste0(traits[i],"_", tmp$LEAD_SNP)
  if (i == 1) { eqtl = tmp } else { eqtl = rbind(eqtl, tmp)}
}

# combine
df = inner_join(eqtl, snplevel, by = 'key')
df = df[order(df$LOCUS_COUNT, df$DISCOV_COUNT, df$PVALUE_RE2),]
df$X.SIGN = rowSums(df[,names(df)[grep('mval_',names(df))]] >= 0.9, na.rm = TRUE)

# create supplementum data.frame
suppl = df[,c('LOCUS_COUNT','trait','LEAD_SNP', 'CHR','BP','id','LEAD_SNP_R2',
              'RSID','gene_id','hgnc_symbol','X.STUDY','X.SIGN',
              'BETA_FE','STD_FE','PVALUE_FE','STAT1_RE2','STAT2_RE2','PVALUE_RE2',
              names(df)[grep('mval_',names(df))])]
suppl$RSID = sub('\\,.*', '', suppl$RSID)

# write results
message('Writing output file.')
write.table(suppl, file = outputFile, row.names = F, quote = F, sep = '\t', na = "NA")
system(paste0('chmod 770 ', outputFile))
message('-- Creating table of colocalized GTEx multi-tissue eQTLs completed. ---')
