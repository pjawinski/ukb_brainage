#!/usr/bin/env Rscript

# ===================================================================
# === Create table of discoveries for main article & supplementum ===
# ===================================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=7) {
  stop(paste0('expected 7 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
noveltyFile=args[1] # noveltyFile="results/combined/snplevel.novelty.txt"
replicationFile=args[2] # replicationFile="results/combined/replicateDiscoveries.txt"
outputFileMain=args[3] # outputFileMain="results/combined/discoveries.main.txt"
outputFileSuppl=args[4] # outputFileSuppl="results/combined/discoveries.suppl.txt"
traits=args[5] # traits="gap_gm,gap_wm,gap_gwm"
traitNamesMain=args[6] # traitNamesMain="GM,WM,GWM"
traitNamesSuppl=args[7] # traitNamesSuppl="grey matter,white matter,grey and white matter"

message(paste0('\n--- Creating table of discoveries for main article & supplementum ---',
               '\nnoveltyFile: ', noveltyFile,
               '\nreplicationFile: ', replicationFile,
               '\noutputFileMain: ', outputFileMain,
               '\noutputFileSuppl: ', outputFileSuppl,
               '\ntraits: ', traits,
               '\ntraitNamesMain: ', traitNamesMain,
               '\ntraitNamesSuppl: ', traitNamesSuppl, '\n'))

# attach packages to current R session
for (pkg in c('dplyr', 'stringr')) { eval(bquote(suppressPackageStartupMessages(require(.(pkg))))) }

# load data
message('Loading Files.')
df = read.delim(noveltyFile, sep = '\t', header = T)
if (replicationFile!="") { repl = read.delim(replicationFile, sep = '\t', header = T) }

# transform variables
traits = str_split(traits, ',')[[1]]
traitNamesMain = str_split(traitNamesMain, ',')[[1]]
traitNamesSuppl = str_split(traitNamesSuppl, ',')[[1]]

# add columns
df$BPstr = format(df$BP, big.mark=",")
df$IDtrunc = gsub("_.*","",df$ID)
df$A1A2 = paste0(str_trunc(df$A1,3,'right', ellipsis=""),'/',str_trunc(df$A2,3,'right', ellipsis=""))
df$A1_FREQround = paste0(format(round(df$A1_FREQ,2), nsmall = 2))
df$BETASE = paste0(format(round(df$BETA,2), nsmall = 2), ' (', format(round(df$SE,2), nsmall = 2), ')')
df$Pstr = paste0(formatC(df$P,0))
df$GeneDistance = paste0(df$NEAREST_GENE)
df$GeneDistance[df$DISTANCE>0] = paste0(df$NEAREST_GENE[df$DISTANCE>0], ' (', round(df$DISTANCE[df$DISTANCE>0]/1000,0), 'kb)')
df$TraitNamesMain = factor(df$TRAIT, levels = traits, labels = traitNamesMain)
df$TraitNamesSuppl = factor(df$TRAIT, levels = traits, labels = traitNamesSuppl)

df$TraitsPerLocus = ""
for (i in 1:nrow(df)) {
  tmp = df[df$LOCUS_COUNT == df$LOCUS_COUNT[i],]
  tmp = tmp[!duplicated(tmp$TraitNamesMain),]
  df$TraitsPerLocus[i] = paste0(tmp$TraitNamesMain, collapse = ',')
}

# merge SBayesRC and SusieR results
df$cssizes = df$sbayes_cssizes
idx = df$CHR == 'X' | df$CHR== 'XY' | df$CHR == 'Y'
df$cssizes[idx] = df$susieR_cssize[idx]

# merge with replication
if (replicationFile!="") {
  repl = repl[,c('DISCOV_TRAITNAME','DISCOV_ID','REPLIC_CONSIST_DIRECTION','REPLIC_CONSIST_REPLICATED')]
  dfmerge = left_join(df,repl,by = c('TraitNamesSuppl' = 'DISCOV_TRAITNAME', 'ID' = 'DISCOV_ID'))

  # add some further columns for output
  dfmerge$ReplMain = 'F' # '-'
  #dfmerge$ReplMain[dfmerge$REPLIC_CONSIST_DIRECTION=="TRUE"] = '+' # '+'
  dfmerge$ReplMain[dfmerge$REPLIC_CONSIST_REPLICATED=="TRUE"] = 'T' # '*'
  dfmerge$LITERATURE_SHORT[df$LITERATURE_SHORT=="-"] = '*'

  # create main data.frame
  main = dfmerge[!duplicated(df$LOCUS_COUNT),c('CHR','BPstr', 'IDtrunc','A1A2','A1_FREQround', 'BETASE','Pstr',
                'susieR_cssizes','prioritized','TraitsPerLocus','ReplMain','LITERATURE_SHORT')]
} else {
  dfmerge = df
  dfmerge$LITERATURE_SHORT[df$LITERATURE_SHORT=="-"] = '*'
  main = dfmerge[!duplicated(df$LOCUS_COUNT),c('CHR','BPstr', 'IDtrunc','A1A2','A1_FREQround', 'BETASE','Pstr',
                'cssizes','prioritized','TraitsPerLocus','LITERATURE_SHORT')]
}

# create supplementum data.frame
dfmerge$`NA` = ""
if (!("Z" %in% names(dfmerge))) { dfmerge$Z = dfmerge$BETA / dfmerge$SE }
if (!("ETA2" %in% names(dfmerge))) { dfmerge$ETA2 = dfmerge$Z^2 / (dfmerge$Z^2+dfmerge$N-2) } # approximate estimate for eta2, neglecting covariates

suppl = dfmerge[,c('LOCUS_COUNT','TraitNamesSuppl','CYTOBAND','CHR','BP',
              'ID','A1','A2','A1_FREQ',
              'BETA','SE','Z','P','ETA2','N',
              'REGION','NEAREST_GENE','NEAREST_GENE_DESCRIPTION','NEAREST_GENE_BIOTYPE','DISTANCE',
              'sbayes_cssizes','susieR_cssizes','finemap_cssize',
              'NA','sbayes_genes','nonsyn_customPthresh', 'SMR_eqtl','SMR_sqtl','GTEx_singleTissue','GTEx_multiTissue','PoPS','prioritized',
              'NA','GWAS_CATALOG','LITERATURE_LONG')]

# write results
message('Writing output file.')
write.table(main, file = outputFileMain, row.names = F, quote = F, sep = '\t', na = "-")
write.table(suppl, file = outputFileSuppl, row.names = F, quote = F, sep = '\t', na = "-")
system(paste0('chmod 770 ', outputFileMain))
system(paste0('chmod 770 ', outputFileSuppl))
message('-- Creating table of discoveries completed. ---')
