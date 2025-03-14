#!/usr/bin/env Rscript

# ==========================================================
# === create phesant summary text file and combined plot ===
# ==========================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop(paste0('expected 4 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
traits=args[1] # traits="rc_auc,rc_phi,rc_range"
phesantSummary=args[2] # phesantSummary="results/rc_auc/phesant/phewas_summary.txt,results/rc_phi/phesant/phewas_summary.txt,results/rc_range/phesant/phewas_summary.txt"
phesantPlot=args[3] # phesantPlot="results/rc_auc/phesant/phewas.png,results/rc_phi/phesant/phewas.png,results/rc_range/phesant/phewas.png"
outputFile=args[4] # outputFile="results/combined/suppl.phewas"

message(paste0('\n--- Create phesant summary text file and combined plot ---',
               '\ntraits: ', traits,
               '\nphesantSummary: ', phesantSummary,
               '\nphesantPlot: ', phesantPlot,
               '\noutputFile: ', outputFile,'\n'))

# attach packages to current R session
for (pkg in c('data.table', 'dplyr', 'ggplot2', 'ggpubr', 'plotly', 'patchwork','magick','stringr')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# transform variables
phesantSummary = str_split(phesantSummary, ',')[[1]]
phesantPlot = str_split(phesantPlot, ',')[[1]]
traits = str_split(traits, ',')[[1]]

# open and merge datasets
for (i in 1:length(traits)) {
  
  # open datasets and keep relevant variables
  tmp = read.delim(phesantSummary[i], sep = '\t', header = T)
  if (i ==1) {
    tmp = tmp[,c('description', 'varName', 'custom_category', 'Path', 'varType', 'resType', 'n', 'ntotal', 'beta', 'lower', 'upper', 'se', 'pvalue', 'fdr', 'z', 'rho')]
  } else { 
    tmp = tmp[,c('varName', 'beta', 'lower', 'upper', 'se', 'pvalue', 'fdr', 'z', 'rho')]
  }
  tmp$rhoAbs = abs(tmp$rho)
  
  # create unique variable names
  names(tmp)[which(names(tmp) %in% c('beta', 'lower', 'upper', 'se', 'pvalue', 'fdr', 'z', 'rho', 'rhoAbs'))] = 
    paste0(traits[i],c('_beta', '_lower', '_upper', '_se', '_pvalue', '_fdr', '_z', '_rho', '_rhoAbs'))
  
  # merge datasets
  if (i ==1) { df = tmp } else { df = left_join(df,tmp, by = 'varName') }
}
  
# get top p-value and top |rho|
df$top_pvalue =  df[,grep('pvalue',names(df))] %>% apply(1, FUN = min)
df$top_fdr = df[,grep('fdr',names(df))] %>% apply(1, FUN = min)
df$top_rhoAbs = df[,grep('rhoAbs',names(df))] %>% apply(1, FUN = max)

# create output
message('Writing txt file.')
df$`NA` = ""
cols = c('varName', 'description', 'ntotal', 'top_rhoAbs','top_pvalue','top_fdr')
for (i in 1:length(traits)) {
  cols = c(cols,'NA',paste0(traits[i],c('_beta', '_se', '_rho', '_pvalue', '_fdr')))
}
cols = c(cols, c('NA', 'n', 'varType', 'resType', 'custom_category', 'Path'))
output = df[,cols]
output = output[order(output$top_fdr,output$top_pvalue),]
write.table(output, paste0(outputFile,'.txt'), sep = '\t', quote = F, row.names = F)

# combine phesant plots
for (i in 1:length(traits)) {
  tmp = image_read(phesantPlot[i])
  tmp = ggplot() + background_image(tmp) + coord_fixed(ratio = image_info(tmp)$height/image_info(tmp)$width)
  if (i ==1) { plot = tmp } else { plot = plot / tmp }
}

# draw plot with annotations
plot = plot + plot_annotation(tag_levels ='a') & 
  theme(plot.tag.position = c(0, 1), plot.tag = element_text(size = 22))

# save plot
message('Writing png file.')
png(width = 8.7, height = 5.0*3, units = "in", res = 600, filename = paste0(outputFile,'.png'))
plot
invisible(dev.off())
system(paste0('chmod 770 ', outputFile, '*'))
message('-- Completed: Creating suppl. table of PheWAS results ---')
