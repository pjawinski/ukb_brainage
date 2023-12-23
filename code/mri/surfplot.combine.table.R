# ===========================================================================
# === create phesant summary text file and combined plot for supplementum ===
# ===========================================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop(paste0('expected 4 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
traits=args[1] # traits="gap_gm,gap_wm,gap_gwm"
surfplotTXT=args[2] # surfplotTXT="results/gap_gm/surfplot/surfplot.txt,results/gap_wm/surfplot/surfplot.txt,results/gap_gwm/surfplot/surfplot.txt"
phesantPlot=args[3] # phesantPlot="results/rc_auc/phesant/phewas.png,results/rc_phi/phesant/phewas.png,results/rc_range/phesant/phewas.png"
outputFile=args[4] # outputFile="results/combined/suppl.phewas"

message(paste0('\n--- Settings ---',
               '\ntraits: ', traits,
               '\nsurfplotTXT: ', surfplotTXT,
               '\nphesantPlot: ', phesantPlot,
               '\noutputFile: ', outputFile,'\n'))

# attach packages to current R session
for (pkg in c('data.table', 'dplyr', 'ggplot2', 'ggpubr', 'plotly', 'patchwork','magick','stringr')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# transform variables
surfplotTXT = str_split(surfplotTXT, ',')[[1]]
phesantPlot = str_split(phesantPlot, ',')[[1]]
traits = str_split(traits, ',')[[1]]

# open and merge datasets
for (i in 1:length(traits)) {
  
  # open datasets and keep relevant variables
  tmp = read.delim(surfplotTXT[i], sep = '\t', header = T)
  tmp$fdr = p.adjust(tmp$pval, method = 'BH')
  if (i == 1) {
    tmp = tmp[,c('id', 'n', 'pval', 'fdr', 'rho', 'rhoabs')]
  } else { 
    tmp = tmp[,c('id', 'pval', 'fdr', 'rho', 'rhoabs')]
  }

  # create unique variable names
  names(tmp)[which(names(tmp) %in% c('pval', 'fdr', 'rho', 'rhoabs'))] = 
    paste0(traits[i],c('_pval', '_fdr', '_rho', '_rhoAbs'))
  
  # merge datasets
  if (i == 1) { df = tmp } else { df = left_join(df,tmp, by = 'id') }
}
  
# get top p-value and top |rho|
df$top_pval =  df[,grep('pval',names(df))] %>% apply(1, FUN = min)
df$top_fdr = df[,grep('fdr',names(df))] %>% apply(1, FUN = min)
df$top_rhoAbs = df[,grep('rhoAbs',names(df))] %>% apply(1, FUN = max)

# get type (cortical / subcortical)'
df$type = 'subcortical'
df$type[grep('ThickAvg',df$id)] = 'cortical'
df$measure = 'volume'
df$measure[grep('ThickAvg',df$id)] = 'cortical thickness'

# create structure column
df$structure = df$id %>% 
  str_replace_all('\\.', ' ') %>%
  str_replace_all('RH_', 'Right ') %>%
  str_replace_all('LH_', 'Left ') %>%
  str_replace_all('ThickAvg_', '') %>%
  tolower() %>%
  str_replace('^\\w{1}', toupper)

# create output
message('Writing txt file.')
df$`NA` = ""
cols = c('structure','type', 'measure', 'n','top_pval','top_fdr','top_rhoAbs')
for (i in 1:length(traits)) {
  cols = c(cols,'NA',paste0(traits[i],c('_pval', '_fdr', '_rho', '_rhoAbs')))
}
output = df[,cols]
output = output[order(output$top_fdr,output$top_pval),]
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
message('-- Creating suppl. table of PheWAS results completed. ---')
