#!/usr/bin/env Rscript

# ======================================
# === create combined manhattan plot ===
# ======================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=7) {
  stop(paste0('expected 7 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
traits = args[1] # traits="gap_gm,gap_wm,gap_gwm"
plotTitles = args[2] # plotTitles="Grey_matter,White_matter,Grey_and_white_matter"
ldscPlots = args[3] # ldscPlots="results/gap_gm/ldsc/ldsc.partitioned.results.png,results/gap_wm/ldsc/ldsc.partitioned.results.png,results/gap_gwm/ldsc/ldsc.partitioned.results.png"
ldscSummary = args[4] # ldscSummary="results/gap_gm/ldsc/ldsc.partitioned.results.summary.txt,results/gap_wm/ldsc/ldsc.partitioned.results.summary.txt,results/gap_gwm/ldsc/ldsc.partitioned.results.summary.txt"
outputFile = args[5] # outputFile="results/combined/ldsc.partitioned"
width = as.numeric(args[6]) # width = 14
height = as.numeric(args[7]) # height = 17.5

message(paste0('\n--- Completed: Combine partitioned ldsc results (figures and summary tables) ---',
               '\ntraits: ', traits,
               '\nplotTitles: ', plotTitles,
               '\nldscPlots: ', ldscPlots,
               '\nldscSummary: ', ldscSummary,
               '\noutputFile: ', outputFile,
               '\nwidth: ', width,
               '\nheight: ', height,'\n'))

# attach required packages
for (pkg in c('dplyr','ggpubr','ggplot2','patchwork','magick','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# transform variables
traits = str_split(traits, ',')[[1]]
plotTitles = str_split(plotTitles, ',')[[1]]
plotTitles = str_replace_all(plotTitles, "_", " ")
ldscPlots = str_split(ldscPlots, ',')[[1]]
ldscSummary = str_split(ldscSummary, ',')[[1]]

# open and merge datasets
for (i in 1:length(traits)) {
  
  # open datasets and keep relevant variables
  message(paste0('Loading ', ldscSummary[i]))
  tmp = read.delim(ldscSummary[i], sep = '\t', header = T)
  if (i ==1) {
    tmp = tmp[,c('annotation', 'Prop._SNPs', 'Prop._h2', 'Enrichment', 'Enrichment_p', 'Enrichment_FDR')]
  } else { 
    tmp = tmp[,c('annotation', 'Prop._h2', 'Enrichment', 'Enrichment_p', 'Enrichment_FDR')]
  }
  
  # create unique variable names
  names(tmp)[which(names(tmp) %in% c('Prop._h2', 'Enrichment', 'Enrichment_p', 'Enrichment_FDR'))] = 
    paste0(traits[i],c('_Prop._h2', '_Enrichment', '_Enrichment_p', '_Enrichment_FDR'))
  
  # merge datasets
  if (i ==1) { df = tmp } else { df = left_join(df,tmp, by = 'annotation') }
}
  
# get top p-value and top |rho|
df$top_Enrichment =  df[,grep('Enrichment',names(df))] %>% apply(1, FUN = max)
df$top_Enrichment_p = df[,grep('Enrichment_p',names(df))] %>% apply(1, FUN = min)
df$top_Enrichment_FDR = df[,grep('Enrichment_FDR',names(df))] %>% apply(1, FUN = min)
df = df[order(df$top_Enrichment_FDR,df$top_Enrichment_p,-df$top_Enrichment),]

# create output
message('Writing txt file.')
df$`NA` = ""
cols = c('annotation', 'Prop._SNPs', 'top_Enrichment','top_Enrichment_p','top_Enrichment_FDR')
for (i in 1:length(traits)) {
  cols = c(cols,'NA',paste0(traits[i],c('_Prop._h2', '_Enrichment', '_Enrichment_p','_Enrichment_FDR')))
}
output = df[,cols]
write.table(output, paste0(outputFile,'.txt'), sep = '\t', quote = F, row.names = F)

# load images
for (i in 1:length(traits)) {
  
  # load manhattan plots
  message(paste0('Loading ', ldscPlots[i]))
  ldsc = image_read(ldscPlots[i])
  ldsc.width = image_info(ldsc)$width
  ldsc.height = image_info(ldsc)$height
  ldsc = ggplot() +
  background_image(ldsc) + coord_fixed(ratio = image_info(ldsc)$height/image_info(ldsc)$width)
  assign(paste0(traits[i],"_manh"), ldsc)
  
  # collect plots
  tmp = ldsc
  if (i == 1) { pl = tmp } else { pl = pl / ldsc }
}

# # set layout
# if (length(traits) > 1) {
#   pl = pl + plot_annotation(tag_levels = list(c(rbind(letters[1:length(traits)],letters[(1+length(traits)):(2*length(traits))])))) & # list(c(rbind(letters,' ')))
#     theme(plot.tag.position = c(0, 1),
#           plot.tag = element_text(size = 22))
# }

# set layout
if (length(plotTitles) > 1) {
  pl = pl + plot_annotation(tag_levels = list(c(plotTitles))) &
    theme(plot.tag.position = c(0.52, 1),
          plot.tag = element_text(size = 12, face = 'bold', hjust = 0.5, vjust = 1.9))
}

# save file
message('Saving plot.')
png(width = width, height = height, units = "in", res = 300, filename = paste0(outputFile,'.png'))
pl
invisible(dev.off())
message('--- Completed: Combine partitioned ldsc results (figures and summary tables) ---')

