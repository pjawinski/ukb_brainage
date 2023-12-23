#!/usr/bin/env Rscript

# ======================================
# === create combined manhattan plot ===
# ======================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop(paste0('expected 5 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
plotTitles = args[1] # plotTitles="Grey_matter,White_matter,Grey_and_white_matter,Cross-trait"
plots = args[2] # plots="results/combined/replicateDiscoveries.gm.png,results/combined/replicateDiscoveries.wm.png,results/combined/replicateDiscoveries.gwm.png,results/combined/replicateDiscoveries.crosstrait.png"
outFile = args[3] # outFile="results/combined/replicateDiscoveries.qq.png"
width = as.numeric(args[4]) # width = 10
height = as.numeric(args[5]) # height = 4

message(paste0('\n--- Combine plots ---',
               '\nplotTitles: ', plotTitles,
               '\nplots: ', plots,
               '\noutFile: ', outFile,
               '\nwidth: ', width,
               '\nheight: ', height,'\n'))

# attach required packages
for (pkg in c('ggpubr','ggplot2','patchwork','magick','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# transform variables
plotTitles = str_split(plotTitles, ',')[[1]]
plotTitles = str_replace_all(plotTitles, "_", " ")
plots = str_split(plots, ',')[[1]]

# load images
for (i in 1:length(plots)) {
  
  # load manhattan plots
  message(sprintf('[%d/%d] Loading %s',i,length(plots),plots[i]))
  tmp = image_read(plots[i])
  tmp.width = image_info(tmp)$width
  tmp.height = image_info(tmp)$height
  tmp = ggplot() +
  background_image(tmp) + coord_fixed(ratio = image_info(tmp)$height/image_info(tmp)$width)
  assign(sprintf('pl%d',i), tmp)

  # collect plots
  if (i == 1) { pl = tmp } else { pl = pl + tmp }
}

# set layout
if (length(plotTitles) == 1) {
  pl = pl + plot_layout(nrow = 1) + plot_annotation(tag_levels = list(c(letters))) & # list(c(rbind(letters,' ')))
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(size = 22))
} else {
  pl = pl + plot_layout(nrow = 1) + plot_annotation(tag_levels = list(plotTitles)) &
    theme(plot.tag.position = c(0.58, 1),
          plot.tag = element_text(size = 12, face = 'bold', hjust = 0.5, vjust = 0.4))
}

# save file
png(width = width, height = height, units = "in", res = 300, filename = outFile)
pl
invisible(dev.off())
message('--- Completed: Combine plots ---')
