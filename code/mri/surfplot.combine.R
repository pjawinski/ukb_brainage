#!/usr/bin/env Rscript

# =============================
# === combine surface plots ===
# =============================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop(paste0('expected 5 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
traits = args[1] # traits="gap_gm,gap_wm,gap_gwm"
surfPlots = args[2] # surfPlots="results/gap_gm/surfplot/surfplot.png,results/gap_wm/surfplot/surfplot.png,results/gap_gwm/surfplot/surfplot.png"
outFile = args[3] # outFile="results/combined/surfplot.png"
width = as.numeric(args[4]) # width = 10
height = as.numeric(args[5]) # height = 8.85

message(paste0('\n--- Settings: Combine surface plots ---',
               '\ntraits: ', traits,
               '\nsurfPlots: ', surfPlots,
               '\noutFile: ', outFile,
               '\nwidth: ', width,
               '\nheight: ', height,'\n'))

# attach required packages
for (pkg in c('dplyr', 'ggpubr','ggplot2','patchwork','magick','stringr')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# transform variables
traits = str_split(traits, ',')[[1]]
surfPlots = str_split(surfPlots, ',')[[1]]

# load images
for (i in 1:length(traits)) {
  
  # load surfplots
  message(paste0('[', i, '/', length(traits), '] Loading ', surfPlots[i]))
  surfplot = image_read(surfPlots[i])
  surfplot = image_trim(surfplot)
  surfplot.width = image_info(surfplot)$width
  surfplot.height = image_info(surfplot)$height
  surfplot = ggplot() +
  annotation_raster(surfplot, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  ylim(c(-0.1,1)) +
  coord_fixed(ratio = surfplot.height/surfplot.width) +
  theme_void() &
  theme(plot.tag.position  = c(0.05, 0.95)) 

  #surfplot = ggplot() +
  #background_image(surfplot) + coord_fixed(ratio = image_info(surfplot)$height/image_info(surfplot)$width)
  #assign(paste0(traits[i],"_surfplot"), surfplot)
  
  # collect plots
  tmp = surfplot + plot_layout(widths = c(surfplot.width))
  if (i == 1) { pl = tmp } else { pl = pl + tmp }
}

# set layout
message('Writing png file.')
if (length(traits) > 1) {
  pl = pl + plot_annotation(tag_levels = list(letters)) + plot_layout(ncol = 2) & # list(c(rbind(letters,' ')))
    theme(plot.tag = element_text(size = 22))
}

# save png
message('Saving plot.')
png(width = width, height = height, units = "in", res = 300, filename = outFile)
pl
invisible(dev.off())
message('--- Completed: combine surface plots ---')


