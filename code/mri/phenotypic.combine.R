#!/usr/bin/env Rscript

# ===================================================
# === combine accuracy, phesant, and surface plot ===
# ===================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop(paste0('expected 4 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
accuracyplot = args[1] # accuracyplot="results/mri/accuracy.gwm.png"
phesantplot = args[2] # phesantplot="results/gap_gwm/phesant/phewas.png"
surfplot = args[3] # surfplot="results/gap_gwm/surfplot/surfplot.png"
outFile = args[4] # outFile="results/combined/phenotypic.png"

message(paste0('\n--- Combine phesant and surface plot ---',
               '\naccuracyplot: ', accuracyplot,
               '\nphesantplot: ', phesantplot,
               '\nsurfplot: ', surfplot,
               '\noutFile: ', outFile,'\n'))

# attach required packages
for (pkg in c('ggpubr','ggplot2','patchwork','magick','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# load phesant
message(sprintf('[1/5] Loading %s', phesantplot))
accuracy = image_read(accuracyplot)
accuracy = image_trim(accuracy)
accuracy.width = image_info(accuracy)$width
accuracy.height = image_info(accuracy)$height
accuracy = ggplot() +
  background_image(accuracy) + coord_fixed(ratio = image_info(accuracy)$height/image_info(accuracy)$width) &
  theme(plot.tag.position  = c(-0.02, 0.98)) 

# load phesant
message(sprintf('[2/5] Loading %s', phesantplot))
phesant = image_read(phesantplot)
phesant = image_trim(phesant)
phesant.width = image_info(phesant)$width
phesant.height = image_info(phesant)$height
phesant = ggplot() +
  background_image(phesant) + coord_fixed(ratio = image_info(phesant)$height/image_info(phesant)$width) &
  theme(plot.tag.position  = c(0.02, 1)) 

# load surfplot
message(sprintf('[3/5] Loading %s', surfplot))
surf = image_read(surfplot)
surf = image_trim(surf)
surf.width = image_info(surf)$width
surf.height = image_info(surf)$height
surf = ggplot() +
  annotation_raster(surf, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  ylim(c(0,0.97)) +
  coord_fixed(ratio = image_info(surf)$height/image_info(surf)$width) +
  theme_void() &
  theme(plot.tag.position  = c(0.03, 1)) 

# combine plots
message('[4/5] Combining plots.')
layout <- c(
  area(t = 1, b = 50,  l = 3, r = 98),
  area(t = 54, b = 110,  l = 1, r = 71),
  area(t = 54, b = 103, l = 67, r = 102)
)
tmp = accuracy + phesant + surf + 
  plot_layout(design = layout)

pl = tmp + plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 14, face = 'bold', hjust = 0, vjust = 0.5))

# save file
message(sprintf('[5/5] Saving %s',outFile))
png(width = 8.1, height = 5.7, units = "in", res = 600, filename = outFile)
pl
invisible(dev.off())
message('--- Completed: Combine phesant and surface plot ---')

