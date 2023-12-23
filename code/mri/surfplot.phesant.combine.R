#!/usr/bin/env Rscript

# ========================================
# === combine phesant and surface plot ===
# ========================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop(paste0('expected 3 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
phesantplot = args[1] # phesantplot="results/gap_gwm/phesant/phewas.png"
surfplot = args[2] # surfplot="results/gap_gwm/surfplot/surfplot.png"
outFile = args[3] # outFile="results/combined/surfplot.phesant.png"

message(paste0('\n--- Combine phesant and surface plot ---',
               '\nphesantplot: ', phesantplot,
               '\nsurfplot: ', surfplot,
               '\noutFile: ', outFile,'\n'))

# attach required packages
for (pkg in c('ggpubr','ggplot2','patchwork','magick','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# load rgSelection
message(sprintf('[1/4] Loading %s', phesantplot))
phesant = image_read(phesantplot)
phesant.width = image_info(phesant)$width
phesant.height = image_info(phesant)$height
phesant = ggplot() +
background_image(phesant) + coord_fixed(ratio = image_info(phesant)$height/image_info(phesant)$width) &
  theme(plot.tag.position  = c(.14, 1)) 

# load rgSelection
message(sprintf('[2/4] Loading %s', surfplot))
surf = image_read(surfplot)
surf.width = image_info(surf)$width
surf.height = image_info(surf)$height
surf = ggplot() +
  annotation_raster(surf, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  ylim(c(-0.17,0.94)) +
  coord_fixed(ratio = image_info(surf)$height/image_info(surf)$width) +
  theme_void() &
  theme(plot.tag.position  = c(0.05, 1)) 

# combine plots
message('[3/4] Combining plots.')
tmp = plot_spacer() + phesant + plot_spacer() + surf + plot_spacer() + plot_layout(widths = c(-925,phesant.width,-1000, surf.width,-650))
#tmp = phesant + surf + plot_layout(widths = c(phesant.width,surf.width))
pl = tmp + plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 13, face = 'bold', hjust = 0, vjust = 0.5))

# save file
message(sprintf('[4/4] Saving %s',outFile))
png(width = 8.1, height = 3.31, units = "in", res = 600, filename = outFile)
pl
invisible(dev.off())
message('--- Completed: Combine phesant and surface plot ---')
               
