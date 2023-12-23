#!/usr/bin/env Rscript

# ===========================================================
# === create combined rg plot (selected traits + volcano) ===
# ===========================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop(paste0('expected 3 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
rgSelection = args[1] # rgSelection="results/combined/rgSelection.png"
rgNeale = args[2] # rgNeale="results/combined/rgNeale.volcano.forest.png"
outFile = args[3] # outFile="results/combined/rgCombined.png"

message(paste0('\n--- Combine rgSelection and rgNeale plots ---',
               '\nrgSelection: ', rgSelection,
               '\nrgNeale: ', rgNeale,
               '\noutFile: ', outFile,'\n'))

# attach required packages
for (pkg in c('ggpubr','ggplot2','patchwork','magick','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# load rgSelection
message(sprintf('[1/4] Loading %s', rgSelection))
plotSelect = image_read(rgSelection)
plotSelect.width = image_info(plotSelect)$width
plotSelect.height = image_info(plotSelect)$height
plotSelect = ggplot() +
background_image(plotSelect) + coord_fixed(ratio = image_info(plotSelect)$height/image_info(plotSelect)$width) &
  theme(plot.tag.position  = c(.2, .98)) 

# load rgSelection
message(sprintf('[2/4] Loading %s', rgNeale))
plotNeale = image_read(rgNeale)
plotNeale.width = image_info(plotNeale)$width
plotNeale.height = image_info(plotNeale)$height
plotNeale = ggplot() +
background_image(plotNeale) + coord_fixed(ratio = image_info(plotNeale)$height/image_info(plotNeale)$width) &
  theme(plot.tag.position  = c(.000, 1.36)) 

# combine plots
message('[3/4] Combining plots.')
tmp = plot_spacer() + plotSelect + plot_spacer() + plotNeale + plot_layout(widths = c(-400,plotSelect.width,-400, plotNeale.width))
pl = tmp + plot_annotation(tag_levels = list(c('a','b\n\n\n\n\n\n\n\n\n\n\n\n\n  c'))) &
  theme(plot.tag = element_text(size = 13, face = 'bold', hjust = 0, vjust = 1.9))

# save file
message(sprintf('[4/4] Saving %s',outFile))
png(width = 7.79, height = 6.20, units = "in", res = 600, filename = outFile)
pl
invisible(dev.off())

