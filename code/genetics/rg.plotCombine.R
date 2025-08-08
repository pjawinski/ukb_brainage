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
rgSelection = args[1] # rgSelection="results/combined/rgSelection"
rgNeale = args[2] # rgNeale="results/combined/rgNeale.volcano.forest"
outFile = args[3] # outFile="results/combined/rgCombined"

message(paste0('\n--- Combine rgSelection and rgNeale plots ---',
               '\nrgSelection: ', rgSelection,
               '\nrgNeale: ', rgNeale,
               '\noutFile: ', outFile,'\n'))

# attach required packages
for (pkg in c('ggpubr','ggplot2','grid','grImport2','patchwork','magick','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# initialize plots
for (fileFormat in c(".png",".svg")) {

  if (fileFormat == ".png") {

    # load rgSelection
    message(sprintf('[PNG 1/4] Loading %s.png', rgSelection))
    plotSelect = image_read(sprintf('%s.png',rgSelection))
    plotSelect.width = image_info(plotSelect)$width
    plotSelect.height = image_info(plotSelect)$height
    plotSelect = ggplot() +
    background_image(plotSelect) + coord_fixed(ratio = image_info(plotSelect)$height/image_info(plotSelect)$width) &
      theme(plot.tag.position  = c(.2, .98)) 

    # load rgNeale
    message(sprintf('[PNG 2/4] Loading %s.png', rgNeale))
    plotNeale = image_read(sprintf('%s.png',rgNeale))
    plotNeale.width = image_info(plotNeale)$width
    plotNeale.height = image_info(plotNeale)$height
    plotNeale = ggplot() +
    background_image(plotNeale) + coord_fixed(ratio = image_info(plotNeale)$height/image_info(plotNeale)$width) &
      theme(plot.tag.position  = c(.000, 1.36)) 

    # combine plots
    message('[PNG 3/4] Combining plots.')
    tmp = plot_spacer() + plotSelect + plot_spacer() + plotNeale + plot_layout(widths = c(-400,plotSelect.width,-400, plotNeale.width))
    pl = tmp + plot_annotation(tag_levels = list(c('a','b\n\n\n\n\n\n\n\n\n\n\n\n\n  c'))) &
      theme(plot.tag = element_text(size = 13, face = 'bold', hjust = 0, vjust = 1.9))

    # save file
    message(sprintf('[PNG 4/4] Saving %s.png\n',outFile))
    png(width = 7.79, height = 6.20, units = "in", res = 600, filename = sprintf('%s.png',outFile))
    pl
    invisible(dev.off())

  } else if (fileFormat == ".svg") {
    
    # function to read SVG as grob
    read_svg_grob <- function(file) {
      pic <- suppressWarnings(readPicture(file))
      grid.grabExpr(grid.picture(pic))
    }

    # load rgSelection
    message(sprintf('[SVG 1/4] Loading %s.svg', rgSelection))
    plotSelect = read_svg_grob(sprintf('%s.svg',rgSelection))
    plotSelect = ggplot() +
      annotation_custom(plotSelect, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) &
      theme_void() + 
      theme(plot.tag.position  = c(.3, 0.995),
            plot.margin = margin(5, 0, 0, 0))

    # load rgNeale
    message(sprintf('[SVG 2/4] Loading %s.svg', rgNeale))
    plotNeale = read_svg_grob(sprintf('%s.svg',rgNeale))
    plotNeale = ggplot() +
      annotation_custom(plotNeale, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) &
      theme_void() + 
      theme(plot.tag.position = c(.000, 1.3925),
            plot.margin = margin(0, 0, 0, 0))

    # combine plots
    message('[SVG 3/4] Combining plots.')
    #tmp <- plotSelect + plotNeale
    tmp = plot_spacer() + plotSelect + plot_spacer() + plotNeale + plot_layout(widths = c(-0.3,1,-0.3,1))
    pl <- tmp + plot_annotation(tag_levels = list(c('a', 'b\n\n\n\n\n\n\n\n\n\n\n\n  c'))) &
      theme(plot.tag = element_text(size = 15, face = 'bold', hjust = 0, vjust = 1.9))

    # save file
    message(sprintf('[SVG 4/4] Saving %s.pdf', outFile))
    ggsave(sprintf('%s.pdf',outFile), pl, width = 8.5, height = 6.40)

  }
}

