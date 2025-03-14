#!/usr/bin/env Rscript

# =============================
# === Combine Z score plots ===
# =============================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=9) {
  stop(paste0('expected 9 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
plotTitles = args[1] # plotTitles="Cz_alpha,Cz_beta,Cz_delta,Cz_theta,Cz_broadband,Occ._alpha,Occ._alpha_freq"
plots = args[2] # plots="results/pwr_cz_alpha/metal/compareZ.png,results/pwr_cz_beta/metal/compareZ.png,results/pwr_cz_delta/metal/compareZ.png,results/pwr_cz_theta/metal/compareZ.png,results/pwr_cz_broadband/metal/compareZ.png,results/pwr_occ_alpha/metal/compareZ.png,results/pwr_occ_alphapeakfreq/metal/compareZ.png"
outputFile = args[3] # outputFile="results/combined/compareZ.png"
width = as.numeric(args[4]) # width = 11.6
height = as.numeric(args[5]) # height = 5.6
ncol = as.numeric(args[6]) # ncol = 3
titleSize = as.numeric(args[7]) # titleSize = 10
titlePosX = as.numeric(args[8]) # titlePosX = 0.53
titlePosY = as.numeric(args[9]) # titlePosY = 1

message(paste0('\n--- Combine Z score plots ---',
               '\nplotTitles: ', plotTitles,
               '\nplots: ', plots,
               '\noutputFile: ', outputFile,
               '\nwidth: ', width,
               '\nheight: ', height,
               '\nncol: ', ncol,
               '\ntitleSize: ', titleSize,
               '\ntitlePosX: ', titlePosX,
               '\ntitlePosY: ', titlePosY,'\n'))

# attach required packages
for (pkg in c('ggpubr','ggplot2','patchwork','magick','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# transform variables
plotTitles = str_split(plotTitles, ',')[[1]]
plotTitles = str_replace_all(plotTitles, "_", " ")
plots = str_split(plots, ',')[[1]]

# load images
for (i in 1:length(plotTitles)) {
  
  # load zPlot plots
  message(paste0('[', i, '/', length(plotTitles), '] Loading ', plots[i]))
  zPlot = image_read(plots[i])
  zPlot.width = image_info(zPlot)$width
  zPlot.height = image_info(zPlot)$height
  zPlot = ggplot() +
  background_image(zPlot) + coord_fixed(ratio = image_info(zPlot)$height/image_info(zPlot)$width)
  assign(paste0(plotTitles[i],"_freq"), zPlot)
  
  # collect plots
  if (i == 1) { pl = zPlot } else { pl = pl + zPlot }
}

# set layout
if (length(plotTitles) > 1) {
  pl = pl + plot_layout(ncol = ncol) + plot_annotation(tag_levels = list(c(plotTitles))) &
    theme(plot.tag.position = c(titlePosX, titlePosY),
          plot.tag = element_text(size = titleSize, face = 'bold', vjust = 1))
}

# save file
message('Saving plot.')
png(width = width, height = height, units = "in", res = 300, filename = outputFile); pl; invisible(dev.off())
message('--- Completed: Combine Z score plots ---')
