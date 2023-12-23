#!/usr/bin/env Rscript

# ======================================
# === combine manhattan and qq-plots ===
# ======================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=7) {
  stop(paste0('expected 7 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
traits = args[1] # traits="pwr_cz_alpha,pwr_cz_beta,pwr_cz_delta,pwr_cz_theta,pwr_cz_broadband,pwr_occ_alpha,pwr_occ_alphapeakfreq"
plotTitles = args[2] # plotTitles="Cz_alpha,Cz_beta,Cz_delta,Cz_theta,Cz_broadband,Occ._alpha,Occ._alpha_peak_frequency"
manhattanPlots = args[3] # manhattanPlots="results/pwr_cz_alpha/manhattan/manhattan.png,results/pwr_cz_beta/manhattan/manhattan.png,results/pwr_cz_delta/manhattan/manhattan.png,results/pwr_cz_theta/manhattan/manhattan.png,results/pwr_cz_broadband/manhattan/manhattan.png,results/pwr_occ_alpha/manhattan/manhattan.png,results/pwr_occ_alphapeakfreq/manhattan/manhattan.png"
qqPlots = args[4] # qqPlots="results/pwr_cz_alpha/qqplot/qqplot.png,results/pwr_cz_beta/qqplot/qqplot.png,results/pwr_cz_delta/qqplot/qqplot.png,results/pwr_cz_theta/qqplot/qqplot.png,results/pwr_cz_broadband/qqplot/qqplot.png,results/pwr_occ_alpha/qqplot/qqplot.png,results/pwr_occ_alphapeakfreq/qqplot/qqplot.png"
outputFile = args[5] # outputFile="results/combined/qqplot.manhattan.png"
width = as.numeric(args[6]) # width = 14
height = as.numeric(args[7]) # height = 17.5

message(paste0('\n--- Combine Manhattan and qq-plots ---',
               '\ntraits: ', traits,
               '\nplotTitles: ', plotTitles,
               '\nmanhattanPlots: ', manhattanPlots,
               '\nqqPlots: ', qqPlots,
               '\noutputFile: ', outputFile,
               '\nwidth: ', width,
               '\nheight: ', height,'\n'))

# attach required packages
for (pkg in c('ggpubr','ggplot2','patchwork','magick','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# transform variables
traits = str_split(traits, ',')[[1]]
plotTitles = str_split(plotTitles, ',')[[1]]
plotTitles = str_replace_all(plotTitles, "_", " ")
manhattanPlots = str_split(manhattanPlots, ',')[[1]]
qqPlots = str_split(qqPlots, ',')[[1]]

# load images
for (i in 1:length(traits)) {
  
  # load manhattan plots
  message(paste0('[', i, '/', length(traits), '] Loading ', manhattanPlots[i], ' and ', qqPlots[i]))
  manh = image_read(manhattanPlots[i])
  manh.width = image_info(manh)$width
  manh.height = image_info(manh)$height
  manh = ggplot() +
  background_image(manh) + coord_fixed(ratio = image_info(manh)$height/image_info(manh)$width)
  assign(paste0(traits[i],"_manh"), manh)
  
  # load qq-plots
  qq = image_read(qqPlots[i])
  qq.width = image_info(qq)$width
  qq.height = image_info(qq)$height
  qq = ggplot() +
  background_image(qq) + coord_fixed(ratio = image_info(qq)$height/image_info(qq)$width)
  assign(paste0(traits[i],"_qq"), qq)
  
  # collect plots
  tmp = manh + qq + plot_layout(widths = c(manh.width,qq.width))
  if (i == 1) { pl = tmp } else { pl = pl / tmp }
}

# set layout
if (length(plotTitles) == 1 & length(traits) > 1) {
  pl = pl + plot_annotation(tag_levels = list(c(rbind(letters[1:length(traits)],letters[(1+length(traits)):(2*length(traits))])))) & # list(c(rbind(letters,' ')))
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(size = 22))
  } else if (length(traits) > 1) {
  pl = pl + plot_annotation(tag_levels = list(c(rbind(plotTitles,' ')))) &
    theme(plot.tag.position = c(0.52, 1),
          plot.tag = element_text(size = 16, face = 'bold', hjust = 0.5, vjust = 1.9))
}

# save file
message(sprintf(' - saving %s',outputFile))
png(width = width, height = height, units = "in", res = 300, filename = outputFile)
pl
invisible(dev.off())
message('--- Completed: Combine Manhattan and qq-plots ---')
