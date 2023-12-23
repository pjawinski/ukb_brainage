#!/usr/bin/env Rscript

# ======================================
# === create combined manhattan plot ===
# ======================================

# set working directory
setwd('/home/groups/markett/ukb_brainage')

# attach required packages
library(ggpubr)
library(ggplot2)
library(patchwork)
library(magick)

# load manhattan plots
for (trait in c('gm', 'wm', 'gwm')) {
  tmp = image_read(paste0('results/gap_',trait, '/manhattan/manhattan.png'))
  tmp = ggplot() +
        background_image(tmp) + coord_fixed(ratio = image_info(tmp)$height/image_info(tmp)$width)
  assign(trait, tmp)
}

# merge them
plot = gm / wm / gwm  + plot_annotation(tag_levels ='a') & 
  theme(plot.tag.position = c(0, 1), plot.tag = element_text(size = 22))

# save file
png(width = 11, height = 4.0*3, units = "in", res = 300, filename = 'results/combined/snplevel.manhattan.png')
plot
invisible(dev.off())
