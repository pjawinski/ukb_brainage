#!/usr/bin/env Rscript

# ===================================================================
# === merge single-model plots (benchmark and feature importance) ===
# ===================================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop(paste0('expected 4 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
benchmark = args[1] # benchmark="results/mri/singlemodel.benchmark.png"
gmImportance = args[2] # gmImportance="results/mri/ml.xgb/singlemodel/xgb_gm_importance.png"
wmImportance = args[3] # wmImportance="results/mri/ml.xgb/singlemodel/xgb_wm_importance.png"
outFile = args[4] # outFile="results/mri/singlemodel.composite.png"

message(paste0('\n--- composite single-model plots (benchmark and feature importance) | settings ---',
               '\nbenchmark: ', benchmark,
               '\ngmImportance: ', gmImportance,
               '\nwmImportance: ', wmImportance,
               '\noutFile: ', outFile,'\n'))

# attach required packages
for (pkg in c('dplyr','magick','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# load images
message(' - loading images')
p1 = image_read(benchmark) %>% image_trim()
p2 = image_read(gmImportance) %>% image_trim() %>% image_scale("x655")
p3 = image_read(wmImportance) %>% image_trim() %>% image_scale("x655")

# Composite the images
message(' - crating composite')
canvas = image_blank(width = 3375, height = 2250, color = "white")
composite_img = canvas %>%
  image_composite(p1, offset = "+50+50") %>%
  image_composite(p2, offset = "+1500+50") %>%
  image_composite(p3, offset = "+1500+790")

# save image
message(sprintf(' - saving %s',outFile))
image_write(composite_img, outFile)
message('--- Completed: composite single-model plots (benchmark and feature importance) ---')
