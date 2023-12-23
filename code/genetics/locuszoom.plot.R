#!/usr/bin/env Rscript

# ========================================
# === combine locuszoom regional plots ===
# ========================================

# get arguments from command line
args = commandArgs(trailingOnly=TRUE)

# set arguments
traitDescription = args[1] # traitDescription = traitDescription="grey matter"
targetDir = args[2] # targetDir = "results/gap_gm/locuszoom"
conditionalFile = args[3] # conditionalFile = "results/gap_gm/conditional/conditional.cleaned.tophits.annovar.txt"

message(paste0('\n--- Combining Locuszoom Plots | Settings ---',
        '\ntraitDescription: ', traitDescription,
        '\ntargetDir: ', targetDir,
        '\nconditionalFile: ', conditionalFile,'\n'))

# attach packages to current R session
for (pkg in c('cowplot', 'dplyr', 'ggplot2', 'magick', 'patchwork')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# get pdf file list and lead snps data (A1_FREQ, INFO, ETA2)
message('Retrieving pdf file list and lead snp data.')
df.cond = read.delim(conditionalFile, sep = '\t', header = T, quote = "")
pdflist = list.files(paste0(targetDir, '/lz.output/'))
snplist = pdflist %>% gsub(pattern='^.+?_(.*)', replace='\\1') %>% sub(pattern='_', replace=':') %>% gsub(pattern='.pdf',replacement='')
df.lz = data.frame(pdf = pdflist, ID = snplist) %>% left_join(y = df.cond, by = 'ID')

# sort by CHR and BP
message('Sorting by chrommosome and base pair.')
df.lz = df.lz %>% mutate(CHR_num = replace(CHR, CHR=='X', 23)) %>%
		  mutate(CHR_num = replace(CHR_num, CHR_num=='Y', 24)) %>%
		  mutate(CHR_num = replace(CHR_num, CHR_num=='XY', 25)) %>%
		  mutate(CHR_num = replace(CHR_num, CHR_num=='MT', 26)) %>%
		  mutate(CHR_num = as.numeric(CHR_num)) %>%
		  arrange(CHR_num, BP) %>%
		  select(-CHR_num)

# replace A1_FREQ by MAF
message('Replacing A1_FREQ by MAF.')
df.lz$MAF[df.lz$MAF > 0.5] = 1 - df.lz$MAF[df.lz$MAF > 0.5] 

# print 2x3 plots on a page
message('Printing 2x3 Locuszoom plots on a single page.')
for (i in seq(1,nrow(df.lz),6)) {
    for (j in 0:5) {
      if (i+j <= nrow(df.lz)) {  
         lz.temp = ggdraw() + 
              draw_image(magick::image_read_pdf(path = paste0(targetDir, '/lz.output/', df.lz$pdf[i+j]), density = 300), scale = 1) +
              draw_label(paste0(traitDescription, '\n', df.lz$ID[i+j],
              	' (MAF = ',  format(round(df.lz$A1_FREQ[i+j], 2), nsmall = 2), 
              	', INFO = ', format(round(df.lz$INFO[i+j], 2), nsmall = 2),
              	', Gene: ', df.lz$NEAREST_GENE[i+j], ')'),
              	 y = 0.97, size = 8, lineheight = 1.1)
      } else {
         lz.temp = plot_spacer()
      }
      assign(paste0('lz.',j), lz.temp) 
    }

    png(file = paste0(targetDir, '/lz.', (i+5)/6, '.png'), width=9.6, height=13.3, units = 'in', res = 300)
	print({
		lz.0 + lz.1 + lz.2 + lz.3 + lz.4 + lz.5 + plot_layout(ncol = 2)
       })
    dev.off()    
    system(paste0('chmod 770 ', targetDir, '/lz.', (i+5)/6, '.png'))
}
message('-- Combining Locuszoom Plots | Finished ---\n')


