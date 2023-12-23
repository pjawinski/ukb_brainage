#!/usr/bin/env Rscript

# =========================
# === create GSMR plots ===
# =========================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=6) {
  stop(paste0('expected 6 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}
# set arguments
gsmrFile = args[1] # gsmrFile="results/gap_gm/gsmr/gsmr.eff_plot.gz"
exposure_levels = args[2] # exposure_levels="mr_als_vanRheenen_2021,mr_bip_mullins_2021,mr_bmi_locke_2015,mr_cud_johnson_2020,mr_edu_okbay_2016,mr_eduPruned_okbay_2016,mr_hdlc_willer_2013,mr_height_wood_2014,mr_ldlc_willer_2013,mr_mdd_wray_2018,mr_scz_trubetskoy_2022,mr_trigl_willer_2013,mr_whr_shungin_2015"
exposure_labels = args[3] # exposure_labels="Amyotrophic_Lateral_Sclerosis,Bipolar_Disorder,Body-Mass-Index_(BMI),Cannabis Use Disorder,Educational Attainment,Educational Attainment (pruned),HDL-C,Height,LDL-C,Major Depression,Schizophrenia,Triglyceride,Waist-Hip-Ratio"
outcome_levels = args[4] # outcome_levels="gap_gm"
outcome_labels = args[5] # outcome_labels="Grey matter brain age gap"
targetDir = args[6] # targetDir="results/gap_gm/gsmr/"

message(paste0('\n--- Creating GSMR plots | Settings ---',
               '\ngsmrFile: ', gsmrFile,
               '\nexposure_levels: ', exposure_levels,
               '\nexposure_labels: ', exposure_labels,
               '\noutcome_levels: ', outcome_levels,
               '\noutcome_labels: ', outcome_labels,
               '\ntargetDir: ', targetDir,'\n'))

# attach packages to current R session
for (pkg in c('cowplot', 'dplyr', 'ggplot2','magick','patchwork','stringr')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# source gsmr plot function and load data
message(' - sourcing gsmr_plot.r and loading data.')
wd = funr::get_script_path()
source(sprintf('%s/gsmr_plot.r', wd))
gsmr_data = read_gsmr_data(gsmrFile)

# transform variables
exposure_levels = str_split(exposure_levels, ',')[[1]]
exposure_labels = str_replace_all(exposure_labels, "_", " ")
exposure_labels = str_split(exposure_labels, ',')[[1]]
outcome_levels = str_split(outcome_levels, ',')[[1]]
outcome_labels = str_replace_all(outcome_labels, "_", " ")
outcome_labels = str_split(outcome_labels, ',')[[1]]

# replace levels with labels
message(' - replacing levels with labels.')
for (i in 1:length(exposure_levels)) {
  gsmr_data$pheno[gsmr_data$pheno==exposure_levels[i]] = exposure_labels[i]
  gsmr_data$bxy_result[gsmr_data$bxy_result[,1]==exposure_levels[i],1] = exposure_labels[i]
  gsmr_data$bxy_result[gsmr_data$bxy_result[,2]==exposure_levels[i],2] = exposure_labels[i]
}
for (i in 1:length(outcome_levels)) {
  gsmr_data$pheno[gsmr_data$pheno==outcome_levels[i]] = outcome_labels[i]
  gsmr_data$bxy_result[gsmr_data$bxy_result[,1]==outcome_levels[i],1] = outcome_labels[i]
  gsmr_data$bxy_result[gsmr_data$bxy_result[,2]==outcome_levels[i],2] = outcome_labels[i]
}

# create individual plots
message(' - creating individual plots.')
gsmrSummary = gsmr_summary(gsmr_data)
for (i in 1:length(exposure_levels)) {
  for (j in 1:length(outcome_levels)) {
    
    # get statistics
    idx = gsmrSummary$Exposure == exposure_labels[i] & gsmrSummary$Outcome == outcome_labels[j]
    bxy = sprintf('%.3f', as.numeric(gsmrSummary$bxy[idx]))
    se = sprintf('%.3f', as.numeric(gsmrSummary$se[idx]))
    p = as.numeric(gsmrSummary$p[idx])
    if (bxy == 'NaN') { next }
    
    # create plot and add title
    pdf(file = sprintf('%s/gsmr.exposure.%d.outcome.%d.pdf',targetDir,i,j), width=5.84, height=5.42) # png(filename = 'code/figures/power.png', width=5.98, height=4.48, units = "in", res = 600)
    par(mar = c(4,4,3,2), mgp = c(2.5, 0.7, 0), lwd=1)
    plot_gsmr_effect(gsmr_data, exposure_labels[i], outcome_labels[j], colors()[75])
    title(main = bquote(atop(bold(.(exposure_labels[i])), atop('',''))))
    
    # add statistics
    if (p >= 0.001) { p = sprintf('%.3f', p) } else { p = sprintf('%.1e', p) }
    title(cex.main = 1, main = bquote(atop('', '(' * italic('b')['xy'] ~ '=' ~ .(bxy) * '; s.e. =' ~ .(se) * ';' ~ italic('p')['xy'] ~ '=' ~ .(p) * ')')))
    dev.off()
  }
}

# print 3x4 plots on a page
message('\n - printing 3x4 GSMR plots on a single page.')
plotCount = 0
pageCount = 0
for (i in 1:length(exposure_levels)) {
  for (j in 1:length(outcome_levels)) {
    
    # get plots with available .pdf file
    if (file.exists(sprintf('%s/gsmr.exposure.%d.outcome.%d.pdf',targetDir,i,j))) {  
      plotCount = plotCount + 1
      pl.temp = ggdraw() + 
        draw_image(magick::image_read_pdf(path = sprintf('%s/gsmr.exposure.%d.outcome.%d.pdf',targetDir,i,j), density = 300), scale = 1)
      assign(paste0('pl.',plotCount), pl.temp) 
    }
    
    # add plot spacers if page is incomplete
    if (i*j == length(exposure_levels)*length(outcome_levels) & plotCount < 12) {
      for (k in (plotCount+1):12) {
        pl.temp = plot_spacer()
        assign(paste0('pl.',k), pl.temp)
      }
    }
    
    # draw page
    if (i*j == length(exposure_levels)*length(outcome_levels) | plotCount == 12) {
      plotCount = 0
      pageCount = pageCount + 1
      message(sprintf(' - writing %s/gsmr.exposure.page%s.png',targetDir,pageCount))
      png(file = sprintf('%s/gsmr.exposure.page%s.png',targetDir,pageCount), width=9.6, height=13.3, units = 'in', res = 300)
      print({
        pl.1 + pl.2 + pl.3 + pl.4 + pl.5 + pl.6 + pl.7 + pl.8 + pl.9 + pl.10 + pl.11 + pl.12 + plot_layout(ncol = 3)
      })
      dev.off()    
      system(sprintf('chmod 770 %s/gsmr.exposure.page%s.png',targetDir,pageCount))
    }
  }
}
system(sprintf('rm -f %s/*.pdf',targetDir))
message('-- Completed: Creating GSMR Plots ---\n')

