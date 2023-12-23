#!/usr/bin/env Rscript

# ===============================
# === create correlation plot ===
# ===============================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=8) {
  stop(paste0('expected 8 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
corrTable = args[1] # corrTable="results/combined/rgSelection.txt"
traits = args[2] # traits="gap_gm,gap_wm,gap_gwm"
traitLabels = args[3] # traitLabels="Grey_matter,White_matter,Grey_&_White"
matchVar = args[4] # matchVar="p2"
outFile = args[5] # outFile="results/combined/rgSelection.png"
scaleLimit = as.numeric(args[6]) # scaleLimit=0.35
width = as.numeric(args[7]) # width = 2.76
height = as.numeric(args[8]) # height = 6.85

message(paste0('\n--- Create genetic correlation plot ---',
               '\ncorrTable: ', corrTable,
               '\ncolNames: ', traits,
               '\ntraitLabels: ', traitLabels,
               '\nmatchVar: ', matchVar,
               '\noutFile: ', outFile,
               '\nscaleLimit: ', scaleLimit,
               '\nwidth: ', width,
               '\nheight: ', height,'\n'))

# attach required packages
for (pkg in c('corrplot','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# transform variables
traits = unlist(stringr::str_split(traits, ','))
traitLabels = stringr::str_replace_all(traitLabels,'_',' ')
traitLabels = unlist(stringr::str_split(traitLabels, ','))

# load dataset
message(' - loading %s', corrTable)
df = read.delim(corrTable, header = T)
studyNames = data.frame(matrix(ncol = 2, byrow = T, c(
  '01_mdd_howard_2019', 'Major Depression',
  '01_bip_mullins_2021', 'Bipolar Disorder',
  '01_scz_trubetskoy_2022', 'Schizophrenia',
  '01_adhd_demontis_2023', 'ADHD',
  '01_asd_grove_2019', 'Autism Spectrum Disorder',
  '01_ptsd_nievergelt_2019', 'Post-Traumatic Stress Disorder',
  '01_an_watson_2019', 'Anorexia Nervosa',
  '02_smkInit_saunders_2022', 'Smoking Initiation',
  '02_cpd_saunders_2022', 'Cigarettes Per Day',
  '02_dpw_saunders_2022', 'Drinks per week',
  '02_cof_zhong_2019', 'Coffee consumption',
  '02_cud_johnson_2020', 'Cannabis Use Disorder',
  '03_ad_wightman_2021', 'Alzheimer\'s Disease',
  '03_stroke_mishra_2022', 'Stroke',
  '03_epiGen_ilae_2022', 'Generalized Epilepsy',
  '03_epiFoc_ilae_2022', 'Focal Epilepsy',
  '03_als_vanRheenen_2021', 'Amyotrophic Lateral Sclerosis',
  '04_well_baselmans_2019', 'Well-being',
  '04_neur_baselmans_2019', 'Neuroticism',
  '04_lone_day_2018', 'Loneliness',
  '04_risk_karlssonlinner_2019', 'Risk-taking tendency',
  '05_ins_jansen_2019', 'Insomnia',
  '05_sd_dashti_2019', 'Sleep Duration',
  '05_chron_jones_2019', 'Chronotype',
  '05_ds_wang_2019', 'Daytime sleepiness',
  '06_ea_okbay_2022', 'Educational Attainment',
  '06_int_savage_2018', 'Intelligence',
  '06_cp_lee_2018', 'Cognitive Performance',
  '06_rt_davies_2018', 'Reaction Time',
  '06_mem_davies_2016', 'Memory performance',
  '07_height_yengo_2022', 'Height',
  '07_bmi_yengo_2018', 'Body-Mass-Index',
  '07_whr_pulit_2019', 'Waist-Hip Ratio',
  '08_cad_aragam_2022', 'Coronary Artery Disease',
  '08_dbp_evangelou_2018', 'Diastolic Blood Pressure',
  '08_sbp_evangelou_2018', 'Systolic Blood Pressure',
  '08_myocardial_hartiala_2021', 'Myocardial Infarction',
  '08_diabetes_xue_2018', 'Type-2 Diabetes')))
names(studyNames) = c(matchVar,'trait')

# get rg, p, and FDR columns
df = dplyr::left_join(studyNames, df, matchVar)
rg = as.matrix(df[,paste0(traits,'_rg')])
p = as.matrix(df[,paste0(traits,'_p')])
fdr = as.matrix(df[,paste0(traits,'_FDR')])
rownames(rg) = rownames(p) = rownames(fdr) = df$trait
colnames(rg) = colnames(p) = colnames(fdr) = traitLabels
  
# set non-significant correlations to 0 (deprecated)
# set nominally significant p-values to 0.049  (for plotting)
# set FDR-significant p-values to 0.0009 (for plotting)
rg[p > 0.05] = 0
p[p < 0.05] = 0.049
p[fdr < 0.05] = 0.0009

# set background color
tlCol = c(rep('#183E60',7), rep('#900000',5), rep('#183E60',5), rep('#900000',4), rep('#183E60',4), rep('#900000',5), rep('#183E60',3), rep('#900000',5))

# draw plot
message(sprintf(' - saving %s',outFile))
png(filename = outFile, width = width, height = height, units = 'in', res = 600)
par(xpd = TRUE)
corrplot(rg, method = 'color', addgrid.col="grey95",
     mar = c(2, 0.1, 0.1, 0.1),
     is.corr = FALSE, #, bg = backColor,
     p.mat = p, sig.level = c(0.001,0.05), insig = 'label_sig', pch.cex = 0.6,  # number.cex = 0.4, addCoef.col ='black',
     tl.cex = 0.6,
     tl.col = tlCol,
     tl.srt = 45,
     col.lim = c(-scaleLimit,scaleLimit),
     col=rev(COL2('RdBu', 200)),
     cl.pos='n')

    # add category label background
    xleft = length(traits)+0.5; xright = length(traits)+2
    rect(xleft = xleft, xright = xright, ybottom = 31.5, ytop = 38.5, col = '#EDF7F9', border = FALSE)
    rect(xleft = xleft, xright = xright, ybottom = 26.5, ytop = 31.5, col = '#FDEFEE', border = FALSE)
    rect(xleft = xleft, xright = xright, ybottom = 21.5, ytop = 26.5, col = '#EDF7F9', border = FALSE)
    rect(xleft = xleft, xright = xright, ybottom = 17.5, ytop = 21.5, col = '#FDEFEE', border = FALSE)
    rect(xleft = xleft, xright = xright, ybottom = 13.5, ytop = 17.5, col = '#EDF7F9', border = FALSE)
    rect(xleft = xleft, xright = xright, ybottom = 8.5, ytop = 13.5, col = '#FDEFEE', border = FALSE)
    rect(xleft = xleft, xright = xright, ybottom = 5.5, ytop = 8.5, col = '#EDF7F9', border = FALSE)
    rect(xleft = xleft, xright = xright, ybottom = 0.5, ytop = 5.5, col = '#FDEFEE', border = FALSE)
    
    # add category labels
    xleft = length(traits)+1.25;
    text(xleft, 35, 'Psychiatric', cex = 0.6, srt = 90, col = '#183E60')
    text(xleft, 29, 'Substance use', cex = 0.6, srt = 90, col = '#900000')
    text(xleft, 24, 'Neurological', cex = 0.6, srt = 90, col = '#183E60')
    text(xleft, 19.5, 'Personality', cex = 0.6, srt = 90, col = '#900000')
    text(xleft, 15.5, 'Sleep', cex = 0.6, srt = 90, col = '#183E60')
    text(xleft, 11, 'Cognition', cex = 0.6, srt = 90, col = '#900000')
    text(xleft, 7, 'Anthrop.', cex = 0.6, srt = 90, col = '#183E60')
    text(xleft, 3, 'Cardiovascular', cex = 0.6, srt = 90, col = '#183E60')
    
    # add color legend and labels
    colorlegend(colbar = rev(COL2('RdBu', 200)), cex = 0.5, xlim=c(-8,length(traits)+1), ylim=c(-1.75,-0.5),  labels = c('',''), vertical=FALSE)
    text(-7.75, -0.8, -scaleLimit, pos = 2, cex = 0.6, col = 'black')
    text(length(traits)+0.75, -0.8, scaleLimit, pos = 4, cex = 0.6, col = 'black')
    text(mean(c(-7.75,length(traits)+1)), -0.8, expression("Genetic correlation (r"["G"]*')'), pos = 1, cex = 0.6, col = 'black')
    
    # add outer lines for categories
    xleft = -9; xright = length(traits)+2
    lines(x = c(xleft,xright), y = c(38.5,38.5), lwd = 1, lty = 3, col = 'grey30')
    lines(x = c(xleft,xright), y = c(31.5,31.5), lwd = 1, lty = 3, col = 'grey30')
    lines(x = c(xleft,xright), y = c(26.5,26.5), lwd = 1, lty = 3, col = 'grey30')
    lines(x = c(xleft,xright), y = c(21.5,21.5), lwd = 1, lty = 3, col = 'grey30')
    lines(x = c(xleft,xright), y = c(17.5,17.5), lwd = 1, lty = 3, col = 'grey30')
    lines(x = c(xleft,xright), y = c(13.5,13.5), lwd = 1, lty = 3, col = 'grey30')
    lines(x = c(xleft,xright), y = c(8.5,8.5), lwd = 1, lty = 3, col = 'grey30')
    lines(x = c(xleft,xright), y = c(5.5,5.5), lwd = 1, lty = 3, col = 'grey30')
    lines(x = c(xleft,xright), y = c(0.5,0.5), lwd = 1, lty = 3, col = 'grey30')
    lines(x = c(0.5,0.5), y = c(0.5,38.5), lwd = 1, lty = 3, col = 'grey30')
    lines(x = c(xright-1.5,xright-1.5), y = c(0.5,38.5), lwd = 1, lty = 3, col = 'grey30')
invisible(dev.off())
message('--- Completed: Create genetic correlation plot ---')
