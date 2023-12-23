#!/usr/bin/env Rscript

# ============================================
# === replace gsmr trait names with labels ===
# ============================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop(paste0('expected 3 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
inputFile = args[1] # inputFile="results/combined/gsmr.txt"
outFile = args[2] # outFile="results/combined/gsmr.labels.txt"
matchVar = args[3] # matchVar="trait"

message(paste0('\n--- Replace gsmr trait names with labels ---',
               '\ninputFile: ', inputFile,
               '\noutFile: ', outFile,
               '\nmatchVar: ', matchVar,'\n'))

# load dataset
message(sprintf(' - loading %s', inputFile))
df = read.delim(inputFile, header = T)
studyNames = data.frame(matrix(ncol = 2, byrow = T, c(
  'mr_bmi_locke_2015','Body-Mass-Index (Locke et al. 2015)',
  'mr_whr_shungin_2015','Waist-Hip-Ratio (Shungin et al. 2015)',
  'mr_dbp_noUKB_evangelou_2018','Diastolic Blood Pressure (Evangelou et al. 2018)',
  'mr_sbp_noUKB_evangelou_2018','Systolic Blood Pressure (Evangelou et al. 2018)',
  'mr_pp_noUKB_evangelou_2018','Pulse Pressure (Evangelou et al. 2018)',
  'mr_ldlc_willer_2013','LDL-c (Willer et al. 2013)',
  'mr_hdlc_willer_2013','HDL-c (Willer et al. 2013)',
  'mr_trigl_willer_2013','Triglyceride (Willer et al. 2013)',
  'mr_t2d_scott_2017','Type-2-Diabetes (Scptt et al. 2017)',
  'mr_cad_nikpay_2015','Coronary Artery Disease (Nikpay et al. 2015)',
  'mr_scz_trubetskoy_2022','Schizophrenia (Trubetskoy et al. 2022)',
  'mr_eduPruned_okbay_2016','Educational attainment (Okbay et al. 2016)')))
  names(studyNames) = c(matchVar,'label')

# add labels and drop matchVar
df = dplyr::left_join(df, studyNames, matchVar)
df = df[,-which(matchVar %in% names(df))]
df = df[,c('label',names(df)[!(names(df) %in% 'label')])]

# write output
message(sprintf(' - writing %s',outFile))
write.table(df, paste0(outFile), sep = '\t', quote = F, row.names = F, na = '')
message('--- Completed: Replace gsmr trait names with labels ---')

