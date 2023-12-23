#!/usr/bin/env Rscript

# ==========================================
# === replace rg trait names with labels ===
# ==========================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop(paste0('expected 3 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
inputFile = args[1] # inputFile="results/combined/rgSelection.txt"
outFile = args[2] # outFile="results/combined/rgSelection.labels.txt"
matchVar = args[3] # matchVar="p2"

message(paste0('\n--- Replace rg trait names with labels ---',
               '\ninputFile: ', inputFile,
               '\noutFile: ', outFile,
               '\nmatchVar: ', matchVar,'\n'))

# load dataset
message(sprintf(' - loading %s', inputFile))
df = read.delim(inputFile, header = T)
studyNames = data.frame(matrix(ncol = 2, byrow = T, c(
  '01_mdd_howard_2019', 'Major Depression (Howard et al. 2019)',
  '01_bip_mullins_2021', 'Bipolar Disorder (Mullins et al. 2021)',
  '01_scz_trubetskoy_2022', 'Schizophrenia (Trubetskoy et al. 2022)',
  '01_adhd_demontis_2023', 'ADHD (Demontis et al. 2023)',
  '01_asd_grove_2019', 'Autism Spectrum Disorder (Grove et al. 2019)',
  '01_ptsd_nievergelt_2019', 'Post-Traumatic Stress Disorder (Nievergelt et al. 2019)',
  '01_an_watson_2019', 'Anorexia Nervosa (Watson et al. 2019)',
  '02_smkInit_saunders_2022', 'Smoking Initiation (Saunders et al. 2022)',
  '02_cpd_saunders_2022', 'Cigarettes Per Day (Saunders et al. 2022)',
  '02_dpw_saunders_2022', 'Drinks per week (Saunders et al. 2022)',
  '02_cof_zhong_2019', 'Coffee consumption (Zhong et al. 2019)',
  '02_cud_johnson_2020', 'Cannabis Use Disorder (Johnson et al. 2020)',
  '03_ad_wightman_2021', 'Alzheimer\'s Disease (Wightman et al. 2021)',
  '03_stroke_mishra_2022', 'Stroke (Mishra et al. 2022)',
  '03_epiGen_ilae_2022', 'Generalized Epilepsy (ILAE consortium 2022)',
  '03_epiFoc_ilae_2022', 'Focal Epilepsy (ILAE consortium 2022)',
  '03_als_vanRheenen_2021', 'Amyotrophic Lateral Sclerosis (van Rheenen et al. 2021)',
  '04_well_baselmans_2019', 'Well-being (Baselmans et al. 2019)',
  '04_neur_baselmans_2019', 'Neuroticism (Baselmans et al. 2019)',
  '04_lone_day_2018', 'Loneliness (Day et al. 2018)',
  '04_risk_karlssonlinner_2019', 'Risk-taking tendency (Karlsson Linner et al. 2019)',
  '05_ins_jansen_2019', 'Insomnia (Jansen et al. 2019)',
  '05_sd_dashti_2019', 'Sleep Duration (Dashti et al. 2019)',
  '05_chron_jones_2019', 'Chronotype (Jones et al. 2019)',
  '05_ds_wang_2019', 'Daytime sleepiness (Wang et al. 2019)',
  '06_ea_okbay_2022', 'Educational Attainment (Okbay et al. 2022)',
  '06_int_savage_2018', 'Intelligence (Savage et al. 2018)',
  '06_cp_lee_2018', 'Cognitive Performance (Lee et al. 2018)',
  '06_rt_davies_2018', 'Reaction Time (Davies et al. 2018)',
  '06_mem_davies_2016', 'Memory performance (Davies et al. 2016)',
  '07_height_yengo_2022', 'Height (Yengo et al. 2022)',
  '07_bmi_yengo_2018', 'Body-Mass-Index (Yengo et al. 2018)',
  '07_whr_pulit_2019', 'Waist-Hip Ratio (Pulit et al. 2019)',
  '08_cad_aragam_2022', 'Coronary Artery Disease (Aragam et al. 2022)',
  '08_dbp_evangelou_2018', 'Diastolic Blood Pressure (Evangelou et al. 2018)',
  '08_sbp_evangelou_2018', 'Systolic Blood Pressure (Evangelou et al. 2018)',
  '08_myocardial_hartiala_2021', 'Myocardial Infarction (Hartiala et al. 2021)',
  '08_diabetes_xue_2018', 'Type-2 Diabetes (Xue et al. 2018)')))
names(studyNames) = c(matchVar,'label')

# add labels and drop matchVar
df = dplyr::left_join(df, studyNames, matchVar)
df = df[,-which(matchVar %in% names(df))]
df = df[,c('label',names(df)[!(names(df) %in% 'label')])]

# write output
message(sprintf(' - writing %s',outFile))
write.table(df, paste0(outFile), sep = '\t', quote = F, row.names = F, na = '')
message('--- Completed: Combine results from genetic correlation analysis ---')

