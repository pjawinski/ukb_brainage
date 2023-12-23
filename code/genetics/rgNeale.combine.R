#!/usr/bin/env Rscript

# ==================================
# === combine rg.results (Neale) ===
# ==================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop(paste0('expected 3 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
traits = args[1] # traits="pwr_cz_alpha,pwr_cz_beta,pwr_cz_delta,pwr_cz_theta,pwr_cz_broadband,pwr_occ_alpha,pwr_occ_alphapeakfreq"
inputFiles = args[2] # inputFiles="results/pwr_cz_alpha/rgNeale/rg.results,results/pwr_cz_beta/rgNeale/rg.results,results/pwr_cz_delta/rgNeale/rg.results,results/pwr_cz_theta/rgNeale/rg.results,results/pwr_cz_broadband/rgNeale/rg.results,results/pwr_occ_alpha/rgNeale/rg.results,results/pwr_occ_alphapeakfreq/rgNeale/rg.results"
outputFile = args[3] # outputFile="results/combined/rgNeale.txt"

message(paste0('\n--- Combine results from genetic correlation analysis (Neale) ---',
               '\ntraits: ', traits,
               '\ninputFiles: ', inputFiles,
               '\noutputFile: ', outputFile,'\n'))

# attach required packages
for (pkg in c('dplyr','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# transform variables
traits = str_split(traits, ',')[[1]]
inputFiles = str_split(inputFiles, ',')[[1]]

# open and merge datasets
for (i in 1:length(traits)) {
  
  # open datasets and keep relevant variables
  message(paste0('Loading ', inputFiles[i]))
  tmp = read.delim(inputFiles[i], sep = '\t', header = T)
  if (i ==1) {
    tmp = tmp[,c('p2','trait_description','category','showcase','h2_obs','h2_obs_se','h2_int','h2_int_se','gcov_int','gcov_int_se','rg','se','z','p')]
  } else { 
    tmp = tmp[,c('p2','rg','se','z','p')]
  }
  
  # calculate FDR and abs(rg)
  tmp$FDR = p.adjust(tmp$p, method = 'BH')
  tmp$rgAbs = abs(tmp$rg)

  # create unique variable names
  names(tmp)[which(names(tmp) %in% c('rg','se','z','p','FDR','rgAbs'))] = 
    paste0(traits[i],c('_rg','_se','_z','_p','_FDR','_rgAbs'))
  
  # merge datasets
  if (i ==1) { df = tmp } else { df = left_join(df,tmp, by = 'p2') }
}
  
# get top p-value and top |rho|
df$top_rgAbs =  df[,grep('_rgAbs',names(df))] %>% apply(1, FUN = max)
df$top_p = df[,grep('_p',names(df))] %>% apply(1, FUN = min)
df$top_FDR = df[,grep('_FDR',names(df))] %>% apply(1, FUN = min)
df = df[order(df$top_FDR,df$top_p,-df$top_rg),]

# get categories
category = str_split_fixed(df$category, ' - ', 4)
tmp = matrix(NA, nrow(category),ncol(category))
for (i in 1:nrow(category)) {
  k = 1
  for (j in ncol(category):1) {
    if (category[i,j] != "") { tmp[i,k] = category[i,j]; k = k + 1 }
  }
} 
category = data.frame(tmp)
names(category) = c('Cat1_Title','Cat2_Title','Cat3_Title','Cat4_Title')
df = cbind(df,category)

# make custom categories in accordance with assigned PHESANT categories
df$custom_category = df$Cat3_Title
df$custom_category[df$Cat2_Title=='Imaging'] = 'Imaging'
df$custom_category[df$Cat2_Title=='Local environment'] = 'Local environment'
df$custom_category[df$Cat2_Title=='Physical activity measurement'] = 'Physical activity'
df$custom_category[df$Cat2_Title=='Cognitive function online'] = 'Cognitive function'
df$custom_category[df$Cat2_Title=='Diet by 24-hour recall'] = 'Diet by 24-hour recall'
df$custom_category[df$Cat2_Title=='Digestive health'] = 'Digestive health'
df$custom_category[df$Cat2_Title=='Mental health'] = 'Mental health'
df$custom_category[df$Cat2_Title=='Work environment'] = 'Work environment'
df$custom_category[df$Cat2_Title=='Baseline characteristics'] = 'Sociodemographics'
df$custom_category[df$Cat2_Title=='Cognitive function'] = 'Cognitive function'
df$custom_category[df$Cat2_Title=='Physical measures'] = 'Physical measures'
df$custom_category[df$Cat3_Title=='Physical activity measurement'] = 'Physical activity'
df$custom_category[df$Cat3_Title=='Baseline characteristics'] = 'Sociodemographics'
df$custom_category[df$Cat3_Title=='Reception'] = 'Sociodemographics'
df$custom_category[df$Cat3_Title=='Employment'] = 'Sociodemographics'
df$custom_category[df$Cat3_Title=='Summary Administration'] = 'Hospital Inpatient - Administration'
df$custom_category[df$Cat3_Title=='Summary Operations'] = 'Operations'
df$custom_category[df$Cat3_Title=='Summary Maternity'] = 'Maternity and sex-specific factors'
df$custom_category[df$Cat3_Title=='Summary Psychiatric'] = 'Mental health'
df$custom_category[df$Cat3_Title=='Sex-specific factors'] = 'Maternity and sex-specific factors'
df$custom_category[df$Cat3_Title=='Health and medical history'] = 'Medical history and conditions'
df$custom_category[df$Cat3_Title=='Medical conditions'] = 'Medical history and conditions'
df$custom_category[df$Cat3_Title=='Asthma outcomes'] = 'Medical history and conditions'
df$custom_category[df$Cat3_Title=='Myocardial infarction outcomes'] = 'Medical history and conditions'
df$custom_category[df$Cat3_Title=='Psychosocial factors'] = 'Mental health'
df$custom_category[df$Cat3_Title=='Family history'] = 'Family history and early life factors'
df$custom_category[df$Cat3_Title=='Early life factors'] = 'Family history and early life factors'

# additional assignments
df$custom_category[df$Cat1_Title=='Biological samples'] = 'Biological samples'
df$custom_category[df$Cat2_Title=='Local environment'] = 'Lifestyle and environment'
df$custom_category[df$Cat3_Title=='Summary Diagnoses'] = 'Medical history and conditions'
df$custom_category[df$Cat1_Title=='N\\A'] = 'Medical history and conditions'
df$custom_category[df$Cat3_Title=='Summary Administration'] = 'Medical history and conditions'
df$custom_category[df$Cat4_Title=='Medication'] = 'Medications'
df$custom_category[df$trait_description=='is_female, based on inferred genetic sex'] = 'Maternity and sex-specific factors'
df$custom_category[df$trait_description=='Number of operations, self-reported'] = 'Medical history and conditions'

# change NAs
df$showcase[df$showcase == 'N/A'] = "NA"
df$category[df$category == 'N\\A'] = "NA"

# create output
message(sprintf('Writing %s',outputFile))
df$`NA` = ""
cols = c('trait_description','custom_category','top_rgAbs','top_p','top_FDR')
for (i in 1:length(traits)) {
  cols = c(cols,'NA',paste0(traits[i],c('_rg','_se','_z','_p','_FDR')))
}
cols = c(cols,'NA','h2_obs','h2_obs_se','h2_int','h2_int_se','gcov_int','gcov_int_se','category','showcase')
output = df[,cols]
write.table(output, outputFile, sep = '\t', quote = F, row.names = F)
message('--- Completed: Combine results from genetic correlation analysis (Neale) ---')

