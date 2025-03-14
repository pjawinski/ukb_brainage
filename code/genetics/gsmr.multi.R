#!/usr/bin/env Rscript

# ====================================================================================================
# === run multiple alternative mendelian randomization analyses on variants extracted by GCTA-GSMR ===
# ====================================================================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop(paste0('expected 3 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
gsmrFile = args[1] # gsmrFile="results/gap_wm/gwama/eur/gsmr/gsmr.eff_plot.gz"
gsmrFileNoHEIDI = args[2] # gsmrFileNoHEIDI="results/gap_wm/gwama/eur/gsmr/gsmr.noHEIDI.eff_plot.gz"
outFile = args[3] # outFile="results/gap_wm/gwama/eur/gsmr/gsmr.mrmulti"

message(paste0('\n--- Run multiple MR analyses | Settings ---',
               '\ngsmrFile: ', gsmrFile,
               '\ngsmrFileNoHEIDI: ', gsmrFileNoHEIDI,
               '\noutFile: ', outFile,'\n'))

# attach packages to current R session
for (pkg in c('doParallel','dplyr','foreach','MendelianRandomization','stringr')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# source gsmr plot function and load data
message(' - sourcing gsmr_plot.r and loading data.')
wd = funr::get_script_path()
source(sprintf('%s/gsmr_plot.r', wd))
gsmr_data = read_gsmr_data(gsmrFile)
gsmrSummary = gsmr_summary(gsmr_data)
gsmrNoHEIDI_data = read_gsmr_data(gsmrFileNoHEIDI)
gsmrSummaryNoHEIDI = gsmr_summary(gsmrNoHEIDI_data)

# run mr analyses
message('\n - running multiple MR analyses.')
cl = parallel::makeForkCluster(nrow(gsmrSummaryNoHEIDI))
doParallel::registerDoParallel(cl)
mrmany = foreach(i=1:nrow(gsmrSummaryNoHEIDI), .combine = 'rbind') %dopar% {

  # create output data
  tmp = matrix(NA, ncol = 37, 1) %>% data.frame()
    names(tmp) = c('exposure','outcome', 'n_snps_total', 'n_snps_HEIDI',
      'gsmr_beta', 'gsmr_se', 'gsmr_p',
      'ivw_beta','ivw_se','ivw_p',
      'median_beta','median_se','median_p',
      'egger_beta','egger_se','egger_p','egger_intbeta','egger_intse','egger_intp',
      'ml_beta','ml_se','ml_p',
      'mbe_beta','mbe_se','mbe_p',
      'cm_beta','cm_se','cm_p',
      'divw_beta','divw_se','divw_p',
      'pivw_beta','pivw_se','pivw_p',
      'lasso_beta','lasso_se','lasso_p')

    # get gsmr statistics
    tmp[1,1:2] = c(gsmrSummaryNoHEIDI$Exposure[i], gsmrSummaryNoHEIDI$Outcome[i])
    if (gsmrSummary$bxy[i] != 'NaN') { 
      bxy = as.numeric(gsmrSummary$bxy[i])
      se = as.numeric(gsmrSummary$se[i])
      p = as.numeric(gsmrSummary$p[i])
      n_snps_HEIDI = as.numeric(gsmrSummary$n_snps[i])
      n_snps_total = as.numeric(gsmrSummaryNoHEIDI$n_snps[i])

      # get snp effects
      snpeffect = gsmr_snp_effect(gsmrNoHEIDI_data, gsmrSummaryNoHEIDI$Exposure[i], gsmrSummaryNoHEIDI$Outcome[i]) %>% as.data.frame()
      MRInputObject = mr_input(
        snps = snpeffect$snp,
        bx = snpeffect$bzx,
        bxse = snpeffect$bzx_se,
        by = snpeffect$bzy,
        byse = snpeffect$bzy_se)

      # run MR analyses
      ivw = mr_ivw(MRInputObject, model = "default", robust = FALSE, penalized = FALSE, correl = FALSE, weights = "simple", psi = 0, distribution = "normal", alpha = 0.05)
      median = mr_median(MRInputObject, weighting = "weighted", distribution = "normal", alpha = 0.05, iterations = 10000, seed = 314159265)
      egger = mr_egger(MRInputObject, robust = FALSE, penalized = FALSE, correl = FALSE, distribution = "normal", alpha = 0.05)
      ml = mr_maxlik(MRInputObject, model = "default", correl = FALSE, psi = 0, distribution = "normal", alpha = 0.05)
      mbe = mr_mbe(MRInputObject, weighting = "weighted", stderror = "delta", phi = 1, seed = 314159265, iterations = 10000, distribution = "normal", alpha = 0.05)
      cm = mr_conmix(MRInputObject, psi = 0, CIMin = NA, CIMax = NA, CIStep = 0.01, alpha = 0.05)
      divw = mr_divw(MRInputObject, over.dispersion = TRUE, alpha = 0.05, diagnostics = FALSE)
      pivw = mr_pivw(MRInputObject, over.dispersion = TRUE, delta = 0, sel.pval = NULL, Boot.Fieller = TRUE, alpha = 0.05)
      ltry = try(mr_lasso(MRInputObject, distribution = "normal", alpha = 0.05, lambda = numeric(0)), silent = T)
      if (!inherits(ltry, "try-error")) { lasso.Estimate = ltry$Estimate; lasso.StdError = ltry$StdError; lasso.Pvalue = ltry$Pvalue } else { lasso.Estimate = 0; lasso.StdError = 1; lasso.Pvalue = 1 }

      # calculate se for contamination mixture model (only if confidence interval contains single range of values)
      if(length(cm$CIUpper) > 1) { cm.StdError = NA } else { cm.StdError = (cm$CIUpper-cm$CILower)/(2*1.959964) }

      # paste results into data frame
      tmp[1,3:ncol(tmp)] = c(
          n_snps_total, n_snps_HEIDI, bxy, se, p,
          ivw$Estimate, ivw$StdError, ivw$Pvalue, 
          median$Estimate, median$StdError, median$Pvalue, 
          egger$Estimate, egger$StdError.Est, egger$Pvalue.Est, egger$Intercept, egger$StdError.Int, egger$Pvalue.Int, 
          ml$Estimate, ml$StdError, ml$Pvalue, 
          mbe$Estimate, mbe$StdError, mbe$Pvalue, 
          cm$Estimate, cm.StdError, cm$Pvalue, 
          divw$Estimate, divw$StdError, divw$Pvalue, 
          pivw$Estimate, pivw$StdError, pivw$Pvalue,
          lasso.Estimate, lasso.StdError, lasso.Pvalue)
    } else {
      tmp[1,3:ncol(tmp)] = NA
    }
    return(tmp)
}
parallel::stopCluster(cl)

# calculate average (median) p
mrmany$avg_p = apply(mrmany[,c('gsmr_p','ivw_p','median_p','egger_p','ml_p','mbe_p','cm_p','divw_p','pivw_p','lasso_p')], 1, median, na.rm=TRUE)
mrmany$sum_p05 = apply(mrmany[,c('gsmr_p','ivw_p','median_p','egger_p','ml_p','mbe_p','cm_p','divw_p','pivw_p','lasso_p')]<0.05, 1, sum, na.rm=TRUE)

# write outfile
message(sprintf(' - saving %s', outFile))
write.table(mrmany, outFile, sep = '\t', col.names = T, row.names = F, quote = F)
system(sprintf('chmod 770 %s',outFile))
message('-- Completed: Run multiple MR analyses ---\n')

