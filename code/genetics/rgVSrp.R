#!/usr/bin/env Rscript

# =======================================================
# === compare genetic against phenotypic correlations ===
# =======================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=7) {
  stop(paste0('expected 7 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
traits = args[1] # traits = "gap_gm,gap_wm,gap_gwm"
rgFile = args[2] # rgFile="results/combined/gwama.eur.rgNeale.txt"
phenoFile = args[3] # phenoFile="results/combined/discovery.phewas.txt"
outFile = args[4] # outFile="results/combined/rgXrp.png"
width = as.numeric(args[5]) # width=7.77
height = as.numeric(args[6]) # height=2.83
ncols = as.numeric(args[7]) # ncols=3

message(paste0('\n--- compare gentic against phenotypic correlations | settings ---',
               '\nrgFile: ', rgFile,
               '\nphenoFile: ', phenoFile,
               '\noutFile: ', outFile,
               '\nwidth: ', width,
               '\nheight: ', height,
               '\nncols: ', ncols,'\n'))

# attach required packages
for (pkg in c('ggplot2','patchwork','dplyr','purrr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# transform variables
traits = stringr::str_split(traits, ',')[[1]]

# load datasets and join them
message(' - loading datasets')
rg = read.delim(rgFile, sep = '\t', header = T)
ph = read.delim(phenoFile, sep = '\t', header = T)
df = inner_join(rg,ph,by = c('trait_description' = 'description'))

# join data frames
# do not include results from multinomial-logistoc (no directions)
# only select traits with significant heritability
df = df[df$resType!='MULTINOMIAL-LOGISTIC' & df$h2_obs/df$h2_obs_se > 1.96,]

# function for correlation plot
corrplot <- function(data, x_var, y_var) {
  
  # Compute correlation coefficient and p-value
  cor_test = cor.test(data[[x_var]], data[[y_var]])
  correlation = round(cor_test$estimate, 2)
  
  # Format p-value based on size
  if (cor_test$p.value < 0.001) {
    p_value = sprintf("p = %s", formatC(cor_test$p.value, format = "e", digits = 1))
  } else {
    p_value = sprintf("p = %s", format.pval(cor_test$p.value, digits = 2))
  }
  
  # Mean absolute difference calculation
  mean_abs_diff = round(mean(abs(data[[x_var]] - data[[y_var]]), na.rm = TRUE), 3)
  

  # Generate the plot
  ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_point(color = "#4F71BE", alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", color = "black", se = TRUE, size = 0.5) +
    labs(
      x = "Genetic correlation",
      y = "Phenotypic correlation",
      title = paste0("r = ", correlation, " | ", p_value, " | MAD = ", mean_abs_diff)
      ) +
    coord_cartesian(xlim = c(-0.25, 0.25), ylim = c(-0.1, 0.1)) +
    theme_bw() +
    theme(plot.title = element_text(size = 8, hjust = 0.5, vjust = 0))
}

# get plots
plots = list()
for (i in 1:length(traits)) {
  assign(traits[i], corrplot(df, paste0(traits[i],"_rg"), paste0(traits[i],"_rho")))
  plots[[i]] <- get(traits[i])
}

# merge them
pl = reduce(plots, `+`) +
  plot_annotation(tag_levels = "a") &  
  plot_layout(ncol = ncols) +
  theme(plot.tag = element_text(size = 16, face = "bold"))

# save file
message(sprintf(' - saving %s',outFile))
png(width = width, height = height, units = "in", res = 300, filename = outFile)
pl
invisible(dev.off())
message('--- Completed: compare gentic against phenotypic correlations ---')
