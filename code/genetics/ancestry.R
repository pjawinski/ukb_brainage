#!/usr/bin/env Rscript

# =======================================================
# ===  plot ancestry components vs. trait of interest ===
# =======================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=13) {
  stop(paste0('expected 13 arguments, but ', length(args), ' arguments provided.'), call.=FALSE)
}

# get arguments from command line
outFile = args[1] # outFile="results/gap_gm/discovery/pcs/pcs.png"
traitFile = args[2] # traitFile="data/traits/gap_gm.txt"
traitName = args[3] # traitName="gap_gm"
traitLabel = args[4] # traitLabel="Grey matter BAG"
covsFile = args[5] # covsFile="data/traits/covs.txt"
covs = args[6] # covs="sex,age,age2,ac1,ac2,array,TIV"
PDadjust = as.numeric(args[7]) # PDadjust=0.375
PDsize = as.numeric(args[8]) # PDsize=0.5
ymin = as.numeric(args[9]) # ymin = -15
ymax = as.numeric(args[10]) # ymax = 15
ysteps = as.numeric(args[11]) # ysteps = 5
width = as.numeric(args[12]) # width=8.77
height = as.numeric(args[13]) # height=5.96

message('\n--- plot ancestry components vs. trait of interest | settings ---',
        '\noutFile: ', outFile,
        '\ntraitFile: ', traitFile,
        '\ntraitName: ', traitName,
        '\ntraitLabel: ', traitLabel,
        '\ncovsFile: ', covsFile,
        '\ncovs: ', covs,
        '\nPDadjust: ', PDadjust,
        '\nPDsize: ', PDsize,
        '\nymin: ', ymin,
        '\nymax: ', ymax,
        '\nysteps: ', ysteps,
        '\nwidth: ', width,
        '\nheight: ', height,'\n')

# attach required packages and create target directory
for (pkg in c('dplyr','ggplot2','ggpointdensity','hexbin','lsr','patchwork','stringr','viridis')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# transform variables
covs = str_split(covs, ',')[[1]]

# load data
message(' - loading data.')
input1 = read.delim(traitFile, header = T)
input2 = read.delim(covsFile, header = T)
df = inner_join(input1, input2, by = c("FID", "IID"))

# residualize variables
df[['trait_res']] = lm(df[,c(traitName,covs)]) %>% residuals()
df[['PC1_res']] = lm(df[,c('PC1',covs)]) %>% residuals()
df[['PC2_res']] = lm(df[,c('PC2',covs)]) %>% residuals()
df[['PC3_res']] = lm(df[,c('PC3',covs)]) %>% residuals()
df[['PC4_res']] = lm(df[,c('PC4',covs)]) %>% residuals()

# Create plots with regression stats
message(' - generating PC1-PC4 plots.')
pcs = paste0("PC", 1:4)
plot_list = lapply(pcs, function(pc) {
  
  # get statistics
  mdl = lm(df[,c(traitName,pc,covs)])
  mdl_summary = summary(mdl)
  beta = coef(mdl_summary)[2, "Estimate"]
  se = coef(mdl_summary)[2, "Std. Error"]
  p_value = coef(mdl_summary)[2, "Pr(>|t|)"]
  p_value_formatted <- ifelse(p_value < 0.001, format(p_value, scientific = TRUE, digits = 1), sprintf("%.3f", p_value))
  eta_sq = etaSquared(mdl)[1,"eta.sq.part"]
  eta_sq_formatted <- ifelse(eta_sq < 0.001, format(eta_sq, scientific = TRUE, digits = 1), sprintf("%.3f", eta_sq))
  stat_text = sprintf("BETA = %.2f ± %.2f\np = %s | ηp² = %s", beta, se, p_value_formatted, eta_sq_formatted)

  # plot relationship of residualized variables
  trait_res = lm(df[,c(traitName,covs)]) %>% residuals()
  pc_res = lm(df[,c(pc,covs)]) %>% residuals()
  plot = ggplot(df, aes(x = !!pc_res, y = !!trait_res)) +  
    geom_pointdensity(adjust = PDadjust, size = PDsize) +
    ggtitle(stat_text) + 
    xlab(pc) +
    geom_smooth(method = "lm", formula = y ~ x) +
    scale_color_viridis() +
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
    scale_y_continuous(breaks = seq(ymin,ymax,ysteps)) +
    theme_bw(base_size = 10) +
    theme(plot.margin = unit(c(5, 5, 5, 0), "mm"),
          panel.border = element_rect(linewidth = 0.25),
          legend.position = 'none',
          line = element_line(linewidth = 0.25),
          plot.title = element_text(size = 7, hjust = 0.5, vjust = 0, face = "plain", lineheight = 1.2),
          axis.title.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = -0.5),
          axis.title.y = element_text(size = 12, angle = 90, hjust = 0.5),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.ticks = element_line(linewidth = 0.25))
  
  # Add y-axis label only for the first plot (PC1)
  if (pc == "PC1") {
    plot <- plot + ylab(traitLabel)
  } else {
    plot <- plot + theme(axis.title.y = element_blank())  # Remove y-axis label
  }
  
})
names(plot_list) <- pcs

# plot pc1 x pc2
message(' - generating PC1xPC2 plots.')
plot_list[["PC1xPC2"]] = ggplot(df, aes(x = PC1_res, y = PC2_res, z = trait_res)) +
  stat_summary_hex(fun = mean, bins = 50) +  # Computes mean in each hexbin
  scale_fill_continuous(name = "BAG\n(years)", na.value = "white", limits = c(-13,13), breaks = seq(-10,10,5)) +
  labs(x = "PC1", y = "PC2") +
  theme_bw(base_size = 10) +
  theme(
    plot.margin = unit(c(0, 0, 5, 0), "mm"),
    panel.border = element_rect(linewidth = 0.25),
    legend.position = 'right',
    legend.margin=margin(t = 50, unit='mm'),
    legend.title = element_text(hjust = 0.5),
    line = element_line(linewidth = 0.25),
    plot.title = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 1, face = 'bold'),
    axis.title.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = -0.5),
    axis.title.y = element_text(size = 12, angle = 90, hjust = 0.5),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.ticks = element_line(linewidth = 0.25)
  )

# plot pc3 x pc4
message(' - generating PC3xPC4 plots.')
plot_list[["PC3xPC4"]] = ggplot(df, aes(x = PC3_res, y = PC4_res, z = trait_res)) +
  stat_summary_hex(fun = mean, bins = 50) +  # Computes mean in each hexbin
  scale_fill_continuous(name = "BAG\n(years)", na.value = "white", limits = c(-13,13), breaks = seq(-10,10,5)) +
  labs(x = "PC3", y = "PC4") +
  theme_bw(base_size = 10) +
  theme(
    plot.margin = unit(c(0, 0, 5, 0), "mm"),
    panel.border = element_rect(linewidth = 0.25),
    legend.position = 'right',
    legend.margin=margin(t = 50, unit='mm'),
    legend.title = element_text(hjust = 0.5),
    line = element_line(linewidth = 0.25),
    plot.title = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 1, face = 'bold'),
    axis.title.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = -0.5),
    axis.title.y = element_text(size = 12, angle = 90, hjust = 0.5),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.ticks = element_line(linewidth = 0.25)
  )

# arrange plots
message(' - arranging plots.')
layout = "
ABCD
ABCD
EEFF
EEFF
EEFF
"
pl = wrap_plots(plot_list, design = layout, guides = 'collect') +
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )

# save file
message(sprintf(' - saving %s',outFile))
png(width = width, height = height, units = "in", res = 300, filename = outFile)
pl
invisible(dev.off())
message('--- Completed: plot ancestry components vs. trait of interest ---')
