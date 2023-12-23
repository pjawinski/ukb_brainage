#!/usr/bin/env Rscript

# =======================================================================
# === draw volcano plot of genetic correlations (Neale et al. traits) ===
# =======================================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=12) {
  stop(paste0('expected 12 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
inputFile = args[1] # inputFile="results/combined/rgNeale.txt"
outFile = args[2] # outFile="results/combined/rgNeale.volcano.forest.png"
rgCol = args[3] # rgCol="gap_gm_rg"
seCol = args[4] # seCol="gap_gm_se"
pCol = args[5] # pCol="gap_gm_p"
multipleTesting = args[6] # 'fdr' or 'bonferroni' or 'both'; multipleTesting = 'fdr'
ylim = as.numeric(args[7]) # upper y axis limit; ylim = 10
ysteps = as.numeric(args[8]) # y axis breaks; ysteps = 2
xlim = as.numeric(args[9]) # upper y axis limit; xlim = 0.5
xsteps = as.numeric(args[10]) # y axis breaks; xsteps = 0.25
width = as.numeric(args[11]) # plot width; width = 7.56
height = as.numeric(args[12]) # plot height, height = 8.58

message(paste0('\n--- Create volcano and forest plot of genetic correlations (Neale et al. traits) ---',
               '\ninputFile: ', inputFile,
               '\noutFile: ', outFile,
               '\nrgCol: ', rgCol,
               '\nseCol: ', seCol,
               '\npCol: ', pCol,
               '\nmultipleTesting: ', multipleTesting,               
               '\nylim: ', ylim,
               '\nysteps ', ysteps,
               '\nxlim: ', xlim,
               '\nxsteps ', xsteps,
               '\nwidth: ', width,
               '\nheight: ', height,'\n'))

# attach required packages
for (pkg in c('dplyr','stringr','ggplot2','ggrepel','patchwork','RColorBrewer')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# load data
df = read.delim(inputFile, header = T, sep = '\t')

# make plot
df$rg = df[,rgCol]
df$se = df[,seCol]
df$p_log10 = -log10(df[,pCol])

# calculate FDR < 0.05
df$fdr = p.adjust(df[,pCol], method = 'BH')
df$fdrLogical = df$fdr < 0.05

# order of custom categories
df$custom_category = factor(df$custom_category,
                            levels = c('Biological samples','Diet by 24-hour recall','Maternity and sex-specific factors','Medical history and conditions','Medications','Physical measures','Mental health','Cognitive function','Family history and early life factors','Lifestyle and environment','Sociodemographics','Work environment'),
                            labels = c('Biological samples', 'Diet by 24-hour recall', 'Maternity and sex-specific factors','Medical history and conditions','Medications','Physical measures','Mental health','Cognitive function','Family history and early life factors','Lifestyle and environment','Sociodemographics','Work environment'))

# ----------------------------
# --- draw rg volcano plot ---
# ----------------------------

# create new data frame for plot
dfvolc = df
dfvolc = dfvolc[order(dfvolc$custom_category, dfvolc[,'trait_description']),]
ncategories = length(unique(dfvolc$custom_category))

# set ylim significance threshold
if (multipleTesting == 'fdr') {
  sig = min(dfvolc$p_log10[dfvolc$fdr < 0.05])
  } else if (multipleTesting == 'bonferroni' | multipleTesting == 'both') {
  sig = -log10(0.05/nrow(dfvolc))
}

# set pvalue of < ylim to ylim
dfvolc$p_log10[dfvolc$p_log10 > ylim] = ylim

# create color palette and shapes
mypalette = c(rbind(brewer.pal(6,"BrBG"),brewer.pal(6,"BrBG"))) # mypalette = c(rep('red',6),rep('#4682B4',6))
myshapes = rep(c(3,4,0,6,19,15),2) # myshapes = c(18,19,3,15,16,17,15,16,17,18,19,3)

# create plot
vp = ggplot(data=dfvolc, aes(x=rg, y=p_log10)) +
  geom_point(aes(col=custom_category, shape=custom_category), cex = 1, show.legend = TRUE, alpha = 1) +
  scale_color_manual(values = mypalette) +
  scale_shape_manual(values = myshapes ) + 
  scale_x_continuous(limits = c(-xlim, xlim), name = expression("Genetic Correlation ("*r[G]* ")")) +
  scale_y_continuous(limits = c(0, ylim), breaks = seq(0,ylim,ysteps), name = expression("-log"[10]*"("*italic(p)*")")) +
  geom_segment(aes(x=-xlim,xend=xlim,y=-Inf,yend=-Inf), colour = "black", size = 0.25) +
  geom_segment(aes(y=0,yend=ylim,x=-Inf,xend=-Inf), colour = "black", size = 0.25) +
  geom_hline(yintercept=sig, lty=3, size = 0.25) +
  theme_bw() +
  theme(plot.margin=margin(5,5,5,5),
        plot.title = element_text(hjust = 1, size = 8),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(size = 12, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 12,margin = margin(t = 0, r = 1, b = 0, l = 0)),
        axis.text.x = element_text(size = 10,angle = 0, hjust=0.5, vjust=1),
        axis.text.y = element_text(size = 10),
        axis.ticks = element_line(size = 0.25),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())

# ---------------------------
# --- draw rg forest plot ---
# ---------------------------

# define example phenotypes
traits = data.frame(trait_description = c('Average total household income before tax',
                                          'Maternal smoking around birth',
                                          'Number of treatments/medications taken',
                                          'Mother\'s age at death',
                                          'Father\'s age at death',
                                          'Frequency of tiredness / lethargy in last 2 weeks',
                                          'Number of self-reported non-cancer illnesses',
                                          'Hearing difficulty/problems with background noise',
                                          'Drive faster than motorway speed limit',
                                          'Medication for cholesterol, blood pressure or diabetes: None of the above',
                                          'Townsend deprivation index at recruitment',
                                          'Neutrophill count',
                                          'Heel bone mineral density (BMD)',
                                          'Vascular/heart problems diagnosed by doctor: None of the above',
                                          'Long-standing illness, disability or infirmity',
                                          'Qualifications: None of the above',
                                          'Age completed full time education',
                                          'Types of transport used (excluding work): Public transport',
                                          'Particulate matter air pollution (pm2.5); 2010',
                                          'White blood cell (leukocyte) count',
                                          'Frequency of unenthusiasm / disinterest in last 2 weeks',
                                          'Frequency of tenseness / restlessness in last 2 weeks',
                                          'Recent restlessness'))

# create data.frame with selected traits
dfforest = left_join(traits, df, by = "trait_description")
dfforest$lower = dfforest$rg + dfforest$se 
dfforest$upper = dfforest$rg - dfforest$se 

# sort by rG
dfforest = dfforest[order(dfforest$rg),]

# get sign for colors
dfforest$trait_description <- factor(dfforest$trait_description, levels=rev(dfforest$trait_description))
dfforest$sign = "sign"
dfforest$sign[dfforest$rg < 0] = "negative"
dfforest$sign[dfforest$rg > 0] = "positive"

# make forrest plot
fp = ggplot(data=dfforest, aes(x=trait_description, y=rg, ymin=lower, ymax=upper, label = trait_description)) +
  geom_errorbar(aes(ymin=lower, ymax=upper),width=0.25,cex=0.25,show.legend = FALSE) +
  geom_point(aes(col=sign), size = 1.8, show.legend = FALSE) +
  scale_color_manual(values=c('steelblue','red')) +
  geom_text(size = 2.5, hjust = 0.5+sign(dfforest$rg)/2, nudge_y =(-sign(dfforest$rg)*0.05)-dfforest$rg) +
  geom_hline(yintercept=0, lty=1, size = 0.25, color = "black") +
  coord_flip(clip = 'off') +
  scale_x_discrete(name = "", position = 'bottom') +
  scale_y_continuous(limits = c(-xlim, xlim), breaks = seq(-xlim,xlim,xsteps), name = expression("Genetic Correlation ("*r[G]* ")")) +
  geom_segment(aes(y=-0.5,yend=0.5,x=-Inf,xend=-Inf), colour = "black", size = 0.25) +
  theme_bw() +
  theme(plot.margin=margin(30,5,5,5),
        plot.title = element_text(hjust = 1, size = 8),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=12,margin = margin(t = 0, r = 1, b = 0, l = 0)),
        axis.text.x = element_text(size=10,angle = 0, hjust=0.5, vjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())

# ---------------------
# --- combine plots ---
# ---------------------

# combine plots
merged = vp/fp + plot_layout(heights = c(2, 3))

# save plot
message(sprintf('Saving %s', outFile))
ggsave(outFile, merged, width = width, height = height, units = "in", dpi = 300)
system(sprintf('chmod 770 %s', outFile))
message('-- Completed: Create volcano and forest plot of genetic correlations (Neale et al. traits) ---')
