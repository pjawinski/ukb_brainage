#!/usr/bin/env Rscript

# ===================================================================
# === draw qq plots of genetic correlations (Neale et al. traits) ===
# ===================================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=10) {
  stop(paste0('expected 10 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
inputFile = args[1] # inputFile="results/combined/rgNeale.txt"
outFile = args[2] # outFile="results/combined/rgNeale.qq.png"
rgCol = args[3] # rgCol="gap_gm_rg"
pCol = args[4] # pCol="gap_gm_p"
ylim = as.numeric(args[5]) # upper y axis limit; ylim = 9
ysteps = as.numeric(args[6]) # y axis breaks; ysteps = 3
xlim = as.numeric(args[7]) # upper y axis limit; xlim = 3
xsteps = as.numeric(args[8]) # y axis breaks; xsteps = 1
width = as.numeric(args[9]) # plot width; width = 4.79
height = as.numeric(args[10]) # plot height, height = 4.57

message(paste0('\n--- Create qq plots of genetic correlations (Neale et al. traits) ---',
               '\ninputFile: ', inputFile,
               '\noutFile: ', outFile,
               '\nrgCol: ', rgCol,
               '\npCol: ', pCol,
               '\nylim: ', ylim,
               '\nysteps ', ysteps,
               '\nxlim: ', xlim,
               '\nxsteps ', xsteps,
               '\nwidth: ', width,
               '\nheight: ', height,'\n'))

# attach required packages
for (pkg in c('dplyr','stringr','ggplot2','patchwork')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# load data
df = read.delim(inputFile, header = T, sep = '\t')

# make plot
df$rg = df[,rgCol]
df$p_log10 = -log10(df[,pCol])

# calculate FDR < 0.05
df$fdr = p.adjust(df[,pCol], method = 'BH')
df$fdrLogical = df$fdr < 0.05

# order of custom categories
df$custom_category = factor(df$custom_category,
                            levels = c('Biological samples','Diet by 24-hour recall','Maternity and sex-specific factors','Medical history and conditions','Medications','Physical measures','Mental health','Cognitive function','Family history and early life factors','Lifestyle and environment','Sociodemographics','Work environment'),
                            labels = c('Biological samples', 'Diet by 24-hour recall', 'Maternity and sex-specific factors','Medical history and conditions','Medications','Physical measures','Mental health','Cognitive function','Family history and early life factors','Lifestyle and environment','Sociodemographics','Work environment'))

# -----------------------
# --- draw rg qq plot ---
# -----------------------

# create new data frame for plot
dfplot = df
dfplot = dfplot[order(dfplot$custom_category, dfplot[,'trait_description']),]
ncategories = length(unique(dfplot$custom_category))

# define qq plot function
qq.plot = function(pvalues, qq_title, ylim, type) {
  #error checking
  if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
  if(!(class(pvalues)=="numeric" || 
       (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
    stop("pvalue vector is not numeric, can't draw plot")
  if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
  if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
  if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
  
  # calculated expected values, transform to -log10 scale
  n = length(pvalues)+1
  exp.x = -log10((rank(pvalues, ties.method="first")-.5)/n) # exp.x2 = -log10((rank(pvalues2, ties.method="first")-.5)/n)
  pvalues = -log10(pvalues) # pvalues2 = -log10(pvalues2)
  
  #this is a helper function to draw the confidence interval
  conf.points=1000
  conf.alpha = 0.05
  conf.points = min(conf.points, n-1);
  mpts<-matrix(nrow=conf.points*2, ncol=2)
  for(i in seq(from=1, to=conf.points)) {
    mpts[i,1]<- -log10((i-.5)/n)
    mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
    mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
    mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
  }
  
  # define dataframes for ggplot
  qqdata = data.frame(pvalues, exp.x) # qqdata = data.frame(pvalues, pvalues2, exp.x, exp.x2)
  confdata = data.frame(x = mpts[,1], y = mpts[,2])
  
  # adjust names
  if (qq_title == "Personality and psychosocial factors") {
    qq_title = "Personality and\npsychosocial factors"
  } else if
  (qq_title == "Sexual and sex-specific factors") {
    qq_title = "Sexual and\nsex-specific factors"
  } else
    qq_title = paste0("\n",qq_title)
  
  # make ggplot
  if (type == "xaxis") {
    qqtheme = theme(plot.margin=margin(0,0,0,0),
                    plot.title = element_text(hjust = 0.5, vjust = 0, size = 8, margin = margin(t = 0, r = 0, b = 5, l = 0)),
                    axis.title = element_text(size = 8),
                    axis.title.x = element_text(size=8, margin = margin(t = 5, r = 0, b = 0, l = 0)),
                    axis.title.y = element_blank(), # axis.title.y = element_text(size=8,margin = margin(t = 0, r = 1, b = 0, l = 0)),
                    axis.text.x = element_text(size=8,angle = 0, hjust=0.5, vjust=1),
                    axis.text.y = element_text(size=8),
                    axis.ticks.length=unit(.1, "cm"),
                    axis.ticks = element_line(size = 0.25),
                    legend.title = element_blank(),
                    legend.key.size = unit(0.5, "cm"),
                    panel.grid.major = element_blank(), #colour = "grey92"
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank()) #colour = "grey92"
  } else if
  (type == "yaxis") {
    qqtheme = theme(plot.margin=margin(0,0,0,0),
                    plot.title = element_text(hjust = 0.5, vjust = 0, size = 8, margin = margin(t = 0, r = 0, b = 5, l = 0)),
                    axis.title = element_text(size = 8),
                    axis.title.x = element_blank(), # axis.title.x = element_text(size=8, margin = margin(t = 5, r = 0, b = 0, l = 0)),
                    axis.title.y = element_text(size=8,margin = margin(t = 0, r = 5, b = 0, l = 0)),
                    axis.text.x = element_text(size=8,angle = 0, hjust=0.5, vjust=1),
                    axis.text.y = element_text(size=8),
                    axis.ticks.length=unit(.1, "cm"),
                    axis.ticks = element_line(size = 0.25),
                    legend.title = element_blank(),
                    legend.key.size = unit(0.5, "cm"),
                    panel.grid.major = element_blank(), #colour = "grey92"
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank()) #colour = "grey92"
  } else if
  (type == "both") {
    qqtheme = theme(plot.margin=margin(0,0,0,0),
                    plot.title = element_text(hjust = 0.5, vjust = 0, size = 8, margin = margin(t = 0, r = 0, b = 5, l = 0)),
                    axis.title = element_text(size = 8),
                    axis.title.x = element_text(size=8, margin = margin(t = 5, r = 0, b = 0, l = 0)),
                    axis.title.y = element_text(size=8,margin = margin(t = 0, r = 5, b = 0, l = 0)),
                    axis.text.x = element_text(size=8,angle = 0, hjust=0.5, vjust=1),
                    axis.text.y = element_text(size=8),
                    axis.ticks.length=unit(.1, "cm"),
                    axis.ticks = element_line(size = 0.25),
                    legend.title = element_blank(),
                    legend.key.size = unit(0.5, "cm"),
                    panel.grid.major = element_blank(), #colour = "grey92"
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank()) #colour = "grey92"
  } else if
  (type == "none") {
    qqtheme = theme(plot.margin=margin(0,0,0,0),
                    plot.title = element_text(hjust = 0.5, vjust = 0, size = 8, margin = margin(t = 0, r = 0, b = 5, l = 0)),
                    axis.title = element_text(size = 8),
                    axis.title.x = element_blank(), # axis.title.x = element_text(size=8, margin = margin(t = 5, r = 0, b = 0, l = 0)),
                    axis.title.y = element_blank(), # axis.title.y = element_text(size=8,margin = margin(t = 0, r = 1, b = 0, l = 0)),
                    axis.text.x = element_text(size=8,angle = 0, hjust=0.5, vjust=1),
                    axis.text.y = element_text(size=8),
                    axis.ticks.length=unit(.1, "cm"),
                    axis.ticks = element_line(size = 0.25),
                    legend.title = element_blank(),
                    legend.key.size = unit(0.5, "cm"),
                    panel.grid.major = element_blank(), #colour = "grey92"
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank()) #colour = "grey92"
  }
  
  ggplot(data=qqdata, aes(x=exp.x, y=pvalues)) +
    ggtitle(qq_title) +
    geom_polygon(data = confdata, aes(x = x, y = y), cex = 1, fill = "grey90", show.legend = FALSE, alpha = 1) +
    geom_line(aes(x = exp.x, y = exp.x)) + 
    geom_point(cex = 0.5, col=c("royalblue3"), show.legend = FALSE, alpha = 1) +
    #geom_text(data = data.frame("NA"), label = "Mental health", x=1.5, y=11, hjust = 0.5, size = 4, fontface='bold') +
    #geom_point(aes(x=exp.x2, y=pvalues2), cex = 0.5, col=c("red"), show.legend = FALSE, alpha = 1) +
    #stat_smooth(aes(x = exp.x, y = exp.x), method="lm",fullrange=TRUE) +
    scale_x_continuous(limits = c(0, 2.8), breaks = seq(0,3,1), name = expression("expected -log"[10]*"("*italic(p)*")")) + # gsub("\\.", " ",
    scale_y_continuous(limits = c(0, ylim), breaks = seq(0,ylim,3), name = expression("observed -log"[10]*"("*italic(p)*")")) + # gsub("\\.", " ",
    geom_segment(aes(x=0,xend=2.8,y=-Inf,yend=-Inf), colour = "black", size = 0.25) +
    geom_segment(aes(y=0,yend=ylim,x=-Inf,xend=-Inf), colour = "black", size = 0.25) +
    theme_bw() +
    qqtheme
}
  

# create qq plots 
qq_biosamples = qq.plot(pvalues = dfplot[dfplot$custom_category=='Biological samples', pCol], qq_title = 'Biological', ylim = 9, type = 'yaxis')
qq_medications = qq.plot(pvalues = dfplot[dfplot$custom_category=='Medications', pCol], qq_title = 'Medications', ylim = 9, type = 'yaxis')
qq_diet = qq.plot(pvalues = dfplot[dfplot$custom_category=='Diet by 24-hour recall', pCol], qq_title = 'Diet by 24-hour', ylim = 9, type = 'none')
qq_maternity = qq.plot(pvalues = dfplot[dfplot$custom_category=='Maternity and sex-specific factors', pCol], qq_title = 'Maternity', ylim = 9, type = 'none')
qq_medical = qq.plot(pvalues = dfplot[dfplot$custom_category=='Medical history and conditions', pCol], qq_title = 'Medical history', ylim = 9, type = 'none')
qq_physical = qq.plot(pvalues = dfplot[dfplot$custom_category=='Physical measures', pCol], qq_title = 'Physical', ylim = 9, type = 'none')
qq_mental = qq.plot(pvalues = dfplot[dfplot$custom_category=='Mental health', pCol], qq_title = 'Mental health', ylim = 9, type = 'none')
qq_cognitive = qq.plot(pvalues = dfplot[dfplot$custom_category=='Cognitive function', pCol], qq_title = 'Cognitive', ylim = 9, type = 'none')
qq_familiy = qq.plot(pvalues = dfplot[dfplot$custom_category=='Family history and early life factors', pCol], qq_title = 'Family history', ylim = 9, type = 'both')
qq_lifestyle = qq.plot(pvalues = dfplot[dfplot$custom_category=='Lifestyle and environment', pCol], qq_title = 'Lifestyle', ylim = 9, type = 'xaxis')
qq_socio = qq.plot(pvalues = dfplot[dfplot$custom_category=='Sociodemographics', pCol], qq_title = 'Sociodemographics', ylim = 9, type = 'xaxis')
qq_work = qq.plot(pvalues = dfplot[dfplot$custom_category=='Work environment', pCol], qq_title = 'Work environment', ylim = 9, type = 'xaxis')

# merge plots
qqplots = qq_biosamples + qq_diet + qq_maternity + qq_medical + 
qq_medications + qq_physical + qq_mental + qq_cognitive +
qq_familiy + qq_lifestyle + qq_socio + qq_work +
   plot_layout(ncol = 4)

# save plot  
message(sprintf('Saving %s', outFile))
ggsave(outFile, qqplots, width = width, height = height, units = "in", dpi = 300)
system(sprintf('chmod 770 %s', outFile))
message('-- Completed: Create qq plots of genetic correlations (Neale et al. traits) ---')


