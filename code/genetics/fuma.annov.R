#!/usr/bin/env Rscript

# =================================================
# === Create figure of FUMA annov.stats results ===
# =================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop(paste0('expected 4 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
annovarStats = args[1] # annovarStats="results/gap_gm/fuma/FUMA_job174874.zip,results/gap_wm/fuma/FUMA_job174875.zip,results/gap_gwm/fuma/FUMA_job174876.zip"
outputFile = args[2] # outputFile="results/combined/fuma.annov.stats.png"
width = as.numeric(args[3]) # width = 8.04
height = as.numeric(args[4]) # height = 7.46

message(paste0('\n--- Settings for FUMA annov.stats Figure ---',
               '\nannovarStats: ', annovarStats,
               '\noutputFile: ', outputFile,
               '\nwidth: ', width,
               '\nheight: ', height,'\n'))

# attach packages to current R session
for (pkg in c('dplyr','ggplot2','ggrepel','patchwork','stringr')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# transform variables
annovarStats = str_split(annovarStats, ',')[[1]]
outputFile = str_split(outputFile, ',')[[1]]

# load data and determine maximum enrichment
for (i in 1:length(annovarStats)) {
  message(paste0('[', i, '/', length(annovarStats), '] ', annovarStats[i], '...'))
  tmp = read.delim(unz(annovarStats[i], 'annov.stats.txt'), sep = '\t', head = TRUE)
  assign(paste0('annov.stats',i),tmp)
  if (i == 1) { maxEnrich = max(tmp$enrichment) } else { maxEnrich = max(c(maxEnrich,tmp$enrichment)) }
}

# define colors for pie chart
cols = c('#68A042','#318466','#2D496B','#4678A9','#4D97B7','#2874A6','#59485D','#81374B','#A61E28','#CB3228','#E95361','#D8611E','#FF64AB')

# make plots and merge them
for (i in 1:length(annovarStats)) {

  # load data
  annov.stats = get(paste0('annov.stats',i))

  # order categories and set as factor
  df = data.frame(annot = c('UTR3', 'UTR5','downstream','exonic','intergenic','intronic','ncRNA_exonic','ncRNA_intronic','ncRNA_splicing','splicing','upstream'))
  df = left_join(df, annov.stats, 'annot')
  df$annot = factor(df$annot, levels = df$annot)

  # create labels
  df$pct = paste0(sprintf("%.1f", df$prop*100),"%")
  df$pct[df$count == 0] = "NA"
  #df$label = df$annot %>% paste("(") %>% paste0(df$pct) %>% paste0(")") %>% str_replace(pattern = "_", replacement = " ")
  df$label = df$annot %>% str_replace(pattern = "_", replacement = " ")
  #df$pielabel = df$annot %>% paste("\n(") %>% paste0(df$pct) %>% paste0(")") %>% str_replace(pattern = "_", replacement = " ")
  #df$pielabel = df$annot %>% str_replace(pattern = "_", replacement = " ")
  df$pielabel = df$pct
  df$enrLabel = "non-significant"
  df$enrLabel[df$enrichment > 1 & df$fisher.P < 0.05] = "enriched"
  df$enrLabel[df$enrichment < 1 & df$fisher.P < 0.05] = "depleted"

  # create description
  df$p = sprintf("%.3f", df$fisher.P)
  df$p[df$fisher.P < 0.001] = sprintf("%.0g", df$fisher.P[df$fisher.P < 0.001])
  df$p = paste("p =", df$p)
  df$p[df$p == "p = 0"] = "p < 5e-324"
  df$description = "FE =" %>% paste(sprintf("%.2f",df$enrichment)) %>% paste0("; ") %>% paste(df$p)
  df$description = df$p
  df$description[df$count == 0] = ""

# plot enrichment
  enrich = ggplot(df, aes(y = annot, x = enrichment, fill = enrLabel)) +
    geom_bar(stat="identity", width = 0.8) +
    scale_fill_manual("enrLabel", values = c(enriched = "#2874A6", depleted = "#81374B", `non-significant` = "#A6A6A6")) +
    geom_vline(xintercept=1, linetype="dashed", color = "black", size = 0.25) +
    coord_cartesian(xlim = c(0, max(maxEnrich*1.03)), expand = FALSE) +
    xlab('Fold Enrichment') +
    scale_y_discrete(limits = rev(levels(df$annot)), labels = rev(df$label)) +
    scale_x_continuous(breaks = seq(0,floor(maxEnrich),1)) +
    #geom_segment(aes(y=0.6,yend=Inf,x=-Inf,xend=-Inf), colour = "black", size = 0.25)+
    #geom_segment(aes(x=0,xend=floor(maxEnrich),y=-Inf,yend=-Inf), colour = "black", size = 0.25) +
    theme_bw(base_family = "sans") +
    theme(plot.margin=margin(15,5,0,0),
          panel.border = element_rect(color = "black"),
          panel.grid.major = element_line(size = 0.25),
          panel.grid.minor = element_line(size = 0.25),
          axis.line = element_blank(), # axis.line = element_line(colour = "black", size = 0.25),
          plot.title = element_text(hjust = 1, size = 8),
          legend.position = 'bottom',
          legend.text = element_text(size=8),
          legend.key.size = unit(1,"line"),
          legend.margin = margin(0,20,0,0),
          axis.title = element_text(size = 10),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size=9,margin = margin(t = 3, r = 0, b = 0, l = 0)),
          axis.text.y = element_text(size=8,angle = 0, hjust=1, vjust=0.5),
          axis.text.x = element_text(size=8),
          axis.ticks.y = element_blank(),
          axis.ticks = element_line(colour = "black", size = 0.25),
          legend.title = element_blank())
  
  # pvals right to enrich plot
  description = ggplot(df, aes(y = annot, x = enrichment)) +
    geom_bar(stat="identity", width = 0.8, fill = c("#FFFFFF")) +
    scale_fill_manual(values = c("FFFFFF")) + 
    coord_cartesian(xlim = c(0, max(df$enrichment*1.03)), expand = FALSE) +
    xlab('Fold Enrichment') +
    scale_y_discrete(limits = rev(levels(df$annot)), labels = rev(df$description)) +
    theme_bw(base_family = "sans") +
    theme(plot.margin=margin(15,0,0,0),
          panel.grid = element_blank(),
          plot.title = element_blank(),
          legend.position = 'none',
          legend.text = element_text(size=8),
          legend.key.size = unit(1,"line"),
          legend.margin = margin(0,20,0,0),
          panel.border = element_blank(),
          axis.title = element_blank(),
          axis.title.x = element_text(size=9,margin = margin(t = 3, r = 0, b = 0, l = 0),color = c("#FFFFFF")),
          axis.text.y = element_text(size=8,angle = 0, hjust=0, vjust=0.5),
          axis.text.x = element_text(size=8,color = c("#FFFFFF")),
          axis.ticks.y =  element_blank(),
          axis.ticks.x =  element_blank(),
          axis.ticks = element_line(size = 0.25),
          legend.title = element_blank())

  # pie chart
  df = df %>% mutate(lab.ypos = 1-(cumsum(prop) - 0.5*prop))
  pie = ggplot(df, aes(x = "", y = prop, fill = annot)) +
    geom_bar(stat = "identity", color = "white") +
    coord_polar(theta = "y", start = 0, direction = -1, clip = 'off') +
    scale_fill_manual(values = cols[1:11], labels=df$label) + 
    geom_text(data = df[df$prop > 0.05,], aes(x = 1.15, y = lab.ypos, label = pielabel), color = "white", size = 2.8) +
    geom_text_repel(
      data          = df[df$count!=0 & df$prop <= 0.05,],
      aes(x=1.44, y = lab.ypos, label = pielabel),
      size          = 2.5,# 1.5
      point.padding = 0,
      segment.size  = 0.2 , # 0.2
      segment.color = "grey30",
      segment.curvature = 0,
      direction = "both", # direction = "y"
      nudge_x = 0.5,
      min.segment.length = 0,
      xlim = c(1.4, NA)
    ) +
    theme_void(base_family = "sans") +
    theme(plot.margin=margin(0,0,0,0),
          axis.title = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size=8),
          legend.key.size = unit(1,"line"),
          legend.key.height = unit(0.75,"line"))
  
  # collect plots
  if (i < length(annovarStats)) { 
    enrich = enrich + guides(fill='none') 
    pie = pie + guides(fill='none')
  }
  tmp = enrich + description + pie + plot_layout(widths = c(13,1,8))
  if (i == 1) { pl = tmp } else { pl = pl / tmp }
  
}

# set layout
if (length(annovarStats) > 1) {
  pl = pl + plot_annotation(tag_levels = list(c(rbind(letters,' ', ' ')))) & 
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(size = 14))
}
pl

# plot together
message(paste0('Saving ', outputFile, ' ...'))
png(filename = paste0(outputFile), width = width, height = height, units = "in", res = 600)
pl
invisible(dev.off())
system(paste0('chmod 770 ', outputFile))
message('-- Creating figure of FUMA annov.stats results completed. ---')


