#!/usr/bin/env Rscript

# ======================================
# === draw gene-based manhattan plot ===
# ======================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 8) {
  stop(paste0('expected 8 arguments, but ', length(args), ' arguments provided.'), call.=FALSE)
}

# get arguments from command line
trait = args[1] # trait='gap_gm'
targetDir = args[2] # targetDir='results/gap_gm/fastbat/'
fastbatResults = args[3] # fastbatResults='results/gap_gm/fastbat/fastbat.0kb.gene.fastbat.annovar.clumped'
annotationThresh = args[4] # annotationThresh = 'bonferroni' # 'fdr'
ylim = as.numeric(args[5]) # ylim = 35
ysteps = as.numeric(args[6]) # ysteps = 5
width = as.numeric(args[7]) # width = 11
height = as.numeric(args[8]) # height = 4

message(paste0('\n--- Manhatten plot settings ---',
               '\ntrait: ', trait,
               '\ntargetDir: ', targetDir,
               '\nfastbatResults: ', fastbatResults,
               '\nannotationThresh: ', annotationThresh,
               '\nylim: ', ylim,
               '\nysteps: ', ysteps,
               '\nwidth: ', width,
               '\nheight: ', height,'\n'))

# load packages
for (pkg in c('data.table', 'dplyr', 'ggplot2', 'ggrepel')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# start loop
system(paste0('mkdir -p ', targetDir))

# import data
message(paste0('Importing data.'))
GWAS = read.delim(fastbatResults, sep='\t', header=T, stringsAsFactors=FALSE)

# prepare for manhattan plot
names(GWAS)[names(GWAS)=='Chr'] = 'CHR'
names(GWAS)[names(GWAS)=='Start'] = 'BP'
names(GWAS)[names(GWAS)=='Pvalue'] = 'P'
#GWAS$CHR[GWAS$CHR=='X'] = 23
#GWAS$CHR[GWAS$CHR=='Y'] = 24
GWAS$CHR[GWAS$CHR==25] = 23 # GWAS$CHR[GWAS$CHR=='XY'] = 23
GWAS$CHR[GWAS$CHR==26] = 25 # GWAS$CHR[GWAS$CHR=='MT'] = 25
GWAS$CHR = as.numeric(GWAS$CHR)

# sort by chromosome
idx = sort(GWAS$CHR, index.return = T)$ix
GWAS = GWAS[idx,]

# get cumulative base pair positions and center positions
# credits to Daniël Roelfs (http://www.danielroelfs.com/coding/manhattan_plots/)
message(paste0('Geting cumulative base pair positions.'))
nCHR = length(unique(GWAS$CHR))
GWAS$BPcum = NA
s = 0
nbp = c()
for (i in unique(GWAS$CHR)){
  nbp[i] = max(GWAS[GWAS$CHR == i,]$BP)
  GWAS[GWAS$CHR == i,'BPcum'] = GWAS[GWAS$CHR == i,'BP'] + s
  s = s + nbp[i]
}

axis.set = GWAS %>% 
  group_by(CHR) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

# edit GWAS findings exceeding plot
GWAS$P_exceeding = GWAS$P
GWAS$BPcum_exceeding = GWAS$BPcum
GWAS$Gene_exceeding = GWAS$Gene
if (sum(!is.na(GWAS$Gene) & GWAS$P < 10^-ylim) > 0) {
  GWAS[!is.na(GWAS$Gene) & GWAS$P < 10^-ylim,]$P_exceeding = 10^-ylim
  GWAS[!is.na(GWAS$Gene) & GWAS$P < 10^-ylim,]$BPcum_exceeding = GWAS[!is.na(GWAS$Gene) & GWAS$P < 10^-ylim,]$BPcum + 60000000
  GWAS[!is.na(GWAS$Gene) & GWAS$P < 10^-ylim,]$Gene_exceeding = paste0(GWAS[!is.na(GWAS$Gene) & GWAS$P < 10^-ylim,]$Gene, '\n(p = ', sprintf('%0.e',GWAS[!is.na(GWAS$Gene) & GWAS$P < 10^-ylim,]$P), ')')
}
GWAS$P[GWAS$P < 10^-ylim] = 10^-ylim

# set ylim significance threshold
sigFDR = max(GWAS$P[GWAS$FDR < 0.05])
sigBonf = 0.05/nrow(GWAS)

# define Genes that shall be annotated
if (annotationThresh == 'fdr') {
  GWAS$annotation[GWAS$P < sigFDR & GWAS$GENE_COUNT == 1] = 1
} else if (annotationThresh == 'bonferroni') {
  GWAS$annotation[GWAS$P < sigBonf & GWAS$GENE_COUNT == 1] = 1
} else {
  annotationThresh = as.numeric(annotationThresh)
  GWAS$annotation[GWAS$P < annotationThresh & GWAS$GENE_COUNT == 1] = 1
}

# create manhattan plot
message('Creating manhattan plot.')
set.seed(8245)
manhplot = ggplot(data = subset(GWAS, P >= 10^-ylim)) + # 
  geom_point(aes(x=BPcum, y=-log10(P), color=as.factor(CHR)), alpha = 1, size = 0.8, stroke = 0, shape = 16) +
  geom_point(data=subset(GWAS, GENE_COUNT == 1 & P <= sigFDR & P > sigBonf), aes(x=BPcum, y=-log10(P)), color='black', shape=1, size=1.5) +
  geom_point(data=subset(GWAS, GENE_COUNT == 1 & P <= sigBonf & P >= 10^-ylim), aes(x=BPcum, y=-log10(P)), color='black', shape=5, size=3) +
  geom_point(data=subset(GWAS, GENE_COUNT == 1 & P <= 10^-ylim), aes(x=BPcum, y=-log10(P_exceeding)), color='red', shape=5, size=3) +
  scale_color_manual(values = rep(c('#282873', '#6e90ca'), nCHR)) +
  scale_x_continuous(expand = c(0.03,0), label = c(1:18,'', 20, '', 22, 'X', 'Y MT', ''), breaks = axis.set$center) + # label = axis.set$CHR % label = c(1:22, 'X', 'Y MT', '') label = c(1:18,'', 20, '', 22, 'X', 'Y MT', '')
  scale_y_continuous(expand = c(0,0), limits = c(0,ylim), breaks = seq(0,ylim,ysteps)) +
  coord_cartesian(clip = 'off') +
  geom_hline(yintercept = -log10(sigBonf), color = 'black', linetype = 'solid', size = 0.25) +
  geom_hline(yintercept = -log10(sigFDR), color = 'black', linetype = 'dashed', size = 0.25) + 
  labs(x = 'Chromosome', y = expression('-log'[10]*'('*italic(p)*')')) + 
  geom_segment(aes(x=min(axis.set$center),xend=max(axis.set$center),y=-Inf,yend=-Inf), colour = 'black', size = 0.25)+
  geom_segment(aes(y=0,yend=ylim,x=-Inf,xend=-Inf), colour = 'black', size = 0.25) +
  geom_text_repel(
    data          = subset(GWAS, !is.na(annotation) & P >= 10^-ylim),
    aes(x=BPcum, y=-log10(P), label=Gene),
    size          = 3.5,# 1.5
    segment.size  = 0.2 , # 0.2
    segment.color = 'grey30',
    direction = 'both', # direction = 'y'
    nudge_y = 1.5
  ) +
  geom_text_repel(
    data          = subset(GWAS, !is.na(annotation) & P < 10^-ylim),
    aes(x=BPcum_exceeding, y=-log10(P_exceeding), label=Gene_exceeding),
    size          = 3.5,# 1.5
    segment.size  = 0.2 , # 0.2
    segment.color = 'grey30',
    direction = 'y', # direction = 'y'
    nudge_x = 80000000,
    nudge_y = -5,
    arrow = arrow(length = unit(0.02, 'npc'))
  ) +
  theme_bw() +
  theme( 
    legend.position = 'none',
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = 'black', size = 0.25),
    axis.ticks.length=unit(.15, 'cm'),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle = 0, size = 12, vjust = 0.5, margin = margin(t = 3, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(angle = 0, size = 12, vjust = 0.5, margin = margin(t = 0, r = 3, b = 0, l = 0)),
    axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 5, b = 0, l = 5)),
    axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    plot.margin=unit(c(0.25,0.75,0,0),'cm')
  )

# save plot
message('Saving manhattan plot.')
ggsave(paste0(targetDir, '/fastbat.manhattan.png'), manhplot, width = width, height = height, units = 'in', dpi = 300)
system(paste0('chmod -R 770 ', targetDir))
message(paste0('--- Manhatten plot finished ---\n'))
