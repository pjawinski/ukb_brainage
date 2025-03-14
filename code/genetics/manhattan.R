#!/usr/bin/env Rscript

# ===========================
# === draw manhattan plot ===
# ===========================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=12) {
  stop(paste0('expected 12 arguments, but ', length(args), ' arguments provided.'), call.=FALSE)
}

# get arguments from command line
trait = args[1] # trait="gap_gm"
targetDir = args[2] # targetDir="results/gap_gm/replicate/mrmega.all"
sumstats = args[3] # sumstats="results/gap_gm/replicate/mrmega.all/mrmega.weights.gz"
nCriterion = args[4] # nCriterion = FALSE | LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5 
annotationFile = args[5] # annotationFile="results/gap_gm/replicate/mrmega.all/manhattan.annotation.txt"
annotationThresh = as.numeric(args[6]) # annotationThresh=1E-100
sig = as.numeric(args[7]) # sig=5E-8
yend = as.numeric(args[8]) # yend=12
ysteps = as.numeric(args[9]) # ysteps=4
width = as.numeric(args[10]) # width = 10.67
height = as.numeric(args[11]) # height = 4
preview = args[12] # preview = FALSE

logInfo = paste0('\n--- Manhatten plot settings ---',
               '\ntrait: ', trait,
               '\ntargetDir: ', targetDir,
               '\nsumstats: ', sumstats,
               '\nnCriterion: ', nCriterion,
               '\nannotationFile: ', annotationFile,
               '\nannotationThresh: ', annotationThresh,
               '\nsig: ', sig,
               '\nyend: ', yend,
               '\nysteps: ', ysteps,
               '\nwidth: ', width,
               '\nheight: ', height,
               '\npreview: ', preview,'\n')
message(logInfo)

# attach required packages
for (pkg in c('data.table','dplyr','ggplot2','ggrepel')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# start loop
system(paste0('mkdir -p ', targetDir))

# import data
message(paste0('Importing data.'))
GWAS_raw = data.frame(fread(cmd=paste0("gzip -dc ", sumstats), tmpdir = getwd(), header=T, stringsAsFactors=FALSE))
if (annotationFile != "") {
  GWAS_annot = data.frame(fread(paste0(annotationFile), tmpdir = getwd(), sep='\t', header=T, stringsAsFactors=FALSE))
    if ('prioritized' %in% names(GWAS_annot)) {
      GWAS_annot = GWAS_annot[GWAS_annot$TRAIT == trait & GWAS_annot$COJO_pJ < sig,c('ID','prioritized')]
      } else {
      GWAS_annot = GWAS_annot[GWAS_annot$SNP_CNT == 1 & GWAS_annot$P < sig & GWAS_annot$LEAD_SNP_pJ < sig, c("ID", "NEAREST_GENE")]
      }
      names(GWAS_annot) = c('ID','GENE_ANNOT')
}

# exclude snps that do not meet N criterion
if (nCriterion == TRUE) {
  message('Excluding variants that do not meet N criterion.')
  maxN = max(GWAS_raw$N)
  minN = quantile(GWAS_raw$N, probs = 0.9)[[1]] / 1.5
  nExcl = sum(GWAS_raw$N < minN)
  nIncl = sum(GWAS_raw$N >= minN)
  GWAS_raw = GWAS_raw[GWAS_raw$N >= minN,]
  logInfo = paste0(logInfo,
    '\n--- calculations ---',
    '\nmaxN: ', maxN,
    '\nminN: ', minN,
    '\nnExcl: ', nExcl,
    '\nnIncl: ', nIncl)
}

# get gene name of lead snps
if (annotationFile != "") { GWAS = left_join(x = GWAS_raw, y = GWAS_annot, by = "ID") } else { GWAS = data.frame(GWAS_raw, GENE_ANNOT = NA) }

# keep only 1% of variations and lead (for faster plotting when playing around with ggplot parameters)
if (preview == TRUE) {
GWAS = GWAS[rep_len(c(TRUE,rep(FALSE,99)), nrow(GWAS)) | !is.na(GWAS$GENE_ANNOT),]
}

# rename chromosomes
GWAS$CHR[GWAS$CHR=='X'] = 23
GWAS$CHR[GWAS$CHR=='Y'] = 24
GWAS$CHR[GWAS$CHR=='XY'] = 23
GWAS$CHR[GWAS$CHR=='MT'] = 25
GWAS$CHR = as.numeric(GWAS$CHR)

# sort by chromosome
idx = sort(GWAS$CHR, index.return = T)$ix
GWAS = GWAS[idx,]

# get cumulative base pair positions and center positions
# credits to Daniël Roelfs (http://www.danielroelfs.com/coding/manhattan_plots/)
message(paste0('Getting cumulative base pair positions.'))
nCHR = length(unique(GWAS$CHR))
GWAS$BPcum = NA
s = 0
nbp = c()
for (i in unique(GWAS$CHR)){
  nbp[i] = max(GWAS[GWAS$CHR == i,]$BP)
  GWAS[GWAS$CHR == i,"BPcum"] = GWAS[GWAS$CHR == i,"BP"] + s
  s = s + nbp[i]
}

axis.set = GWAS %>% 
  group_by(CHR) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

# truncate GENE_ANNOT at comma character
GWAS$GENE_ANNOT[!is.na(GWAS$GENE_ANNOT)] = sub(",.*", "", GWAS$GENE_ANNOT[!is.na(GWAS$GENE_ANNOT)])

# edit GWAS findings exceeding plot
GWAS$P_exceeding = GWAS$P
GWAS$BPcum_exceeding = GWAS$BPcum
GWAS$GENE_ANNOT_exceeding = GWAS$GENE_ANNOT

if (sum(!is.na(GWAS$GENE_ANNOT) & GWAS$P < 10^-yend) > 0) {
  GWAS[!is.na(GWAS$GENE_ANNOT) & GWAS$P < 10^-yend,]$P_exceeding = 10^-yend
  GWAS[!is.na(GWAS$GENE_ANNOT) & GWAS$P < 10^-yend,]$BPcum_exceeding = GWAS[!is.na(GWAS$GENE_ANNOT) & GWAS$P < 10^-yend,]$BPcum + 60000000
  GWAS[!is.na(GWAS$GENE_ANNOT) & GWAS$P < 10^-yend,]$GENE_ANNOT_exceeding = paste0(GWAS[!is.na(GWAS$GENE_ANNOT) & GWAS$P < 10^-yend,]$GENE_ANNOT, '\n(p = ', sprintf("%0.e",GWAS[!is.na(GWAS$GENE_ANNOT) & GWAS$P < 10^-yend,]$P), ')')
}

# create labels (prevent overplotting)
label = c(1:18,'', 20,'',22,'X','Y MT','')
if (length(label) != length(axis.set$center)) {
  label = as.character(factor(axis.set$CHR, levels = c(1:25), labels = c(1:22,'X','Y','MT')))
  label[which(label %in% c('19','21'))] = ''
}

# create manhattan plot
message(paste0('Creating manhattan plot.'))
set.seed(8245)

manhplot = ggplot(data = subset(GWAS, P >= 10^-yend)) +
  geom_point(aes(x=BPcum, y=-log10(P), color=as.factor(CHR)), alpha = 1, size = 0.8, stroke = 0, shape = 16) +
  geom_point(data=subset(GWAS, !is.na(GENE_ANNOT) & P < sig & P >= 10^-yend), aes(x=BPcum, y=-log10(P)), color="black", shape=5, size=3) +
  scale_color_manual(values = rep(c("#282873", "#6e90ca"), nCHR)) +
  scale_x_continuous(expand = expansion(mult = c(0.03,0.03), add = c(0,0)), label = label, breaks = axis.set$center) + # label = axis.set$CHR % as.character(as.factor(axis.set$CHR, levels = c(1:25), labels = c(1:22,'X','Y','M'))) % label = c(1:22, 'X', 'Y MT', '') label = c(1:18,'', 20, '', 22, 'X', 'Y MT', '')
  scale_y_continuous(expand = expansion(mult = c(0,0), add = c(0,0.5)), limits = c(0,yend), breaks = seq(0,yend,ysteps)) +
  geom_hline(yintercept = -log10(sig), color = "black", linetype = "solid", size = 0.25) + 
  labs(x = "Chromosome", y = expression("-log"[10]*"("*italic(p)*")")) + 
  geom_segment(aes(x=min(axis.set$center),xend=max(axis.set$center),y=-Inf,yend=-Inf), colour = "black", size = 0.25)+
  geom_segment(aes(y=0,yend=yend,x=-Inf,xend=-Inf), colour = "black", size = 0.25) +
  geom_text_repel(
    data          = subset(GWAS, !is.na(GENE_ANNOT) & P < annotationThresh & P >= 10^-yend),
    aes(x=BPcum, y=-log10(P), label=GENE_ANNOT),
    size          = 3.5,# 1.5
    segment.size  = 0.2 , # 0.2
    segment.color = "grey30",
    direction = "both", # direction = "y"
    nudge_y = 1.7
  ) +
  geom_text_repel(
    data          = subset(GWAS, !is.na(GENE_ANNOT) & P < annotationThresh & P < 10^-yend),
    aes(x=BPcum_exceeding, y=-log10(P_exceeding), label=GENE_ANNOT_exceeding),
    size          = 3.5,# 1.5
    segment.size  = 0.2 , # 0.2
    segment.color = "grey30",
    direction = "y", # direction = "y"
    nudge_x = 80000000,
    nudge_y = -5,
    arrow = arrow(length = unit(0.02, 'npc'))
  ) +
  theme_bw() +
  theme( 
    plot.margin=unit(c(0.25,0.75,0,0),"cm"),
    legend.position = "none",
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "black", size = 0.25),
    axis.ticks.length=unit(.15, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 5, b = 0, l = 5)),
    axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.x = element_text(angle = 0, size = 12, vjust = 0.5, margin = margin(t = 3, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(angle = 0, size = 12, vjust = 0.5, margin = margin(t = 0, r = 3, b = 0, l = 0))
  )

# save plot
message(paste0('Saving manhattan plot.'))
ggsave(paste0(targetDir, '/manhattan.png'), manhplot, width = width, height = height, units = "in", dpi = 300)
system(paste0('chmod -R 770 ', targetDir, '/manhattan.png'))

# save log file
sink(paste0(targetDir, '/manhattan.log')) 
sink(stdout(), type = "message")
message(logInfo)
sink()
system(paste0('chmod -R 770 ', targetDir, '/manhattan.log'))
message(paste0('--- Manhatten plot finished ---\n'))

