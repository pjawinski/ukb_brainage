#!/usr/bin/env Rscript

# =============================================
# === plot partitioned heritability results ===
# =============================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 6) {
  stop(paste0('expected 6 arguments, but ', length(args), ' arguments provided.'), call.=FALSE)
}

# get arguments from command line
ldscFile = args[1] # ldscFile="results/gap_gm/ldsc/ldsc.partitioned.results"
ymin = as.numeric(args[2]) # ymin=-10
ymax = as.numeric(args[3]) # ymax=25
ysteps = as.numeric(args[4]) # ysteps=10
width = as.numeric(args[5]) # width=8.1
height = as.numeric(args[6]) # height=3.7
message(paste0('\n--- Plot partitioned heritability results ---',
               '\nldscFile: ', ldscFile,
               '\nymin: ', ymin,
               '\nymax: ', ymax,
               '\nysteps: ', ysteps,
               '\nwidth: ', width,
               '\nheight: ', height, '\n'))

# load required pages
for (pkg in c('dplyr','ggplot2')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# load data
message(sprintf('(1/4) Loading %s...',ldscFile))
df = read.delim(ldscFile, sep = '\t', head = TRUE)

# define relevant categories (in accordance with Gazal et al. 2018, suppl. Table S1)
main = data.frame(matrix(ncol = 2, byrow = TRUE, data = c(
  'Coding_UCSCL2_0', 'Coding', 
  'Intron_UCSCL2_0', 'Intron', 
  'UTR_3_UCSCL2_0', '3\'-UTR', 
  'UTR_5_UCSCL2_0', '5\'-UTR',
  'Promoter_UCSCL2_0', 'Promoter', 
  'PromoterFlanking_HoffmanL2_0', 'Promoter Flanking', 
  'Enhancer_HoffmanL2_0', 'Enhancer', 
  'WeakEnhancer_HoffmanL2_0', 'Weak Enhancer', 
  'SuperEnhancer_HniszL2_0', 'Super Enhancer (Hnisz)', 
  'Enhancer_AnderssonL2_0', 'FANTOM5 Enhancer',
  'Transcr_HoffmanL2_0', 'Transcribed', 
  'TSS_HoffmanL2_0', 'Transcription Start Site', 
  'CTCF_HoffmanL2_0', 'CCCTC Binding Factor', 
  'Repressed_HoffmanL2_0', 'Repressed',
  'DHS_TrynkaL2_0', 'DNase I hypersensitive sites', 
  'FetalDHS_TrynkaL2_0', 'Fetal DNase I hypersensitive sites', 
  'BivFlnkL2_0', 'Flanking Bivalent TSS/Enhancer',
  'H3K27ac_HniszL2_0', 'H3K27ac (Hnisz)', 
  'H3K27ac_PGC2L2_0', 'H3K27ac (PGC2)', 
  'H3K4me1_TrynkaL2_0', 'H3K4me1', 
  'H3K4me3_TrynkaL2_0', 'H3K4me3', 
  'H3K9ac_TrynkaL2_0', 'H3K9ac', 
  'TFBS_ENCODEL2_0', 'Transcription Factor Binding Site', 
  'DGF_ENCODEL2_0', 'Digital Genomic Footprint', 
  'non_synonymousL2_0', 'Non-synonymous', 
  'synonymousL2_0', 'Synonymous', 
  'GERP.RSsup4L2_0', 'Conserved (GERP RS >= 4)', 
  'Conserved_LindbladTohL2_0', 'Conserved (Lindblad-Toh)', 
  'Conserved_Vertebrate_phastCons46wayL2_0', 'Conserved (Vertebrate)', 
  'Conserved_Mammal_phastCons46wayL2_0', 'Conserved (Mammal)', 
  'Conserved_Primate_phastCons46wayL2_0', 'Conserved (Primate)')))
names(main) = c('Category', 'annotation')

# merge main with df set annotation variable as factor
main = left_join(main, df, 'Category')
main$annotation = sprintf('%s: %0.1f', main$annotation, main$Prop._SNPs*100)
main = main[order(main$Prop._SNPs),]
main$annotation = factor(main$annotation, levels = main$annotation)
main$Enrichment_FDR = p.adjust(main$Enrichment_p, method = 'BH')
main$significance = 0
main$significance[main$Enrichment_FDR < 0.05] = 1
main$significance = factor(main$significance, levels = c(0,1), labels = c('n.s.', 'FDR < 0.05'))

# plot enrichment
message('(2/4) Make partitioned heritability plot...')

enrich = ggplot(main, aes(x = annotation, y = Enrichment, color = significance, shape = significance)) + # abs(tvalue)
  geom_errorbar(aes(x = annotation, ymin = Enrichment - Enrichment_std_error, ymax = Enrichment + Enrichment_std_error),
                width=0.25, colour="black", alpha=1, size=0.25) +
  geom_hline(yintercept=1, linetype="dashed", color = "black", size = 0.25) +
  geom_point(size = 1.8) +
  scale_color_manual(values=c('black','red')) +
  scale_shape_manual(values=c(16,5)) +
  labs(color = "Significance level", shape = "Significance level") +
  coord_cartesian(ylim = c(ymin,ymax), expand = TRUE) + # gsub("\\.", " ",
  scale_y_continuous(breaks = seq(ymin,ymax,ysteps)) +
  ylab(expression("Enrichment of "*italic(h)^2)) + # ylab(expression(atop("Enrichment", "proportion of "*italic(h)^2*" / proportion of SNPs"))) +
  xlab("Functional annotation (ordered by proportion of SNPs)") +
  theme_bw() +
  theme(plot.margin=margin(15,15,5,15),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(hjust = 1, size = 8),
        axis.title = element_text(size = 10),
        axis.title.x = element_blank(), # element_text(size=9, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=11, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=9),
        axis.ticks = element_line(size = 0.25),
        #axis.ticks.length=unit(.1, "cm"),
        axis.line = element_line(size = 0.1),
        legend.title = element_text(size=8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        legend.position = c(0.92, 0.85))

# save plot
message(sprintf('(3/4) Saving  %s.png...',ldscFile))
png(filename = paste0(ldscFile, '.png'), width=8.1, height=3.7, units = "in", res = 600)
enrich
invisible(dev.off())
system(paste0('chmod 770 ', ldscFile, '.png'))

# save results
message(sprintf('(4/4) Saving  %s.summary.txt...',ldscFile))
output = data.frame(annotation = sub(':.*','', main$annotation),
    Prop._SNPs = main$Prop._SNPs,
    Prop._h2 = main$Prop._h2,
    Enrichment = main$Enrichment,
    Enrichment_p = main$Enrichment_p,
    Enrichment_FDR = main$Enrichment_FDR)
write.table(output, file = paste0(ldscFile,'.summary.txt'), sep = '\t', col.names = T, row.names = F, quote = F)
message('--- Completed: Plot partitioned heritability results ---')

