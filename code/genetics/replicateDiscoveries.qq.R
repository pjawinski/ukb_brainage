#!/usr/bin/env Rscript

# ====================
# === draw qq plot ===
# ====================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=14) {
  stop(paste0('expected 14 arguments, but ', length(args), ' arguments provided.'), call.=FALSE)
}

# get arguments from command line
inputFile = args[1] # inputFile="results/combined/replicateDiscoveries.txt"
outFile = args[2] # outFile="results/combined/replicateDiscoveries.gm.png"
pCol = args[3] # pCol="REPLIC_P"
criterionCols = args[4] # criterionCols="DISCOV_TRAITNAME"
criterionVals = args[5] # criterionVals="grey matter"
duplicatedCol = args[6] # duplicatedCol="DISCOV_LOCUS_COUNT"
prune = args[7] # prune = TRUE
drawCI = args[8] # drawCI = FALSE
xend = as.numeric(args[9]) # xend = 8
xsteps = as.numeric(args[10]) # xsteps = 2
yend = as.numeric(args[11]) # yend = 12
ysteps = as.numeric(args[12]) # ysteps = 2
width = as.numeric(args[13]) # width = 3
height = as.numeric(args[14]) # height = 4

logInfo = paste0('\n--- qq-plot settings ---',
               '\ninputFile: ', inputFile,
               '\noutFile: ', outFile,
               '\npCol: ', pCol,
               '\nncriterionCols: ', criterionCols,
               '\nncriterionVals: ', criterionVals,
               '\nnduplicateCol: ', duplicatedCol,
               '\nprune: ', prune,
               '\ndrawCI: ', drawCI,
               '\nxend: ', xend,
               '\nxsteps: ', xsteps,
               '\nyend: ', yend,
               '\nysteps: ', ysteps,
               '\nwidth: ', width,
               '\nheight: ', height,'\n')
message(logInfo)

# attach required packages
for (pkg in c('data.table','dplyr','ggplot2','stringr')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# transform variables
criterionCols = str_split(criterionCols, ',')[[1]]
criterionVals = str_split(criterionVals, ',')[[1]]

# import data
message(paste0('[1/3] Importing data.'))
if (str_sub(inputFile,-3,-1) == '.gz') {
  df = data.frame(fread(cmd=paste0("gzip -dc ", inputFile), tmpdir = getwd(), header=T, stringsAsFactors=FALSE))
} else {
  df = read.delim(inputFile, header = T, sep = '\t')
}

# extract data
if (criterionCols!='') {
  for (i in 1:length(criterionCols)) {
    df = df[df[,criterionCols[i]]==criterionVals[i],]
  }
}
if (duplicatedCol!='') {
  df = df[!duplicated(df[,duplicatedCol]),]
}

# create qq-plot function
# credits to Kamil Slowikowski (https://slowkow.com/notes/ggplot2-qqplot/)
gg_qqplot = function(pvals, ci = 0.90, drawCI = TRUE, truncateLim = Inf, pointshape = 16, pointsize = 2, pointcolor = "#0072BD", pointalpha = 1, prune = FALSE) {
  n  = length(pvals)
  df = data.frame(
    observed = -log10(sort(pvals)),
    expected = -log10(ppoints(n))
  )
  
  # Truncate
  df$observed[df$observed > truncateLim] = truncateLim

  # set axis labels
  log10Pe = expression("Expected -log"[10]*"("*italic(p)*")")
  log10Po = expression("Observed -log"[10]*"("*italic(p)*")")

  # draw confidence intervals, or don't.
  if (drawCI == TRUE) {
  df$clower = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1))
  df$cupper = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  }

  # remove overlapping points (pvals > 0.01)
  if (prune == TRUE) {
  df$expected[df$expected<2] = round(df$expected[df$expected<2],2)
  df = df[-which(df$expected<2 & duplicated(df$expected)),]
  }

  # draw qq-plot
  if (drawCI == TRUE) {
  qqplot = ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1)
  } else {
  qqplot = ggplot(df)
  }
  qqplot +
    geom_point(aes(expected, observed), shape = pointshape, size = pointsize, colour = pointcolor, alpha = pointalpha) + #, shape = 1, size = 3
    geom_abline(intercept = 0, slope = 1, alpha = 1, color = "black", linetype="solid", size = 0.25) + # geom_segment(aes(x = 0, xend = max(expected), y = 0, yend = max(expected)), alpha = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}

# calculate lambdagc
lambdaGC = qchisq(median(df[,pCol]), 1, lower.tail = F)/qchisq(0.5, 1, lower.tail = F)

# create plot
message(paste0('[2/3] Creating qq plot.'))
qqplot = gg_qqplot(df[,pCol], prune = prune, drawCI = drawCI, truncateLim = yend) +
    geom_segment(aes(x=0,xend=xend,y=-Inf,yend=-Inf), colour = "black", size = 0.25) +
    geom_segment(aes(y=0,yend=yend,x=-Inf,xend=-Inf), colour = "black", size = 0.25) +
    theme_bw() +
        scale_x_continuous(expand = expansion(mult = c(0.03,0), add = c(0,0)), limits = c(0,xend), breaks = seq(0,xend,xsteps)) +
        scale_y_continuous(expand = expansion(mult = c(0,0), add = c(0,0.5)), limits = c(0,yend), breaks = seq(ysteps,yend,ysteps)) +
        theme(plot.margin=unit(c(0.25,0.75,0,0),"cm"),
              panel.border = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_line(colour = "black", size = 0.25),
              axis.ticks.length=unit(.15, "cm"),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 5, b = 0, l = 5)),
              axis.title.x = element_text(size = 14, margin = margin(t = 10, r = 0, b = 0, l = 0)),
              axis.text.x = element_text(angle = 0, size = 12, vjust = 0.5, margin = margin(t = 3, r = 0, b = 0, l = 0)),
              axis.text.y = element_text(angle = 0, size = 12, vjust = 0.5, margin = margin(t = 0, r = 3, b = 0, l = 0)),
              line = element_line(size=0.25)) #+
    #annotate("text", x = -Inf, y = Inf, hjust = -0.25, vjust = 2.5, label = sprintf("lambda[GC] == '%0.3f'", lambdaGC), size = 4, parse = TRUE)

# save plot
message(paste0('[3/3] Saving qq plot.'))
ggsave(outFile, qqplot, width = width, height = height, units = "in", dpi = 300)
system(sprintf('chmod -R 770 %s', outFile))
message(paste0('--- qq-plot completed ---\n'))


