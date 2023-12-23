#!/usr/bin/env Rscript

# ===================================
# === draw quantile-quantile plot ===
# ===================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=14) {
  stop(paste0('expected 14 arguments, but ', length(args), ' arguments provided.'), call.=FALSE)
}

# get arguments from command line
trait = args[1] # trait="gap_gm"
targetDir = args[2] # targetDir="results/gap_gm/replicate/EUR/"
sumstats = args[3] # sumstats="results/gap_gm/replicate/EUR/sumstats.txt.gz"
pCol = args[4] # pCol="Pvalue" 
nCriterion = args[5] # nCriterion = 1 | LDSC-like method to exclude snps with N < quantile(N, 0.9) / 1.5 
prune = args[6] # prune = TRUE
drawCI = args[7] # drawCI = FALSE
drawLambda = args[8] # drawLambda = FALSE
xend = as.numeric(args[9]) # xend = 8
xsteps = as.numeric(args[10]) # xsteps = 2
yend = as.numeric(args[11]) # yend = 35
ysteps = as.numeric(args[12]) # ysteps = 5
width = as.numeric(args[13]) # width = 3
height = as.numeric(args[14]) # height = 4

logInfo = paste0('\n--- qq-plot settings ---',
               '\ntrait: ', trait,
               '\ntargetDir: ', targetDir,
               '\nsumstats: ', sumstats,
               '\npCol: ', pCol,
               '\nnCriterion: ', nCriterion,
               '\nprune: ', prune,
               '\ndrawCI: ', drawCI,
               '\ndrawLambda: ', drawLambda,
               '\nxend: ', xend,
               '\nxsteps: ', xsteps,
               '\nyend: ', yend,
               '\nysteps: ', ysteps,
               '\nwidth: ', width,
               '\nheight: ', height,'\n')
message(logInfo)

# attach required packages
for (pkg in c('data.table','dplyr','ggplot2')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# create folder
system(paste0('mkdir -p ', targetDir))

# import data
message(paste0('[1/3] Importing data.'))
if (stringr::str_sub(sumstats,-3,-1) == '.gz') {
  GWAS = data.frame(fread(cmd=sprintf("gzip -dc %s", sumstats), tmpdir = getwd(), header=T, stringsAsFactors=FALSE))
} else {
  GWAS = data.frame(fread(file = sumstats, tmpdir = getwd(), header=T, stringsAsFactors=FALSE))
}

# exclude snps that do not meet N criterion
if (nCriterion == TRUE) {
  message('Excluding variants that do not meet N criterion.')
  maxN = max(GWAS$N)
  minN = quantile(GWAS$N, probs = 0.9)[[1]] / 1.5
  nExcl = sum(GWAS$N < minN)
  nIncl = sum(GWAS$N >= minN)
  GWAS = GWAS[GWAS$N >= minN,]
  logInfo = paste0(logInfo,
    '\n--- calculations ---',
    '\nmaxN: ', maxN,
    '\nminN: ', minN,
    '\nnExcl: ', nExcl,
    '\nnIncl: ', nIncl)
}

# create qq-plot function
# credits to Kamil Slowikowski (https://slowkow.com/notes/ggplot2-qqplot/)
gg_qqplot = function(pvals, ci = 0.90, drawCI = TRUE, truncateLim = Inf, pointshape = 16, pointsize = 1, pointcolor = "#0072BD", pointalpha = 1, prune = FALSE) {
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
lambdaGC = qchisq(median(GWAS[,pCol]), 1, lower.tail = F)/qchisq(0.5, 1, lower.tail = F)

# create plot
message(paste0('[2/3] Creating qq plot.'))
qqplot = gg_qqplot(GWAS[,pCol], prune = prune, drawCI = drawCI, truncateLim = yend) +
    geom_segment(aes(x=0,xend=xend,y=-Inf,yend=-Inf), colour = "black", size = 0.25) +
    geom_segment(aes(y=0,yend=yend,x=-Inf,xend=-Inf), colour = "black", size = 0.25) +
    theme_bw() +
        scale_x_continuous(expand = expansion(mult = c(0.03,0), add = c(0,0)), limits = c(0,xend), breaks = seq(0,xend,xsteps)) +
        scale_y_continuous(expand = expansion(mult = c(0,0), add = c(0,0.5)), limits = c(0,yend), breaks = seq(0,yend,ysteps)) +
        theme(plot.margin=unit(c(0.25,0.75,0,0),"cm"),
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
              axis.text.y = element_text(angle = 0, size = 12, vjust = 0.5, margin = margin(t = 0, r = 3, b = 0, l = 0)))

# add Lambda
if (drawLambda == TRUE) {
qqplot = qqplot +
    annotate("text", x = -Inf, y = Inf, hjust = -0.25, vjust = 2.5, label = sprintf("lambda[GC] == '%0.3f'", lambdaGC), size = 4, parse = TRUE)
}

# save plot
message(paste0('[3/3] Saving qq plot.'))
ggsave(paste0(targetDir, '/qqplot.png'), qqplot, width = width, height = height, units = "in", dpi = 300)
system(paste0('chmod -R 770 ', targetDir, '/qqplot.png'))

# save log file
sink(paste0(targetDir, '/qqplot.log')) 
sink(stdout(), type = "message")
message(logInfo)
sink()
system(paste0('chmod -R 770 ', targetDir, '/qqplot.log'))
message(paste0('--- qq-plot completed ---\n'))


