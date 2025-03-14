#!/usr/bin/env Rscript

# =========================================================
# ===  Compare Z values derived from different analyses ===
# =========================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=20) {
  stop(paste0('expected 20 arguments, but ', length(args), ' arguments provided.'), call.=FALSE)
}

# get arguments from command line
outFile = args[1] # outFile="results/combined/discovery.sex.zscores"
input1 = args[2] # input1="/Users/philippe/Desktop/research/ukb_brainage/results/combined/discovery.sex.phewas.txt"
input2 = args[3] # input2="/Users/philippe/Desktop/research/ukb_brainage/results/combined/discovery.sex.phewas.txt"
input1Col = args[4] # input1Col="a_gap_wm_pvalue"
input2Col = args[5] # input2Col="b_gap_wm_pvalue"
input1SignCol = args[6] # input1SignCol="a_gap_wm_beta"
input2SignCol = args[7] # input2SignCol="b_gap_wm_beta"
matchCol = args[8] # matchCol="varName"
conversion = args[9] # conversion="none" # either 'none' or 'p-to-z' or 'r-to-z'
PDadjust = as.numeric(args[10]) # PDadjust=0.375
PDsize = as.numeric(args[11]) # PDsize=0.5
xlabel = args[12] # xlabel="Males (z-scores)"
ylabel = args[13] # ylabel="Females (z-scores)"
xmin = as.numeric(args[14]) # xmin = -13
xmax = as.numeric(args[15]) # xmax = 13
ymin = as.numeric(args[16]) # ymin = -12
ymax = as.numeric(args[17]) # ymax = 12
width = as.numeric(args[18]) # width = 2.9
height = as.numeric(args[19]) # height = 2.7
preview = args[20] # preview = FALSE

message('\n--- Compare Z scores ---',
                 '\noutFile: ', outFile,
                 '\ninput1: ', input1,
                 '\ninput2: ', input2,
                 '\ninput1Col: ', input1Col,
                 '\ninput2Col: ', input2Col,
                 '\ninput1SignCol: ', input1SignCol,
                 '\ninput2SignCol: ', input2SignCol,
                 '\nmatchCol: ', matchCol,
                 '\nconversion: ', conversion,
                 '\nPDadjust: ', PDadjust,
                 '\nPDsize: ', PDsize,
                 '\nxlabel: ', xlabel,
                 '\nylabel: ', ylabel,
                 '\nxmin: ', xmin,
                 '\nxmax: ', xmax,
                 '\nymin: ', ymin,
                 '\nymax: ', ymax,
                 '\nwidth: ', width,
                 '\nheight: ', height,
                 '\npreview: ', preview,'\n')

# attach required packages and create target directory
for (pkg in c('data.table','dplyr','ggExtra','ggplot2','ggpointdensity','stringr','viridis')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# import data
message('Importing data.')
if (endsWith(input1,'.gz')) {
  df1 = data.frame(data.table::fread(cmd=paste0('gzip -dc ', input1), tmpdir = getwd(), header=T, stringsAsFactors=FALSE))
} else {
  df1 = data.frame(data.table::fread(input1, tmpdir = getwd(), header=T, stringsAsFactors=FALSE))
}
if (endsWith(input1,'.gz')) {
  df2 = data.frame(data.table::fread(cmd=paste0('gzip -dc ', input2), tmpdir = getwd(), header=T, stringsAsFactors=FALSE))
} else {
  df2 = data.frame(data.table::fread(input2, tmpdir = getwd(), header=T, stringsAsFactors=FALSE))
}

# calculate Z scores
if (conversion=='p-to-z') {
  message('Performing p-to-z conversion.')
  df1$Z = qnorm(df1[,input1Col]/2)*-sign(df1[,input1SignCol])
  df2$Z = qnorm(df2[,input2Col]/2)*-sign(df2[,input2SignCol])
} else if (conversion=='r-to-z') {
  message(paste0('Performing r-to-z conversion.'))
  df1$Z = atanh(df1[,input1Col])
  df2$Z = atanh(df2[,input2Col])
} else {
  message(paste0('No conversion performed.'))
  df1$Z = df1[,input1Col]
  df2$Z = df2[,input2Col]
}

# join data frames
message('Matching by identifier.')
df = inner_join(df1[,c(matchCol,'Z')],df2[,c(matchCol,'Z')], by = matchCol)
message('Associations after merging: ', nrow(df))

# Preview: keep only 1% of variants
if (preview == TRUE) {
  df = df[rep_len(c(TRUE,rep(FALSE,99)), nrow(df)),]
}

# draw plot
message('Plotting variables.')
pl = ggplot(df, aes(x = Z.x, y = Z.y)) +
  geom_pointdensity(adjust = PDadjust, size = PDsize) +
  scale_color_viridis() +
  scale_x_continuous(limits = c(xmin,xmax)) +
  scale_y_continuous(limits = c(ymin,ymax)) +
  xlab(xlabel) +
  ylab(ylabel) +
  theme_bw(base_size=10) +
  theme(plot.margin = unit(c(5, 5, 5, 0), "mm"),
        panel.border = element_rect(linewidth = 0.25),
        legend.position = 'none',
        line = element_line(linewidth = 0.25),
        plot.title = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 1, face = 'bold'),
        axis.title.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = -0.5),
        axis.title.y = element_text(size = 12, angle = 90, hjust = 0.5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.ticks = element_line(linewidth = 0.25))

# save plot
message(sprintf('Saving %s.',outFile))
png(outFile, width = width, height = height, units = 'in', res = 300); pl; invisible(dev.off())
system(sprintf('chmod -R 770 %s', outFile))
message('--- Completed: Compare Z values ---')
