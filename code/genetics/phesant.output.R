#!/usr/bin/env Rscript

# ==============================================
# === create phesant output table and figure ===
# ==============================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=10) {
  stop(paste0('expected 10 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
trait = args[1] # trait="gap_gm"
phesantResults= args[2] # phesantResults="results/gap_gm/phesant/phesant.output/results-combined.txt
targetDir = args[3] # targetDir="results/gap_gm/phesant/"
imaging = args[4] # TRUE or FALSE
multipleTesting = args[5] # 'fdr' or 'bonferroni' or 'both'
ylim = as.numeric(args[6]) # upper y axis limit
ybreaks = as.numeric(args[7]) # y axis breaks
width = as.numeric(args[8]) # plot width (8.7 inch)
height = as.numeric(args[9]) # plot height (5.3 inch)
repel_nudge_y = as.numeric(args[10]) # plot height (5.3 inch)

message(paste0('\n--- Settings ---',
        '\ntrait: ', trait,
        '\nphesantResults: ', phesantResults,
        '\ntargetDir: ', targetDir,
        '\nmultipleTesting: ', multipleTesting,
        '\nimaging: ', imaging,
        '\nylim: ', ylim,
        '\nybreaks: ', ybreaks,
        '\nwidth: ', width,
        '\nheight: ', height,
        '\nrepel_nudge_y: ', repel_nudge_y,'\n'))

# attach packages to current R session
for (pkg in c('data.table', 'dplyr', 'ggplot2', 'plotly', 'ggrepel')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# get phesant output 'combined results'
df = data.frame(fread(phesantResults, sep = '\t', header = T, quote = ""))

# remove imaging variables and those without pvalue
if (imaging == FALSE) {
  df = df[df$Cat2_Title!= "Imaging",]
}
df = df[!is.na(df$pvalue),]

# remove quotes in 'description'
df$description = gsub('"' , '', df$description, fixed = TRUE)

# calculate z-value (and se)
df$z = qnorm(df$p/2)*-sign(df$beta)
df$z[df$pvalue==0] = df[df$pvalue==0,] %>% 
                     mutate(z = beta/(abs(lower-upper)/(qnorm(1-(0.05/2))*2))) %>%
                     pull(., var = 'z')
df$se = abs(df$lower-df$upper)/(qnorm(1-(0.05/2))*2)
df$se[is.na(df$se)] = df$beta[is.na(df$se)]/df$z[is.na(df$se)]
df$se[df$beta==-999]=-999

# determine total sample size per phenotype
df$ntotal = df$n %>%
  gsub(pattern = ".*[(]", replacement = "") %>% 
  gsub(pattern = "[)].*", replacement = "") %>%
  gsub(pattern = ".*[/]", replacement = "") %>%
  as.numeric()

# calculate rho and absolute value of rho
df = df %>%
  mutate(rho = sqrt(z^2/(z^2 + (ntotal-2-0)))*sign(beta)) %>%
  mutate(rho_abs = abs(rho))

# set pvalue to 3e-308 if 0
df$pvalue[df$pvalue == 0] = 3e-308

# make custom categories
df$custom_category = df$Cat3_Title
df$custom_category[df$Cat2_Title=='Imaging'] = 'Imaging'
df$custom_category[df$Cat2_Title=='Local environment'] = 'Local environment'
df$custom_category[df$Cat2_Title=='Physical activity measurement'] = 'Physical activity'
df$custom_category[df$Cat2_Title=='Cognitive function online'] = 'Cognitive function'
df$custom_category[df$Cat2_Title=='Diet by 24-hour recall'] = 'Diet by 24-hour recall'
df$custom_category[df$Cat2_Title=='Digestive health'] = 'Digestive health'
df$custom_category[df$Cat2_Title=='Mental health'] = 'Mental health'
df$custom_category[df$Cat2_Title=='Work environment'] = 'Work environment'
df$custom_category[df$Cat2_Title=='Baseline characteristics'] = 'Sociodemographics'
df$custom_category[df$Cat2_Title=='Cognitive function'] = 'Cognitive function'
df$custom_category[df$Cat2_Title=='Physical measures'] = 'Physical measures'
df$custom_category[df$Cat3_Title=='Physical activity measurement'] = 'Physical activity'
df$custom_category[df$Cat3_Title=='Baseline characteristics'] = 'Sociodemographics'
df$custom_category[df$Cat3_Title=='Reception'] = 'Sociodemographics'
df$custom_category[df$Cat3_Title=='Employment'] = 'Sociodemographics'
df$custom_category[df$Cat3_Title=='Summary Administration'] = 'Hospital Inpatient - Administration'
df$custom_category[df$Cat3_Title=='Summary Operations'] = 'Operations'
df$custom_category[df$Cat3_Title=='Summary Maternity'] = 'Maternity and sex-specific factors'
df$custom_category[df$Cat3_Title=='Summary Psychiatric'] = 'Mental health'
df$custom_category[df$Cat3_Title=='Sex-specific factors'] = 'Maternity and sex-specific factors'
df$custom_category[df$Cat3_Title=='Health and medical history'] = 'Medical history and conditions'
df$custom_category[df$Cat3_Title=='Medical conditions'] = 'Medical history and conditions'
df$custom_category[df$Cat3_Title=='Asthma outcomes'] = 'Medical history and conditions'
df$custom_category[df$Cat3_Title=='Myocardial infarction outcomes'] = 'Medical history and conditions'
df$custom_category[df$Cat3_Title=='Psychosocial factors'] = 'Mental health'
df$custom_category[df$Cat3_Title=='Family history'] = 'Family history and early life factors'
df$custom_category[df$Cat3_Title=='Early life factors'] = 'Family history and early life factors'

# calculate fdr
df$fdr = p.adjust(df$pvalue, method = 'BH')

# ----------------------------------
# --- draw phewas manhattan plot ---
# ----------------------------------

# create new data frame for plot
dfplot = df
dfplot = dfplot[order(dfplot$custom_category, dfplot$varName),]
ncategories = length(unique(dfplot$custom_category))

# set x axis positions for trait-phenotype associations
k = 0
for (category in unique(dfplot$custom_category)) {
    dfplot$pos[dfplot$custom_category == category] = ppoints(sum(dfplot$custom_category == category)) + k
    k = k +1
}

# set position for x axis breaks (centered)
axis.set = dfplot %>% 
  group_by(custom_category) %>% 
  summarize(center = mean(pos))

# set ylim significance threshold
if (multipleTesting == 'fdr') {
  sig = max(dfplot$pvalue[dfplot$fdr < 0.05])
} else if (multipleTesting == 'bonferroni' | multipleTesting == 'both') {
  sig = 0.05/nrow(dfplot)
}

# add variables for plotly tooltip
dfplot$variable = dfplot$description
dfplot$`field id` = dfplot$varName
dfplot$p = sprintf("%.6f", as.numeric(dfplot$pval))
dfplot$p[as.numeric(dfplot$pvalue) < 0.000001] = sprintf("%.2g", as.numeric(dfplot$pvalue[as.numeric(dfplot$pvalue) < 0.000001]))
dfplot$r = sprintf("%.6f", as.numeric(dfplot$rho))

# set pvalue of< ylim to ylim
dfplot$pvalue[dfplot$pvalue < 10^-ylim] = 10^-ylim

# get top hits per category for annotation
top = data.table(dfplot[order(dfplot$pvalue),], key = 'custom_category')
top = top[, head(.SD, 1), by=custom_category]
top = top[top$pvalue < sig,]

# set plot colors 
gg_color_hue = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
desat = function(cols, sat=0.5) {
  X = diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
  hsv(X[1,], X[2,], X[3,])
}

set.seed(1234)
colcolors = gg_color_hue(ncategories)
colcolors = colcolors[sample(ncategories, replace = F)]
fillcolors = desat(colcolors, sat = 0.25)

# create manhattan plot
phewas_plot = ggplot2::ggplot(data = dfplot[dfplot$pvalue >= sig,], aes(textvariable = variable, textfield = `field id`, textn = n, textr = r, textp = p)) +
  geom_point(aes(x=pos, y=-log10(pvalue), fill=as.factor(custom_category)), alpha = 1, size = 1.5, stroke = 0, shape = 21) +
  geom_point(data = dfplot[dfplot$pvalue < sig,], aes(x=pos, y=-log10(pvalue), color=as.factor(custom_category)), fill = NA,  alpha = 1, size = 1.5, stroke = 0, shape = 16) +
  scale_fill_manual(values = fillcolors) +
  scale_color_manual(values = colcolors[unique(dfplot$custom_category) %in% dfplot$custom_category[dfplot$pvalue < sig]]) +
  scale_x_continuous(expand = c(0,0), limits = c(-0.5,max(axis.set$center)+1), label = unique(dfplot$custom_category), breaks = axis.set$center) + # label = axis.set$CHR % label = c(1:22, 'X', 'Y MT', '') label = c(1:18,'', 20, '', 22, 'X', 'Y MT', '')
  scale_y_continuous(expand = c(0,0), limits = c(0,ylim), breaks = seq(0,ylim,ybreaks)) +
  geom_hline(yintercept = -log10(sig), color = "black", linetype = "solid", size = 0.25) + 
  geom_segment(aes(x=min(axis.set$center),xend=max(axis.set$center),y=0,yend=0), colour = "black", size = 0.25) +
  geom_segment(aes(y=0,yend=ylim,x=-0.5,xend=-0.5), size = 0.25) +
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "black", size = 0.25),
    axis.ticks.length=unit(.15, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle = -45, size = 9, vjust = 0.5, hjust = 0, margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(angle = 0, size = 9, vjust = 0.5, margin = margin(t = 0, r = 3, b = 0, l = 0)),
    axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 5, b = 0, l = 5)),
    axis.title.x = element_text(size = 14, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    plot.margin=unit(c(0.25,1.0,-0.5,0),"cm"))

# add fdr line
if (multipleTesting == 'both') {
  phewas_plot = phewas_plot + geom_hline(yintercept = -log10(max(dfplot$pvalue[dfplot$fdr < 0.05])), color = "black", linetype = "dashed", size = 0.25)
}

# draw ggplot
phewas_ggplot = phewas_plot +
    labs(x = "", y = expression("-log"[10]*"("*italic(p)*")")) + 
    geom_text_repel(data = top, aes(x=pos, y=-log10(pvalue), label=description),
                    max.time = 5,
                    max.iter = 1000000,
                    box.padding = 0.5,
                    nudge_x = 0.5,
                    nudge_y = repel_nudge_y, # 10
                    direction = 'both',
                    size = 2.5,
                    angle = 0,
                    hjust = 0,
                    segment.size = 0.25,
                    segment.color = "grey10",
                    segment.square  = T,
                    segment.squareShape  = 1,
                    segment.angle     = -90,
                    segment.curvature = 1,
                    segment.ncp = 1,
                    segment.inflect = FALSE)

# draw plotly
phewas_plotly = plotly::ggplotly(phewas_plot + labs(x = ""), tooltip = c("textvariable", "textfield", "textn", "textr", "textp"), dynamicTicks = FALSE, width = width*100, height = height*100)
phewas_plotly = phewas_plotly %>% config(displayModeBar = F) %>% layout(xaxis=list(fixedrange=TRUE)) %>% layout(yaxis=list(fixedrange=TRUE))

# ==========================================
# === saves figure and output data frame ===
# ==========================================

# save ggplot as png
message('(1/3) Saving ggplot as png.')
system(paste0('mkdir -p ', targetDir))
png(width = width, height = height, units = "in", res = 600, filename = paste0(targetDir, '/phesant.png'))
phewas_ggplot
invisible(dev.off())
system(paste0('chmod 770 ', targetDir,'/phesant.png'))

# save plotly as interactive html
# change height and width: phewas_plotly = phewas_plotly %>% layout(height = 800, width = 1200)
message('(2/3) Saving interactive plotly as html.')
htmlwidgets::saveWidget(phewas_plotly, paste0(targetDir, '/phesant.html'), title = paste0('PheWAS ', trait), selfcontained = TRUE)
system(paste0('chmod 770 ', targetDir,'/phesant.html'))
system(paste0('rm -rf ',targetDir,'/phesant_files'))
      
# export result table
message('(3/3) Saving plain text file with results summary.')
output = df[,c('description', 'varName', 'custom_category', 'Path', 'varType', 'resType', 'n',  'beta', 'lower', 'upper', 'se', 'pvalue', 'fdr', 'z', 'ntotal', 'rho')]
write.table(output, file = paste0(targetDir,'/phesant.summary.txt'), sep = '\t', quote = F, row.names = F, col.names = T)
system(paste0('chmod 770 ', targetDir,'/phesant.summary.txt'))
message('-- Creating PHESANT output table and figure completed. ---')
