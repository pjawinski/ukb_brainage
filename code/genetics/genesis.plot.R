#!/usr/bin/env Rscript

# =================================================================
# === plot results of genetic effect size distribution analysis ===
# =================================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=7) {
  stop(paste0('expected 7 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}                                                                                                                                                    

# set arguments
traits = args[1] # traits="gap_gm,gap_wm,gap_gwm"
traitLabels = args[2] # traitLabels="Grey_matter,White_matter,Grey_and_white"
genesisFit = args[3] # genesisFit="results/gap_gm/genesis/genesis.fit3.Rds,results/gap_wm/genesis/genesis.fit3.Rds,results/gap_gwm/genesis/genesis.fit3.Rds"
reftraits = args[4] # reftraits="neur,height"
reftraitLabels = args[5] # reftraitLabels="Neuroticism,Height"
refgenesisFit = args[6] # refgenesisFit="data/sumstats/genesis/04_neur_baselmans_2019/genesis.fit3.Rds,data/sumstats/genesis/mr_height_wood_2014/genesis.fit3.Rds"
outFile = args[7] # outFile="results/combined/genesis"

logInfo = paste0('\n--- plot results of genetic effect size distribution analysis ---',
                 '\ntraits: ', traits,
                 '\ntraitLabels: ', traitLabels,
                 '\ngenesisFit: ', genesisFit,
                 '\nreftraits: ', reftraits,
                 '\nreftraitLabels: ', reftraitLabels,
                 '\nrefgenesisFit: ', refgenesisFit,
                 '\noutFile: ', outFile,'\n')
message(logInfo)

# check if package genesis is available and conda environment name is 'genesis'
is_genesis_available <- suppressWarnings(require('GENESIS'))
env = system('echo $CONDA_DEFAULT_ENV', intern = T)
if (!is_genesis_available ) {
  message(' - GENESIS not available. Checking conda environment.')
  if (env != 'genesis') {
    message(' - conda environment name is not genesis. Please load environment genesis.'); stop()
  } else {
    message(' - conda environment name is genesis. Downloading missing package genesis.')
    options(download.file.method = "curl")
    devtools::install_github('yandorazhang/GENESIS', upgrade = 'never')
  }
}

# attach required packages
for (pkg in c('dplyr','GENESIS','ggplot2','MASS','stringr','patchwork')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# transform variables
traits = str_split(traits, ',')[[1]]
traitLabels = str_split(traitLabels, ',')[[1]]
traitLabels = str_replace_all(traitLabels, "_", " ")
genesisFit = str_split(genesisFit, ',')[[1]]
reftraits = str_split(reftraits, ',')[[1]]
reftraitLabels = str_split(reftraitLabels, ',')[[1]]
reftraitLabels = str_replace_all(reftraitLabels, "_", " ")
refgenesisFit = str_split(refgenesisFit, ',')[[1]]

# import trait data
for (i in 1:length(traits)) {
  
  # load data
  message(sprintf('[%d/%d] Loading %s',i,length(traits)+length(reftraits),genesisFit[i]))
  tmp = readRDS(genesisFit[i])
  
  # assign individual  variables
  assign(sprintf('%s.fit',traits[i]),tmp)
}

# import reftrait data
for (i in 1:length(reftraits)) {
  
  # load data
  message(sprintf('[%d/%d] Loading %s',i+length(traits),length(traits)+length(reftraits),refgenesisFit[i]))
  tmp = readRDS(refgenesisFit[i])
  
  # assign individual  variables
  assign(sprintf('%s.fit',reftraits[i]),tmp)
}

# make effect size distribution plot (density plot)
message(' - making density plot of joint effect sizes')

  # compute y dots for traits
  x_seq = seq(-0.02,0.02,length.out = 1000); 
  for (i in 1:length(traits)) {
    
    # compute y dots for density plot
    tmp = get(sprintf('%s.fit',traits[i]))
    tmpParams = tmp$estimates[grep('Parameter',names(tmp$estimates))] %>% unlist() %>% as.vector()
    tmpSeq = apply(matrix(x_seq,ncol=1),1,function(t) dmixssnp(t,tmpParams))
    
    # create data frame for density plot
    if (i==1) { data = data.frame(x_seq,tmpSeq) } else { data = data.frame(cbind(data,tmpSeq)) }
    names(data)[ncol(data)] = traits[i]
  }

  # compute y dots for reference traits
  for (i in 1:length(reftraits)) {
    
    # compute y dots for density plot
    tmp = get(sprintf('%s.fit',reftraits[i]))
    tmpParams = tmp$estimates[grep('Parameter',names(tmp$estimates))] %>% unlist() %>% as.vector()
    tmpSeq = apply(matrix(x_seq,ncol=1),1,function(t) dmixssnp(t,tmpParams))
    
    # create data frame for density plot
    data = data.frame(cbind(data,tmpSeq))
    names(data)[ncol(data)] = reftraits[i]
  }

  # make effect size distribution plot (density plot)
  dfdensity = reshape2::melt(data, id="x_seq", variable.name="trait", value.name="y_seq", na.rm = F)
  dfdensity$trait = factor(dfdensity$trait, levels = c(reftraits,traits), labels = c(reftraitLabels,traitLabels))
  lineTypes = c('solid','dashed','dotted','dotdash','longdash','twodash')
  lineTypes = c(lineTypes[2:(length(reftraits)+1)],lineTypes[1:length(traits)])
  lineColours = c(rep('grey50',length(reftraits)),rep('royalblue4',length(traits))) # lineColours = c(rev(RColorBrewer::brewer.pal(length(c(traits)), 'PuBu')),'grey70','grey70')

  pldensity = ggplot() +
    geom_line(data=dfdensity, aes(x = x_seq, y = y_seq, color = trait, linetype = trait), linewidth = 0.50) +
    scale_colour_manual(values = lineColours) + # c('grey10','grey60','grey90')
    scale_linetype_manual(values = lineTypes) + # c('grey10','grey60','grey90')
    labs(x = "Joint effect size", y = "Probability density") +
    geom_segment(aes(x=-0.02,xend=0.02,y=-Inf,yend=-Inf), , linewidth = 0.25, colour = 'black') + #colour = "grey70"
    geom_segment(aes(y=0,yend=200,x=-Inf,xend=-Inf),  linewidth = 0.25, colour = 'black') + # colour = "grey70",
    theme_light() +
    theme( 
      legend.position = c(0.8,0.6),
      legend.key = element_blank(),
      legend.key.size = unit(0.5, 'cm'),
      legend.background=element_blank(),
      legend.title = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks.length=unit(.2, "cm"),
      axis.ticks = element_line(colour = 'black', linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(angle = 0, size = 10, vjust = 0.5, margin = margin(t = 3, r = 0, b = 0, l = 0)),
      axis.text.y = element_text(angle = 0, size = 10, vjust = 0.5, margin = margin(t = 0, r = 3, b = 0, l = 0)),
      axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 5, b = 0, l = 5)),
      axis.title.x = element_text(size = 12, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      plot.margin=unit(c(0,0,0,0),"cm")
    )

# project number of discoveries and explained variance
message(' - project number of discoveries and explained explained')

  # set lowerlimit, upperlimit, and number of data points that shall be created
  lowerlimit = 0
  upperlimit = 10000000
  datapoints = 1000

  # Get predicitions for lowerlimit to upperlimit
  GVpercentage = Numdiscoveries = matrix(NA,nrow=datapoints+1,ncol=length(traits)+length(reftraits)+1)
  k = 0
  for (i in seq(lowerlimit, upperlimit, length.out = datapoints+1)){
    k = k+1
    GVpercentage[k,1] = Numdiscoveries[k,1] = i
    for (j in 1:length(traits)) {
      tmp = get(sprintf('%s.fit',traits[j]))
      est = tmp$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates`
      v = tmp$estimates$`Covariance matrix of parameter estimates`
      project = projection(est,v,n=i,CI=FALSE);
      GVpercentage[k,j+1] = project$GVpercentage[1]
      Numdiscoveries[k,j+1] = project$Numdiscoveries[1]
    }
    for (j in 1:length(reftraits)) {
      tmp = get(sprintf('%s.fit',reftraits[j]))
      est = tmp$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates`
      v = tmp$estimates$`Covariance matrix of parameter estimates`
      project = projection(est,v,n=i,CI=FALSE);
      GVpercentage[k,j+length(traits)+1] = project$GVpercentage[1]
      Numdiscoveries[k,j+length(traits)+1] = project$Numdiscoveries[1]
    }
  }
  colnames(GVpercentage) = colnames(Numdiscoveries) = c('N',traits,reftraits)

  # plot expected number of discoveries as a function of sample size
  dfdiscoveries = reshape2::melt(data.frame(Numdiscoveries), id="N", variable.name="trait", value.name="ndiscov", na.rm = F)
  dfdiscoveries$trait = factor(dfdiscoveries$trait, levels = c(reftraits,traits), labels = c(reftraitLabels,traitLabels))
  pldiscoveries = ggplot() +
    geom_line(data=dfdiscoveries, aes(x = N/1000000, y = ndiscov, color = trait, linetype = trait), linewidth = 0.50) +
    scale_colour_manual(values = lineColours) + # c('grey10','grey60','grey90')
    scale_linetype_manual(values = lineTypes) + # c('grey10','grey60','grey90')
    labs(x = 'Sample size', y =  'Expected discoveries') +
    geom_segment(aes(x=0,xend=2,y=-Inf,yend=-Inf), colour = 'black', linewidth = 0.25) +
    geom_segment(aes(y=0,yend=5000,x=-Inf,xend=-Inf), colour = 'black', linewidth = 0.25) +
    scale_x_continuous(labels=c('0', '0.5M', '1M', '1.5M', '2M'), breaks=seq(0, 2, 0.5), limits = c(0,2)) +
    scale_y_continuous(labels=c('0','1k','2k','3k','4k','5k'), breaks=seq(0, 5000, 1000), limits = c(0,5000)) + # scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE), limits = c(0,5000)) +
    theme_light() +
    theme( 
      legend.position = 'none',
      legend.key = element_blank(),
      legend.key.size = unit(0.5, 'cm'),
      legend.background=element_blank(),
      legend.title = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks.length=unit(.2, "cm"),
      axis.ticks = element_line(colour = 'black', linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(angle = 0, size = 10, vjust = 0.5, margin = margin(t = 3, r = 0, b = 0, l = 0)),
      axis.text.y = element_text(angle = 0, size = 10, vjust = 0.5, margin = margin(t = 0, r = 3, b = 0, l = 0)),
      axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 5, b = 0, l = 5)),
      axis.title.x = element_text(size = 12, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      plot.margin=unit(c(0,0.75,0,0),"cm")
    )

  # plot genetic variance explained by discoveries as a function of sample size
  dfvariance = reshape2::melt(data.frame(GVpercentage), id="N", variable.name="trait", value.name="variance", na.rm = F)
  dfvariance$trait = factor(dfvariance$trait, levels = c(reftraits,traits), labels = c(reftraitLabels,traitLabels))
  plvariance = ggplot() +
    geom_line(data=dfvariance, aes(x = N/1000000, y = variance, color = trait, linetype = trait), linewidth = 0.50) +
    scale_colour_manual(values = lineColours) + # c('grey10','grey60','grey90')
    scale_linetype_manual(values = lineTypes) + # c('grey10','grey60','grey90')
    labs(x = "Sample size", y = "Genetic variance explained (%)") +
    geom_segment(aes(x=0,xend=2,y=-Inf,yend=-Inf), colour = 'black', linewidth = 0.25) +
    geom_segment(aes(y=0,yend=100,x=-Inf,xend=-Inf), colour = 'black', linewidth = 0.25) +
    scale_x_continuous(labels=c('0', '0.5M', '1M', '1.5M', '2M'), breaks=seq(0, 2, 0.5), limits = c(0,2)) +
    scale_y_continuous(breaks=seq(0, 100, 20))  +
    theme_light() +
    theme( 
      legend.position = 'none',
      legend.key = element_blank(),
      legend.key.size = unit(0.5, 'cm'),
      legend.background=element_blank(),
      legend.title = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks.length=unit(.2, "cm"),
      axis.ticks = element_line(colour = 'black', linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(angle = 0, size = 10, vjust = 0.5, margin = margin(t = 3, r = 0, b = 0, l = 0)),
      axis.text.y = element_text(angle = 0, size = 10, vjust = 0.5, margin = margin(t = 0, r = 3, b = 0, l = 0)),
      axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 5, b = 0, l = 5)),
      axis.title.x = element_text(size = 12, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      plot.margin=unit(c(0,0.75,0,0),"cm")
    )

# save file
message(sprintf(' - saving %s.png',outFile))
options(warn=-1)
pl = pldensity + pldiscoveries + plvariance + 
  plot_layout(widths = c(3,1,1)) + 
  plot_annotation(tag_levels ='a') & theme(plot.tag = element_text(size = 18, hjust = 0, vjust = 0, face = 'bold'))
png(width = 10, height = 3.5, units = "in", res = 300, filename = sprintf('%s.png',outFile))
pl
invisible(dev.off())
options(warn=0)

# save log file
message(sprintf(' - writing %s.log', outFile))
sink(sprintf('%s.log', outFile))
sink(stdout(), type = "message")
message(logInfo)
sink()
system(sprintf('chmod -R 770 %s*', outFile))
message('--- Completed: plot results of genetic effect size distribution analysis ---')

