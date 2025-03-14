#!/usr/bin/env Rscript

# ======================================
# === Create feature importance plot ===
# ======================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop(paste0('expected 4 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
modelFile=args[1] # modelFile="results/mri/ml.xgb/singlemodel/xgb_gm_models.RData"
outFile=args[2] # outFile="results/mri/ml.xgb/singlemodel/xgb_gm_importance.png"
width=as.numeric(args[3]) # width=7.7
height=as.numeric(args[4]) # height=3.0

message(paste0('\n--- Create feature importance plot | Settings ---',
               '\nmodelFile: ', modelFile,
               '\noutFile: ', outFile,
               '\nwidth: ', width,
               '\nheight: ', height,'\n'))

# attach packages to current R session
for (pkg in c('xgboost','data.table', 'dplyr', 'ggplot2','scales')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# load model
load(modelFile)
df = xgb.importance(model = model_tree)
df$Feature = gsub('V','PC',df$Feature)
df$Feature_num = as.numeric(gsub("PC", "", df$Feature))
df = df[order(df$Feature_num),]

# Create equally sized groups and aggregate gain
group_size = 10
df = df %>%
  mutate(Feature_group = paste0("PC", ((Feature_num - 1) %/% group_size) * group_size + 1, 
                           "-", ((Feature_num - 1) %/% group_size + 1) * group_size))
df.agg <- df %>%
  group_by(Feature_group) %>%
  summarise(Total_Gain = sum(Gain, na.rm = TRUE)) %>%
  mutate(Group_num = as.numeric(gsub("PC", "", gsub("-.*", "", Feature_group)))) %>%
  arrange(Group_num) %>%
  mutate(Feature_group = factor(Feature_group, levels = Feature_group)) 

# Custom labeling logic (clean numeric labels)
df.agg = df.agg %>%
  mutate(Group_label = ifelse(
    Group_num == 1 |
      (Group_num - 1) %% 50 == 40 |
      Group_num == max(Group_num),
    gsub("PC", "", Feature_group), ""))
labelled_data = df.agg %>%
  filter(Group_label != "")

# Visualize the aggregated gains
pl = ggplot(df.agg, aes(x = Feature_group, y = Total_Gain)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Principal Components",
       y = "Feature importance (gain)") +
  scale_x_discrete(
    labels = labelled_data$Group_label,
    breaks = labelled_data$Feature_group) + 
  scale_y_continuous(labels = function(x) sprintf("%.0f%%", x * 100), limits = c(0,0.52), expand = c(0, 0)) +
  geom_segment(aes(x = 0, xend = 50, y = 0, yend = 0), color = "black", size = 0.25) +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 0.5),color = "black", size = 0.25) +
  theme_minimal() +
  theme(plot.margin = margin(10, 20, 10, 10),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.ticks = element_line(color = "black", size = 0.25),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(angle = 0, hjust = 0.5))
# save plot
message(sprintf(' - saving  %s', outFile))
png(filename = outFile, width=width, height=height, units = "in", res = 600)
pl
invisible(dev.off())
system(sprintf('chmod 770 %s', outFile))
message('--- Completed: Create feature importance plot ---')
