###################################
# Script to generate individual spearmans correlation heatmaps
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(here)
library(ggplot2)
library(corrplot)
library(glue)

data.list = readRDS(here("data", "input_stats_by_dataset.rds"))

cor.mat = list()
for (i in 1:length(data.list)) {
  SVG.combined <- do.call(cbind, data.list[[i]])
  
  SVG.sel <- SVG.combined[complete.cases(SVG.combined),]
  colnames(SVG.sel) <- methods
  cor.mat[[i]] = cor(SVG.sel, use="complete.obs", method="spearman")
  
  fn = here("plots", "spearmans_heatmap", paste0(names(data.list)[i], ".pdf"))
  
  pdf(file=fn, width=6, height=6)
  corrplot(cor.mat[[i]], method="circle", addCoef.col="gray70", number.cex = 1, title=names(data.list)[i])
  dev.off()
}