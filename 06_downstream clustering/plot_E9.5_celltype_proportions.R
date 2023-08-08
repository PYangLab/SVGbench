###################################
# Plot ground truth cell type proportions for E9.5
# Carissa Chen, updated Jul 2023
###################################

library(SeuratDisk)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)

seuratObj <- LoadH5Seurat("../../../SpatialData/MOSTA/E9.5.MOSTA.h5Seurat")
sce <- Seurat::as.SingleCellExperiment(seuratObj)

cl_full <- as.character(sce$annotation)

# create plotting dataframe
dfprop = data.frame(
  cty = cl_full,
  group = 1
)

set.seed(1)
my_col = sample(RColorBrewer::brewer.pal(12, "Paired"), 12, replace = FALSE)
my_col = sample(ggthemes::ptol_pal()(12), 12, replace = FALSE)
my_col = ggthemes::ptol_pal()(12)

p = dfprop %>%
  ggplot(aes(x = group, fill = cty)) + 
  geom_bar(position = "fill") + 
  scale_fill_manual(values = my_col) + 
  theme_classic()

fn = here("plots", "E9.5_ct_prop_plot.pdf")
ggsave(p, filename = fn, width = 5, height = 5)
