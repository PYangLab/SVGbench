###################################
# Script to generate proportional plot of SVG overlaps
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(ComplexHeatmap)
library(dplyr)
library(here)

sigSVG = readRDS(here("data", "input_sigSVGs_per_method.rds"))
csv = read.csv(here("data", "dataset_metadata.csv"), head=T)

# binarize the overlapping SVGs
binary.mat <- list()
for (i in 1:length(sigSVG)) {
  binary.mat[[i]] = ComplexHeatmap::list_to_matrix(sigSVG[[i]])
}

# make a combination matrix
comb.mat <- list()
for (i in 1:length(binary.mat)) {
  comb.mat[[i]] = ComplexHeatmap::make_comb_mat(binary.mat[[i]], mode ="distinct")
}

# create plotting dataframe
df = list()
for (i in 1:length(comb.mat)) {
  n = names(comb_degree(comb.mat[[i]]))
  df[[i]] <- data.frame(intersections = n,
                        size=comb_size(comb.mat[[i]])[n],
                        degree=comb_degree(comb.mat[[i]])[n],
                        dataset=names(sigSVG)[i])
}

# plot proportional bar plot
df.prop = do.call(rbind, df)
df.prop$dataset = factor(df.prop$dataset, levels=csv$dataset)

ggplot(df.prop, aes(x=dataset, y=size, fill=as.character(degree))) +
  geom_bar(position="fill", stat="identity", width=0.8) +
  scale_fill_brewer(palette = "Greys") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

ggsave(filename=here("plots", "overlaps_proportions_plot.pdf"), width=10, height=8)
