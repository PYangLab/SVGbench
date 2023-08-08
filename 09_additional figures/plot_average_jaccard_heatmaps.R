###################################
# Create jaccard index heatmaps quantifying SVG overlap
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(jaccard)
library(pheatmap)
library(paletteer)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(corrplot)
library(here)

# load jaccard pairwise matrices
load(here("data", "res_top_average jaccard.RData"))

pal = paletteer::palettes_d$ggthemes$Jewel_Bright[-3]
names(pal) = c("giotto_kmeans","giotto_rank", "seurat_moransi", "meringue", "nnSVG", "somde", "sparkx", "spatialde")

# create matrices for plotting

#### top 200
top200 = lapply(1:length(jaccard.200), function(i){
  mat = jaccard.200[[i]]
  dat.long = reshape2::melt(mat)
  dat.long$dataset = names(jaccard.200)[i]
  dat.long$comp = paste(dat.long$Var1, dat.long$Var2, sep=".")
  return(dat.long)
})
top200 = do.call(rbind, top200)

# calculate mean pairwise jaccard across all datasets
top200.mean = top200 %>%
  group_by(comp) %>%
  summarise(jaccard=mean(value, na.rm=T)) %>%
  mutate(comp1 = sapply(strsplit(comp, "\\."), "[[", 1),
         comp2 = sapply(strsplit(comp, "\\."), "[[", 2))

# pivot into wide matrix
top200.wide = top200.mean %>%
  pivot_wider(id_cols=comp1, names_from=comp2, values_from = jaccard)
top200.wide = as.matrix(top200.wide[,-1])
rownames(top200.wide) = colnames(top200.wide)

#### top 500
top500 = lapply(1:length(jaccard.500), function(i){
  mat = jaccard.500[[i]]
  dat.long = reshape2::melt(mat)
  dat.long$dataset = names(jaccard.500)[i]
  dat.long$comp = paste(dat.long$Var1, dat.long$Var2, sep=".")
  return(dat.long)
})
top500 = do.call(rbind, top500)

# calculate mean pairwise jaccard across all datasets
top500.mean = top500 %>%
  group_by(comp) %>%
  summarise(jaccard=mean(value, na.rm=T)) %>%
  mutate(comp1 = sapply(strsplit(comp, "\\."), "[[", 1),
         comp2 = sapply(strsplit(comp, "\\."), "[[", 2))

# pivot into wide matrix
top500.wide = top500.mean %>%
  pivot_wider(id_cols=comp1, names_from=comp2, values_from = jaccard)
top500.wide = as.matrix(top500.wide[,-1])
rownames(top500.wide) = colnames(top500.wide)

#### top 1000
top1000 = lapply(1:length(jaccard.1000), function(i){
  mat = jaccard.1000[[i]]
  dat.long = reshape2::melt(mat)
  dat.long$dataset = names(jaccard.1000)[i]
  dat.long$comp = paste(dat.long$Var1, dat.long$Var2, sep=".")
  return(dat.long)
})
top1000 = do.call(rbind, top1000)

# calculate mean pairwise jaccard across all datasets
top1000.mean = top1000 %>%
  group_by(comp) %>%
  summarise(jaccard=mean(value, na.rm=T)) %>%
  mutate(comp1 = sapply(strsplit(comp, "\\."), "[[", 1),
         comp2 = sapply(strsplit(comp, "\\."), "[[", 2))

# pivot into wide matrix
top1000.wide = top1000.mean %>%
  pivot_wider(id_cols=comp1, names_from=comp2, values_from = jaccard)
top1000.wide = as.matrix(top1000.wide[,-1])
rownames(top1000.wide) = colnames(top1000.wide)

#### all sig SVGs
topall = lapply(1:length(jaccard.all), function(i){
  mat = jaccard.all[[i]]
  dat.long = reshape2::melt(mat)
  dat.long$dataset = names(jaccard.all)[i]
  dat.long$comp = paste(dat.long$Var1, dat_long$Var2, sep=".")
  return(dat.long)
})
topall = do.call(rbind, topall)

# calculate mean pairwise jaccard across all datasets
topall.mean = topall %>%
  group_by(comp) %>%
  summarise(jaccard=mean(value, na.rm=T)) %>%
  mutate(comp1 = sapply(strsplit(comp, "\\."), "[[", 1),
         comp2 = sapply(strsplit(comp, "\\."), "[[", 2))

# pivot into wide matrix
topall.wide = topall.mean %>%
  pivot_wider(id_cols=comp1, names_from=comp2, values_from = jaccard)
topall.wide = as.matrix(topall.wide[,-1])
rownames(topall.wide) = colnames(topall.wide)

#### visualise heatmaps

pheatmap(t(top200.wide), col=colorRampPalette(brewer.pal(9, "RdPu"))(100), cluster_rows = T, cluster_cols = T, cellwidth = 20, cellheight = 20, filename=here("plots", "top200_jaccard.pdf"))

pheatmap(t(top500.wide), col=colorRampPalette(brewer.pal(9, "BuPu"))(100), cluster_rows = T, cluster_cols = T, cellwidth = 20, cellheight = 20, filename=here("plots", "top500_jaccard.pdf"))

pheatmap(t(top1000.wide), col=colorRampPalette(brewer.pal(9, "GnBu"))(100), cluster_rows = T, cluster_cols = T, cellwidth = 20, cellheight = 20, filename=here("plots", "top1000_jaccard.pdf"))

pheatmap(t(topall.wide), col=colorRampPalette(brewer.pal(9, "Greys"))(100), cluster_rows = T, cluster_cols = T, cellwidth = 20, cellheight = 20, filename=here("plots", "topall_jaccard.pdf"))
