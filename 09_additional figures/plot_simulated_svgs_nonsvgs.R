#################
# Plot simulated spatial and non-spatial genes
# Last updated: Carissa Chen Jul 2023
#################

# load libraries
library(viridis)
library(SingleCellExperiment)
library(ggplot2)
library(scran)
library(parallel)
library(dplyr)
library(reshape2)
library(scales)
library(tidyr)
library(here)

path="../Simulation/test_scDesign3/"
f = list.files(path)

lapply(6, function(x) {
  
  print(x)
  load(paste0(path, f[[x]]), verbose = TRUE)
  print(dim(simu_sce))
  
  de_idx = rownames(simu_sce)[rowData(simu_sce)$is.de == TRUE]
  
  if (length(de_idx) < 60) {
    n = length(de_idx)
    n_col = 8
  } else {
    n = 64
    n_col = 8
  }
  de_genes = de_idx[sample(1:length(de_idx), n)]
  loc = colData(simu_sce)[,c("spatial1","spatial2")]
  expre = lapply(de_genes, function(x){
    curr = as.matrix(counts(simu_sce)[x,])
    curr = log1p(curr)
    return(rescale(curr))
  })
  long = do.call(rbind, expre)
  long = as.data.frame(long)
  colnames(long) <- "Expression"
  long$gene = do.call(c, lapply(de_genes, function(x){rep(x,dim(expre[[1]])[1])}))
  long$x = rep(loc[,1],n)
  long$y = rep(loc[,2],n)
  p = as_tibble(long, rownames = "Cell") %>% 
    ggplot(aes(x = x, y = y, color = Expression)) +
    ggrastr::geom_point_rast(size = 0.2) +
    facet_wrap(~gene, ncol = n_col) + 
    scale_colour_gradientn(colors = viridis_pal(option = "magma")(10), limits=c(0, 1)) + 
    coord_fixed(ratio = 1) + theme(axis.text.x = element_text(angle = 45))
  
  fn = here("plots", paste0("simulated_svgs_", gsub(".RData", "", f[[x]]), ".pdf"))
  ggsave(p, filename = fn, width = 15, height = 15)
  
  non_de_idx = rownames(simu_sce)[rowData(simu_sce)$is.de != TRUE]
  non_de_genes = non_de_idx[sample(1:length(non_de_idx), n)]
  loc = colData(simu_sce)[,c("spatial1","spatial2")]
  expre = lapply(non_de_genes, function(x){
    curr = as.matrix(counts(simu_sce)[x,])
    curr = log1p(curr)
    return(rescale(curr))
  })
  long = do.call(rbind, expre)
  long = as.data.frame(long)
  colnames(long) <- "Expression"
  long$gene = do.call(c, lapply(non_de_genes, function(x){rep(x,dim(expre[[1]])[1])}))
  long$x = rep(loc[,1],n)
  long$y = rep(loc[,2],n)
  p = as_tibble(long, rownames = "Cell") %>% 
    ggplot(aes(x = x, y = y, color = Expression)) +
    ggrastr::geom_point_rast(size = 0.4) +
    facet_wrap(~gene, ncol = n_col)+
    scale_colour_gradientn(colors = viridis_pal(option = "magma")(10), limits=c(0, 1)) + coord_fixed(ratio = 1) +
    theme(axis.text.x = element_text(angle = 45))
  
  fn = here("plots", paste0("simulated_nonsvgs_", gsub(".RData", "", f[[x]]), ".pdf"))
  ggsave(p, filename = fn, width = 15, height = 15)
})
