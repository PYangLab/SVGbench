###################################
# Script to preprocess DBiT-seq data
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(tidyverse)
library(glue)

f = list.files("../../../SpatialData/DBiT-seq/", pattern="*.tsv", recursive = TRUE)
outdir = "../../../SpatialData/DBiT-seq/"

geneCellStat <- matrix(NA, nrow = length(f), ncol=4)
data.idx = sapply(strsplit(f, "\\.|/"), "[[", 1)

for(i in 1:length(f)) {
  
  test = read.table(f[i], sep="\t", head=T)
  geneCellStat[i,1:2] <- dim(t(test))
  
  x = sapply(strsplit(test[,1], "x"), "[[", 1)
  y = sapply(strsplit(test[,1], "x"), "[[", 2)
  meta = data.frame(x,y)
  rownames(test) = 1:nrow(test)
  test = test[,-1]
  
  # count genes expressed in 30 or more cells
  gene_idx = apply(test, 2, function(x){sum(x!=0) >= 30})
  geneCellStat[i,3] <- sum(gene_idx)
  # count cells that negate "where more than 50% of counts are contributed by top 50 genes"
  cell_idx = apply(test, 1, function(x){
    sum(sort(as.numeric(x), decreasing=TRUE)[1:50]) < (sum(as.numeric(x)) / 2)
  })
  geneCellStat[i,4] <- sum(cell_idx)
  
  if (geneCellStat[i,3] > 2000 & geneCellStat[i,4] > 200) {
    
    test_filtered = test[cell_idx, gene_idx]
    meta_filtered = meta[cell_idx,]
    
    dir.create(glue::glue(outdir, data.idx[i]))
    write.table(test_filtered, file = glue::glue(outdir, data.idx[i], "/Spatial_count.csv"), sep =",")
    write.table(meta_filtered, file = glue::glue(outdir, data.idx[i], "/Locations.csv"), sep =",")
    
  }
  
}