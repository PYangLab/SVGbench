###################################
# Script to preprocess slideseqV1 data
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(tidyverse)
library(glue)

f = list.files("../../../SpatialData/Slideseq/")
outdir = "../../../SpatialData/Slideseq/"

geneCellStat <- matrix(NA, nrow = length(f), ncol=4)
data.idx = sapply(strsplit(f, "\\.|/"), "[[", 1)

i=2
  
test = read.csv(file=paste0(f[i], "/MappedDGEForR.csv"), head=T)
meta = read.csv(file=paste0(f[i], "/BeadLocationsForR.csv"), head=T)

rownames(test) = test[,1]
test = test[,-1]
test = t(test)
rownames(test) = 1:nrow(test)

x = meta[,2]
y = meta[,3]
meta = data.frame(x,y)
rownames(meta) = 1:nrow(meta)

geneCellStat[i,1:2] <- dim(t(test))

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
