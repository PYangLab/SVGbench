###################################
# Script to preprocess benchmarking datasets
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(tidyverse)
library(glue)

path = "../../../SpatialData/SpatialBenchmark/DataUpload/"
outdir = "../../../SpatialData/SpatialBenchmark/DataUploadSubset/"

f <- list.files(path)
geneCellStat <- matrix(NA, nrow = length(f), ncol=4)

# run for each raw file
for(i in 1:length(f)) {
  
  print(i)
  test <- read.delim2(paste(path, f[i], "/Spatial_count.txt", sep=""), sep = "\t", header = T)
  meta <- read.delim2(paste(path, f[i], "/Locations.txt", sep=""), sep = "\t", header = T)
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
    
    dir.create(glue::glue(outdir, f[[i]]))
    write.table(test_filtered, file = glue::glue(outdir, f[[i]], "/Spatial_count.csv"), sep =",")
    write.table(meta_filtered, file = glue::glue(outdir, f[[i]], "/Locations.csv"), sep =",")
    
  }
  
}