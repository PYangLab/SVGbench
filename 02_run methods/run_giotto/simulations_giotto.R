###################################
# call SVGs for simulation datasets using Giotto
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(glue)
library(parallel)
library(SpatialExperiment)
library(BiocParallel)
library(here)

source(here("02_run methods", "functions", "call_variable_genes.R"))

path = "../../../SpatialData/Simulations/"
out_path = "../SVG/simulation/giotto/"

f = list.dirs(path, recursive=FALSE)

#### RUN GIOTTO
n_cores = 1

set.seed(1)
mclapply(1:length(f), function(i) {
  
  data_idx = gsub(path, "", f[[i]])
  data_idx = gsub("/", "", data_idx)
  print(glue(i, '........................', data_idx))
  
  countMat = t(read.table(glue(f[[i]], "/Spatial_count.csv"), sep = ",", header = TRUE))
  spatial_coord =  read.table(glue(f[[i]], "/Locations.csv"), sep = ",")
  rownames(spatial_coord) = colnames(countMat) = paste0("cell_", 1:ncol(countMat))
  meta = spatial_coord
  
  # filter
  countMat <- countMat[!duplicated(rownames(countMat)),]
  colzero = colSums(countMat) == 0
  rowzero = rowSums(countMat) == 0
  countMat = countMat[, !colzero]
  countMat = countMat[!rowzero, ]
  meta = meta[colnames(countMat),]
  spatial_coord = spatial_coord[colnames(countMat),]
  
  # create objects
  data_giotto = createGobject(countMat, spatial_coord, meta)
  
  # run giotto
  svg_giotto_kmeans = callSVG.Giotto(data_giotto, svg_method = "kmeans", runHVG = TRUE, n_cores=32)
  svg_giotto_rank = callSVG.Giotto(data_giotto, svg_method = "rank", runHVG = TRUE, n_cores=32)
  
  write.table(svg_giotto_kmeans, file=glue(out_path2, data_idx, "_svg_giotto_kmeans.tsv"), sep="\t")
  write.table(svg_giotto_rank, file=glue(out_path2, data_idx, "_svg_giotto_rank.tsv"), sep="\t")
  rm(countMat, data_giotto)
  
}, mc.cores = n_cores, mc.set.seed = TRUE)