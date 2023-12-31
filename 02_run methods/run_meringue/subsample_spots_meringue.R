###################################
# Script to subsample spots and call SVGs using MERINGUE
# Carissa Chen, updated Oct 2023
###################################

# load libraries
library(glue)
library(parallel)
library(BiocParallel)
library(scuttle)
library(here)

source(here("02_run methods", "functions", "call_variable_genes.R"))

path = "../../../SpatialData/SpatialBenchmark/DataUpload/subsample/"
out_path = "../SVG/subsample_spots/meringue/"

f = list.dirs(path, recursive=FALSE)

cp = 0.8

#### RUN MERINGUE
n_cores = 5

RNGkind("L'Ecuyer-CMRG")
set.seed(2022)
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
  
  # subsample by percentage
  sub_cols = sample(1:ncol(countMat), size=ncol(countMat)*cp, replace=FALSE)
  countMat = countMat[, sub_cols]
  spatial_coord = spatial_coord[sub_cols, ]
  meta = meta[sub_cols, ]
  
  # normalize
  normMat <- normalizeCounts(counts = countMat, log=FALSE, verbose=TRUE)
  
  svg_meringue = callSVG.MERINGUE(normMat, spatial_coord)
  write.table(svg_meringue, file=glue(out_path, data_idx, "_svg_meringue_subsample_spots_", cp, ".tsv"), sep="\t")
  rm(countMat, normMat)
  
}, mc.cores = n_cores, mc.set.seed = TRUE)
