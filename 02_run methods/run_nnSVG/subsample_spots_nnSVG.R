###################################
# Script to subsample spots and call SVGs using nnSVG
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(glue)
library(parallel)
library(SpatialExperiment)
library(BRISC)
library(BiocParallel)
library(here)

source(here("02_run methods", "functions", "call_variable_genes.R"))

path = "../../../SpatialData/SpatialBenchmark/DataUploadSubset/"
out_path = "../SVG/subsample_spots/nnSVG/"

f = list.dirs(path, recursive=FALSE)

cp = 0.8

#### RUN nnSVG
n_cores = 1

RNGkind("L'Ecuyer-CMRG")
set.seed(2022)
mclapply(1:length(f), function(i) {

  #set.seed(r)
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

  # create objects
  data_spe = createSPEObject(countMat, spatial_coord, meta)

  svg_nnSVG = callSVG.nnSVG(data_spe, n_threads = 32, k=k)

  write.table(svg_nnSVG, file=glue(out_path, data_idx, "_svg_nnSVG_subsample_spots_", cp, ".tsv"), sep="\t")
  rm(countMat, data_spe)

}, mc.cores = n_cores, mc.set.seed = TRUE)
