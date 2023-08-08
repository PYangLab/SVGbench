###################################
# Script to subsample spots and call SVGs
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(glue)
library(parallel)
library(SpatialExperiment)
library(BiocParallel)

# load functions to call SVGs
source(here("02_run methods", "functions", "call_variable_genes.R"))

path = "../../../SpatialData/SpatialBenchmark/DataUpload/subsample/"
out_path = "../SVG/subsample_spots/giotto/"

f = list.dirs(path, recursive=FALSE)
f = f[grep("Dataset21|Dataset41", f)]

cp = 0.8

#### RUN GIOTTO
n_cores = 1

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

  # run giotto
  data_giotto = createGobject(countMat, spatial_coord, meta)

  svg_giotto_kmeans = callSVG.Giotto(data_giotto, svg_method = "kmeans", runHVG = TRUE, n_cores=32)
  svg_giotto_rank = callSVG.Giotto(data_giotto, svg_method = "rank", runHVG = TRUE, n_cores=32)

  write.table(svg_giotto_kmeans, file=glue(out_path1, data_idx, "_svg_giotto_kmeans_subsample_spots_", cp, ".tsv"), sep="\t")
  write.table(svg_giotto_rank, file=glue(out_path1, data_idx, "_svg_giotto_rank_subsample_spots_",cp, ".tsv"), sep="\t")
  rm(countMat, data_giotto)
}, mc.cores = n_cores, mc.set.seed = TRUE)
