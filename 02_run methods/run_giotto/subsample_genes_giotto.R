###################################
# Script to subsample genes and call SVGs using Giotto
# Carissa Chen, updated Oct 2023
###################################

# load libraries
library(SingleCellExperiment)
library(glue)
library(parallel)
library(SpatialExperiment)
library(BiocParallel)
library(here)

# load functions to call SVGs
source(here("02_run methods", "functions", "call_variable_genes.R"))

path = "../../../SpatialData/SpatialBenchmark/DataUploadSubset/"
out_path1 = "../SVG/subsample_genes/giotto/"

set.seed(1)
f = list.dirs(path, recursive=FALSE)

gp = 0.5

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

    # subset by percentage
    sub_rows = sample(1:nrow(countMat), size=nrow(countMat)*gp, replace=FALSE)
    countMat = countMat[sub_rows, ]

    # run giotto
    data_giotto = createGobject(countMat, spatial_coord, meta)

    svg_giotto_kmeans = callSVG.Giotto(data_giotto, svg_method = "kmeans", runHVG = TRUE, n_cores=32)

    svg_giotto_rank = callSVG.Giotto(data_giotto, svg_method = "rank", runHVG = TRUE, n_cores=32)

    write.table(svg_giotto_kmeans, file=glue(out_path1, data_idx, "_svg_giotto_kmeans_subset_genes_", gp, ".tsv"), sep="\t")
    write.table(svg_giotto_rank, file=glue(out_path1, data_idx, "_svg_giotto_rank_subset_genes_", gp,".tsv"), sep="\t")
    rm(countMat, data_giotto)

  }, mc.cores = n_cores, mc.set.seed = TRUE)