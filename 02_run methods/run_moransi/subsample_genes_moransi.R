###################################
# Script to subsample genes and call SVGs using SEURAT Moran's I
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(SingleCellExperiment)
library(glue)
library(parallel)
library(SpatialExperiment)
library(BiocParallel)
library(Seurat)
library(here)

source(here("02_run methods", "functions", "call_variable_genes.R"))

path = "../../../SpatialData/SpatialBenchmark/DataUploadSubset/"
out_path = "../SVG/subsample_genes/moransi/"

f = list.dirs(path, recursive=FALSE)

gp = 0.5

#### RUN SEURAT Moran's I
n_cores = length(f)

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

    # subset by percentage
    sub_rows = sample(1:nrow(countMat), size=nrow(countMat)*gp, replace=FALSE)
    countMat = countMat[sub_rows, ]

    # create objects
    data_seu = createSeuObject(countMat)

    seu.moransi = callSVG.Seurat(data_seu[["SCT"]]@scale.data, method="moransi", spatial.loc=spatial_coord)

    write.table(seu.moransi, file=glue(out_path, data_idx, "_svg_seurat_moransi_subset_genes_", gp,".tsv"), sep="\t")

    rm(countMat, data_seu)
  }, mc.cores = n_cores, mc.set.seed = TRUE)
