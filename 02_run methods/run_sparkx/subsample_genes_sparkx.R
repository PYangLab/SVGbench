###################################
# Script to subsample genes and call SVGs using SPARK-X
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(SingleCellExperiment)
library(glue)
library(parallel)
library(SpatialExperiment)
library(BiocParallel)
library(here)

source(here("02_run methods", "functions", "call_variable_genes.R"))

path = "../../../SpatialData/SpatialBenchmark/DataUploadSubset/"
out_path = "../SVG/subsample_genes/sparkx/"

f = list.dirs(path, recursive=FALSE)
f = f[grep("Dataset21|Dataset41", f)]

gp = 0.5

#### RUN SPARK
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

    # create objects
    data_spe = createSPEObject(countMat, spatial_coord, meta)

    # RUN SPARK
    svg_spark = callSVG.SPARKX(data_spe, n_cores = 32)

    write.table(svg_spark, file=glue(out_path, data_idx, "_svg_spark_subset_genes_", gp, ".tsv"), sep="\t")
    rm(countMat, data_spe)
  }, mc.cores = n_cores, mc.set.seed = TRUE)

