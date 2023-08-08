###################################
# Script to subsample genes and call SVGs using MERINGUE
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(SingleCellExperiment)
library(glue)
library(parallel)
library(SpatialExperiment)
library(BiocParallel)
library(scuttle)
library(here)

source(here("02_run methods", "functions", "call_variable_genes.R"))

path = "../../../SpatialData/SpatialBenchmark/DataUploadSubset/"
out_path1 = "../SVG/subsample_genes/meringue"

f = list.dirs(path, recursive=FALSE)
f = f[grep("Dataset21|Dataset41", f)]

gp = 0.5

#### RUN MERINGUE
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
    spatial_coord = spatial_coord[colnames(countMat),]
    
    # subset by percentage
    sub_rows = sample(1:nrow(countMat), size=nrow(countMat)*gp, replace=FALSE)
    countMat = countMat[sub_rows, ]
    
    # normalize
    normMat <- normalizeCounts(counts = countMat, log=FALSE, verbose=TRUE)
    
    svg_meringue = callSVG.MERINGUE(normMat, spatial_coord)
    
    write.table(svg_meringue, file=glue(out_path1, data_idx, "_svg_meringue_subset_genes_", gp,".tsv"), sep="\t")
    rm(countMat, normMat)
    
  }, mc.cores = n_cores, mc.set.seed = TRUE)
