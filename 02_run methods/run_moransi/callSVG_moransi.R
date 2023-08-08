###################################
# call SVGs for benchmarking datasets using Moran's I
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
out_path1 = "../SVG/moransi/"

f = list.dirs(path, recursive=FALSE)

#### RUN SEURAT Moran's I
n_cores = length(f)

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
  data_seu = createSeuObject(countMat)
  
  # run method
  gc1 <- gc(reset = TRUE)
  time = system.time({seu.moransi = callSVG.Seurat(data_seu[["SCT"]]@scale.data, method="moransi", spatial.loc=spatial_coord)})
  gc2 <- gc()

  computation = data.frame(
    time = time[["elapsed"]],
    n_cell = ncol(data_seu),
    n_gene = nrow(data_seu),
    dataset = data_idx,
    Peak_RAM_Used_MiB = sum(gc2[,6]-gc1[,6]),
    method = "seurat_moransi"
  )
  write.table(computation, file=glue(out_path1, data_idx, "_time_seurat_moransi.tsv"), sep="\t")
  write.table(seu.moransi, file=glue(out_path2, data_idx, "_svg_seurat_moransi.tsv"), sep="\t")

  rm(countMat, data_seu)

}, mc.cores = n_cores, mc.set.seed = TRUE)
