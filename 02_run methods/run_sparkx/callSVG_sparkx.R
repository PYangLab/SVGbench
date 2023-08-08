###################################
# call SVGs for benchmarking datasets using SPARK-X
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(SingleCellExperiment)
library(glue)
library(parallel)
library(SpatialExperiment)
library(BiocParallel)

source(here("02_run methods", "functions", "call_variable_genes.R"))
path = "../../../SpatialData/SpatialBenchmark/"

out_path = "../SVG/sparkx/"
f = list.dirs(path, recursive=FALSE)

#### RUN SPARK
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
  data_spe = createSPEObject(countMat, spatial_coord, meta)

  # RUN SPARK-X
  gc1 <- gc(reset = TRUE)
  time = system.time({svg_spark = callSVG.SPARKX(data_spe, n_cores = 32)})
  gc2 <- gc()

  computation = data.frame(
    time = time[["elapsed"]],
    n_cell = ncol(data_spe),
    n_gene = nrow(data_spe),
    dataset = data_idx,
    Peak_RAM_Used_MiB = sum(gc2[,6]-gc1[,6]),
    method = "spark"
  )
  write.table(computation, file=glue(out_path, data_idx, "_time_spark.tsv"), sep="\t")

  write.table(svg_spark, file=glue(out_path, data_idx, "_svg_spark.tsv"), sep="\t")
  rm(countMat, data_spe)

}, mc.cores = n_cores, mc.set.seed = TRUE)