###################################
# call SVGs for benchmarking datasets using Giotto
# Carissa Chen, updated Jul 2023
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
out_path = "../SVG/giotto/"

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

  # run method
  gc1 <- gc(reset = TRUE)
  time = system.time({svg_giotto_kmeans = callSVG.Giotto(data_giotto, svg_method = "kmeans", runHVG = TRUE, n_cores=32)})
  gc2 <- gc()

  computation = data.frame(
    time = time[['elapsed']],
    n_cell = ncol(data_giotto@raw_exprs),
    n_gene = nrow(data_giotto@raw_exprs),
    dataset = data_idx,
    Peak_RAM_Used_MiB = sum(gc2[,6]-gc1[,6]),
    method = "giotto_kmeans"
  )
  write.table(computation, file=glue(out_path, data_idx, '_time_giotto_kmeans.tsv'), sep='\t')

  # run method
  gc1 <- gc(reset = TRUE)
  time = system.time({svg_giotto_rank = callSVG.Giotto(data_giotto, svg_method = "rank", runHVG = TRUE, n_cores=32)})
  gc2 <- gc()

  computation = data.frame(
    time = time[['elapsed']],
    n_cell = ncol(data_giotto@raw_exprs),
    n_gene = nrow(data_giotto@raw_exprs),
    dataset = data_idx,
    Peak_RAM_Used_MiB = sum(gc2[,6]-gc1[,6]),
    method = "giotto_rank"
  )
  write.table(computation, file=glue(out_path, data_idx, '_time_giotto_rank.tsv'), sep='\t')

  write.table(svg_giotto_kmeans, file=glue(out_path, data_idx, "_svg_giotto_kmeans.tsv"), sep="\t")
  write.table(svg_giotto_rank, file=glue(out_path, data_idx, "_svg_giotto_rank.tsv"), sep="\t")
  rm(countMat, data_giotto)

}, mc.cores = n_cores, mc.set.seed = TRUE)