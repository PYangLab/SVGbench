###################################
# call SVGs for benchmarking datasets using MERINGUE
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(glue)
library(parallel)
library(BiocParallel)
library(scuttle)
library(here)

# load functions to call SVGs
source(here("02_run methods", "functions", "call_variable_genes.R"))

path = "../../../SpatialData/SpatialBenchmark/DataUploadSubset/"
out_path = "../SVG/meringue/"

f = list.dirs(path, recursive=FALSE)

#### RUN MERINGUE
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
  
  # normalize
  normMat <- normalizeCounts(counts = countMat, log=FALSE, verbose=TRUE)
  
  # run method
  gc1 <- gc(reset = TRUE)
  time = system.time({ svg_meringue = callSVG.MERINGUE(normMat, spatial_coord)})
  gc2 <- gc()
  
  computation = data.frame(
    time = time[["elapsed"]],
    n_cell = ncol(normMat),
    n_gene = nrow(normMat),
    dataset = data_idx,
    Peak_RAM_Used_MiB = sum(gc2[,6]-gc1[,6]),
    method = "meringue"
  )
  write.table(computation, file=glue(out_path, data_idx, "_time_meringue.tsv"), sep="\t")
  write.table(svg_meringue, file=glue(out_path, data_idx, "_svg_meringue.tsv"), sep="\t")
  rm(normMat, countMat)
  
}, mc.cores = n_cores, mc.set.seed = TRUE)
