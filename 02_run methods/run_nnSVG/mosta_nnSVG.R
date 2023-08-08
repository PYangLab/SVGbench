###################################
# call SVGs for E9.5 dataset using nnSVG
# Carissa Chen, updated Jul 2023
###################################

library(SingleCellExperiment)
library(glue)
library(parallel)
library(SpatialExperiment)
library(BRISC)
library(BiocParallel)
library(zellkonverter)
library(here)

source(here("02_run methods", "functions", "call_variable_genes.R"))
path = "../../../SpatialData/MOSTA/"

f=list.files(path, pattern="^E9.5")

out_path = "../SVG/MOSTA/"

#### RUN nnSVG
i=1
  
print(glue(i, '........................', f[[i]]))
data_idx = gsub(".MOSTA.h5ad", "", f[[i]])

# read in data
mosta <- zellkonverter::readH5AD(file=glue(path, f[[i]]), reader="python")
# filter any genes in fewer than 30 cells
mosta <- mosta[rowData(mosta)$n_cells > 30, ]
# filter any cells where more than 50% of counts are contributed by top 50 genes
mosta <- mosta[,colData(mosta)$pct_counts_in_top_50_genes < 50]

mosta <- mosta[, !duplicated(colnames(mosta))]
mosta <- mosta[!duplicated(rownames(mosta)),]

countMat <- assay(mosta, "count")
colnames(countMat) <- colnames(mosta)
rownames(countMat) <- rownames(mosta)

meta = as.data.frame(reducedDim(mosta, "spatial"))
meta$annotation = mosta$annotation
rownames(meta) <- colnames(countMat)
colnames(meta) <- c("x", "y", "annotation")
spatial_coord = meta[, 1:2]

# create objects
data_spe = createSPEObject(countMat, spatial_coord, meta)

svg_nnSVG = callSVG.nnSVG(data_spe, n_threads = 32)

write.table(svg_nnSVG, file=glue(out_path, data_idx, "_svg_nnSVG.tsv"), sep="\t")
