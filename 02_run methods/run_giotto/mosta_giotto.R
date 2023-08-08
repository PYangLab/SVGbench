###################################
# call SVGs for E9.5 dataset using Giotto
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(SingleCellExperiment)
library(glue)
library(parallel)
library(SpatialExperiment)
library(BiocParallel)
library(zellkonverter)
library(here)

# load functions to call SVGs
source(here("02_run methods", "functions", "call_variable_genes.R"))

path = "../../../SpatialData/MOSTA/"
out_path = "../../../SpatialData/MOSTA/"

f=list.files(path, pattern="^E9.5")

#### RUN GIOTTO

i=1

print(glue(i, '........................', f[[i]]))
data_idx = gsub(".MOSTA.h5ad", "", f[[i]])

# read in data
mosta <- zellkonverter::readH5AD(file=glue(path, f[[i]]), reader="R")
colnames(mosta) = mosta$cell_name
rownames(mosta) = rowData(mosta)$gene_short_name

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

# run giotto
data_giotto = createGobject(countMat, spatial_coord, meta)

svg_giotto_kmeans = callSVG.Giotto(data_giotto, svg_method = "kmeans", runHVG = TRUE, n_cores=32)

svg_giotto_rank = callSVG.Giotto(data_giotto, svg_method = "rank", runHVG = TRUE, n_cores=32)

write.table(svg_giotto_kmeans, file=glue(out_path, data_idx, "_svg_giotto_kmeans.tsv"), sep="\t")
write.table(svg_giotto_rank, file=glue(out_path, data_idx, "_svg_giotto_rank.tsv"), sep="\t")
