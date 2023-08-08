###################################
# call SVGs for E9.5 dataset using MERINGUE
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(SingleCellExperiment)
library(glue)
library(parallel)
library(BiocParallel)
library(zellkonverter)
library(scuttle)
library(here)

source(here("02_run methods", "functions", "call_variable_genes.R"))
path = "../../../SpatialData/MOSTA/"

f=list.files(path, pattern="^E9.5")

out_path = "../SVG/MOSTA/meringue/"

#### RUN MERINGUE

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

# normalize
normMat <- normalizeCounts(counts = countMat, log=FALSE, verbose=TRUE)

svg_meringue = callSVG.MERINGUE(normMat, spatial_coord)

write.table(svg_meringue, file=glue(out_path, data_idx, "_svg_meringue.tsv"), sep="\t")
