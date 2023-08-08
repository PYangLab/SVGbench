###################################
# Script to preprocess slideseqV2 data
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(tidyverse)
library(Seurat)
library(biomaRt)

f = list.files("../../../SpatialData/SlideseqV2")
outdir = "../../../SpatialData/SlideseqV2"

geneCellStat <- matrix(NA, nrow = length(f), ncol=4)
data.idx = sapply(strsplit(f, "\\.|/"), "[[", 1)

# get gene annotations
mart <- useMart("ensembl", "hsapiens_gene_ensembl")
ensembl.gene <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), mart=mart)
ensembl.human.gene <- ensembl.gene[ensembl.gene[,2] != "",]
rownames(ensembl.human.gene) <- make.names(ensembl.human.gene[,1], unique = TRUE)

for (i in 1:2) {
  
  seu = readRDS(f[i])
  test = GetAssayData(object = seu, slot = "counts")
  
  # convert human ENSEMBL ID to gene ID
  geneid = ensembl.human.gene[rownames(test), "hgnc_symbol"]
  rownames(test) = geneid
  test = t(as.matrix(test))
  rownames(test) = 1:nrow(test)
  
  x = Embeddings(object = seu, reduction = "spatial")[,1]
  y = Embeddings(object = seu, reduction = "spatial")[,2]
  meta = data.frame(x,y)
  rownames(meta) = 1:nrow(meta)
  
  geneCellStat[i,1:2] <- dim(t(test))
  
  # count genes expressed in 30 or more cells
  gene_idx = apply(test, 2, function(x){sum(x!=0) >= 30})
  geneCellStat[i,3] <- sum(gene_idx)
  # count cells that negate "where more than 50% of counts are contributed by top 50 genes"
  cell_idx = apply(test, 1, function(x){
    sum(sort(as.numeric(x), decreasing=TRUE)[1:50]) < (sum(as.numeric(x)) / 2)
  })
  geneCellStat[i,4] <- sum(cell_idx)
  
  if (geneCellStat[i,3] > 2000 & geneCellStat[i,4] > 200) {
    
    test_filtered = test[cell_idx, gene_idx]
    meta_filtered = meta[cell_idx,]
    
    dir.create(glue::glue(outdir, data.idx[i]))
    write.table(test_filtered, file = glue::glue(outdir, data.idx[i], "/Spatial_count.csv"), sep =",")
    write.table(meta_filtered, file = glue::glue(outdir, data.idx[i], "/Locations.csv"), sep =",")
    
  }
}