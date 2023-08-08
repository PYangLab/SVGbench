###################################
# Script to create simulation datasets
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(viridis)
library(scDesign3)
library(SingleCellExperiment)
library(ggplot2)
library(Seurat)
library(SeuratObject)
library(scran)
library(parallel)
library(DESeq2)
library(BiocParallel)
library(PseudotimeDE)
library(dplyr)
library(tradeSeq)
library(reshape2)
library(spatialDE)
library(SPARK)
library(scales)
library(tidyr)

path = "../../../SpatialData/SpatialBenchmark/DataUploadSubset/"
f = list.files(path)
f = f[!grepl("zip", f)]

# check dimensions of datasets
dimList = mclapply(1:length(f), function(x) {
  
  print(x)
  count = read.delim(paste0(path, f[[x]], "/", list.files(paste0(path, f[[x]], "/"), pattern = "Spatial_count")))
  return(list(cells = nrow(count), genes = ncol(count)))
  rm(count)
  gc()
  
}, mc.cores = 5)
names(dimList) = f

sapply(dimList, function(x) x$genes)
sapply(dimList, function(x) x$cells)

# subset for datasets that have > 2000 & 10000 cells
dataset = names(which(sapply(dimList, function(x) x$cells) > 2000 & sapply(dimList, function(x) x$cells) < 10000))

# simulate for each individual dataset
lapply(1:length(dataset), function(x) {
  
  print(x)
  
  count = read.delim(paste0(path, dataset[[x]], "/", 
                            list.files(paste0(path, dataset[[x]], "/"), pattern = "Spatial_count")))
  loc = read.delim(paste0(path, dataset[[x]], "/", list.files(paste0(path, dataset[[x]], "/"),
                                                              pattern = "Locations")))
  colnames(loc) = c("spatial1", "spatial2")
  
  example_sce = SingleCellExperiment::SingleCellExperiment(assay = list(counts = t(count)),
                                                           colData = loc)
  example_sce = scater::logNormCounts(example_sce)
  
  example_seu = Seurat::as.Seurat(example_sce)
  example_seu <- NormalizeData(example_seu,
                               normalization.method = "LogNormalize",
                               scale.factor = 10000)
  example_seu <- FindVariableFeatures(example_seu,
                                      selection.method = "vst",
                                      nfeatures = 2000)
  all.genes <- rownames(example_seu)
  example_seu <- ScaleData(example_seu, features = all.genes)
  example_seu <- RunPCA(example_seu, features = VariableFeatures(object = example_seu))
  example_seu <- FindNeighbors(example_seu, dims = 1:10)
  example_seu <- FindClusters(example_seu, resolution = 0.5)
  example_sce$cell_type = example_seu$originalexp_snn_res.0.5
  colnames(example_sce) = paste0("cell", 1:ncol(example_sce))
  
  set.seed(1)
  n_features = 2000
  example_seu <- FindVariableFeatures(example_seu,
                                      selection.method = "vst",
                                      nfeatures = n_features)
  example_sce = example_sce[VariableFeatures(object = example_seu),]
  
  mt_idx<- grep("mt-|mt.",rownames(example_sce))
  if(length(mt_idx)!=0){
    example_sce   <- example_sce[-mt_idx,]
  }
  dim(example_sce)
  
  set.seed(1)
  example_data <- construct_data(
    sce = example_sce,
    assay_use = "counts",
    celltype = "cell_type",
    pseudotime = NULL,
    spatial = c("spatial1", "spatial2"),
    other_covariates = NULL,
    corr_by = "1"
  )
  
  example_marginal <- fit_marginal(
    data = example_data,
    predictor = "gene",
    mu_formula = "s(spatial1, spatial2, bs = 'gp', k= 50)",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 14,
    usebam = FALSE
  )
  
  set.seed(1)
  example_copula <- fit_copula(
    sce = example_sce,
    assay_use = "counts",
    marginal_list = example_marginal,
    family_use = "nb",
    copula = "gaussian",
    n_cores = 14,
    new_covariate = NULL,
    input_data = example_data$dat
  )
  
  example_para <- extract_para(
    sce = example_sce,
    marginal_list = example_marginal,
    n_cores = 14,
    family_use = "nb",
    new_covariate = NULL,
    data = example_data$dat
  )
  
  dev_explain <- sapply(example_marginal, function(x){
    sum = summary(x$fit)
    return(sum$dev.expl)
  })
  
  dev_ordered <- order(dev_explain, decreasing = TRUE)
  
  num_de <- round(length(dev_ordered)/10)
  ordered <- dev_explain[dev_ordered]
  de_idx <- names(ordered)[1:num_de]
  non_de_idx <- names(ordered)[-(1:num_de)]
  non_de_mat <- apply(example_para$mean_mat[,non_de_idx], 2, function(x){
    avg <- (max(x)+min(x))/2
    new_mean <- rep(avg, length(x))
    return(new_mean)
  })
  example_para$mean_mat[,non_de_idx] <- non_de_mat
  
  set.seed(1)
  example_newcount <- simu_new(
    sce = example_sce,
    mean_mat = example_para$mean_mat,
    sigma_mat = example_para$sigma_mat,
    zero_mat = example_para$zero_mat,
    quantile_mat = NULL,
    copula_list = example_copula$copula_list,
    n_cores = 14,
    family_use = "nb",
    input_data = example_data$dat,
    new_covariate = example_data$newCovariate,
    important_feature = rep(TRUE, dim(example_sce)[1])
  )
  logcounts(example_sce) <- log1p(counts(example_sce))
  simu_sce <- example_sce 
  counts(simu_sce) <- example_newcount
  logcounts(simu_sce) <- log1p(counts(simu_sce))
  
  rowData(simu_sce) = data.frame(
    is.de = rownames(simu_sce) %in% de_idx
  )
  
  save(simu_sce, file = paste0("../Simulation/test_scDesign3/", dataset[[x]], ".RData"))
})


# save spatial locations and counts as csv files
path="../Simulation/test_scDesign3/"
path_to_save="../Simulation/test_scDesign3/simulation_scDesign3/"
f = list.files(path)
  
lapply(1:length(f), function(x) {
  
  print(x)
  load(paste0(path, f[[x]]), verbose = TRUE)
  print(dim(simu_sce))
  countMatrix = t(counts(simu_sce))
  loc = colData(simu_sce)[, c("spatial1", "spatial2")]
  
  if (isTRUE(dir.create(paste0(path_to_save, gsub(".RData", "", f[[x]]))))) {
    dir.create(paste0(path_to_save, gsub(".RData", "", f[[x]])))
    write.table(loc, paste0(path_to_save, gsub(".RData", "", f[[x]]), "/Locations",".csv"), sep = ",")
    write.table(countMatrix, paste0(path_to_save, gsub(".RData", "", f[[x]]), "/Spatial_count",".csv"), sep = ",")
  }
  
})
