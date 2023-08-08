createGobject <- function(countMat, 
                          spatial_locs, 
                          cell_metadata, 
                          filter = FALSE, 
                          normalize = TRUE) {
  
  require(Giotto)
  
  python_path = NULL
  if(is.null(python_path)) {
    installGiottoEnvironment()
  }
  
  giotto = createGiottoObject(raw_exprs = countMat,
                                    spatial_locs = spatial_locs,
                                    cell_metadata = cell_metadata)
  
  if (filter == TRUE) {
    message("performing filtering...............")
    giotto <- filterGiotto(gobject = giotto,
                                 expression_threshold = 1,
                                 gene_det_in_min_cells = 50,
                                 min_det_genes_per_cell = 1000,
                                 expression_values = c('raw'),
                                 verbose = T)
    
  }
  if (normalize == TRUE) {
    message("performing normalization...............")
    giotto <- normalizeGiotto(gobject = giotto, scalefactor = 6000, verbose = FALSE)
    giotto <- addStatistics(gobject = giotto)
  }
  
  return(giotto)
}

createSPEObject <- function(countMat, 
                            spatial_locs, 
                            cell_metadata,
                            normalize = TRUE) {
  
  require(SpatialExperiment)
  require(S4Vectors)
  require(scran)
  
  rd = S4Vectors::DataFrame(symbol = rownames(countMat))
  
  data_spe = SpatialExperiment(
    assays = list(counts = countMat),
    colData = S4Vectors::DataFrame(cell_metadata),
    rowData = rd,
    spatialCoords = as.matrix(spatial_locs)
  )
  
  if (normalize == TRUE) {
    message("performing normalization...............")
    # set.seed(123)
    # qclus <- quickCluster(data_spe)
    # data <- computeSumFactors(data_spe, cluster = qclus)
    data_spe <- scater::logNormCounts(data_spe)
  }
  
  return(data_spe)
}

createSeuObject <- function(countMat, normalize=TRUE) {

  require(Seurat)

  seu = CreateSeuratObject(counts=countMat, assay="Spatial")

  if (normalize==TRUE) {
    message("performing normalization...............")
    seu = SCTransform(seu, assay = "Spatial", verbose = FALSE, return.only.var.genes=FALSE)
  }

  return(seu)
}


callSVG.Seurat <- function(logMat, method = "moransi", spatial.loc=spatial.loc) {
  
  require(Seurat)

  if(method == "moransi") {
    svf.info <- FindSpatiallyVariableFeatures(logMat, selection.method = "moransi",
                                          spatial.location = spatial.loc,
                                          nfeatures = nrow(seu))
    
    svf.info <- svf.info[order(svf.info[, 2], -abs(svf.info[, 1])), , drop = FALSE]
    return(svf.info)
  }
}

callSVG.MERINGUE = function(normMat, spatial_locs) {
  
  require(MERINGUE)
  
  w <- getSpatialNeighbors(spatial_locs, filterDist = NA)
  I <- getSpatialPatterns(normMat, w)
  
  return(I)
}


callSVG.Giotto <- function(gobject,
                           svg_method = c("kmeans", "rank"),
                           runHVG = FALSE,
                           runPCA = TRUE, n_cores = 1L) {
  
  if (runHVG == TRUE) {
    gobject <- Giotto::calculateHVG(gobject = gobject)
  }
  
  if (runPCA == TRUE) {
    gobject <- runPCA(gobject =  gobject, center = TRUE, scale_unit = TRUE)
  }
  
  gobject = createSpatialNetwork(gobject = gobject, minimum_k = 0)
  
  if (svg_method == "kmeans") {
    kmtest = Giotto::binSpect(gobject, do_parallel=TRUE, cores=n_cores)
    return(kmtest)
  } else {
    ranktest = Giotto::binSpect(gobject, bin_method = "rank", do_parallel=TRUE, cores=n_cores)
    return(ranktest)
  }
  
}

callSVG.SPARKX <- function(spe,
                          n_cores = 1L) {
  
  require(SPARK)
  require(SpatialExperiment)
  
  sp_count = counts(spe)
  spatial_loc = spatialCoords(spe)
  
  res <- sparkx(sp_count,spatial_loc, 
                numCores = n_cores,
                option="mixture")
  
  return(res$res_mtest)
  
}

callSVG.nnSVG <- function(spe,
                          seed = 1,
                          return_spe = FALSE,
                          n_threads = 1L, k=10) {
  
  set.seed(seed)
  spe <- nnSVG(spe, n_threads = n_threads, n_neighbors=k)
  
  if (return_spe == FALSE) {
    return(rowData(spe))
  } else {
    return(spe)
  }
  
}
