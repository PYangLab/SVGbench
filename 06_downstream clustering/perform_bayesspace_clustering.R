###################################
# Run bayesspace clustering
# Carissa Chen, updated Oct 2023
###################################

library(BayesSpace)
library(parallel)
library(Seurat)
library(Giotto)
library(SeuratDisk)
library(SingleCellExperiment)
library(scran)
library(ggplot2)
library(dplyr)
library(smfishHmrf)

load(file=here("data", "clustering_input.RData"))

# perform BayesSpace spatial clustering
topList = seq(100,2000,100)

for (i in 1:10) {
  
  print(i)
  lapply(1:length(statsList), function(method_idx) {
    
    print(method_idx)
    np = 20 
    
    statsvalue = statsList[[method_idx]]
    pvalue = pvalueList[[method_idx]]
    method = names(statsList)[[method_idx]]
    
    df_top = do.call(rbind, mclapply(topList, function(top) {
      
      print(top)
      set.seed(i)
      s <- sample(colnames(seuratObj), size = round(ncol(seuratObj)*0.8))
      cl <- cl_full[s]
      
      ARI <- FMI <- NMI <- Pur <- c()
      
      seu = seuratObj[names(sort(statsvalue, decreasing = T)[1:top]),names(cl)]
      Seurat::VariableFeatures(seu) <- rownames(seu)
      seu = Seurat::ScaleData(seu)
      seu = Seurat::RunPCA(seu, features = rownames(seu))
      
      #Convert to SCE
      diet.seurat = Seurat::DietSeurat(seu, graphs = "pca")
      sce = Seurat::as.SingleCellExperiment(diet.seurat) 
      coord = seu@reductions$spatial@cell.embeddings
      colnames(coord) = c("row", "col")
      colData(sce) = cbind(colData(sce), coord) 
      mainExpName(sce) = "Spatial"
      
      set.seed(i)
      sce <- spatialPreprocess(sce, platform = "Visium",
                               skip.PCA = TRUE,
                               log.normalize = FALSE)
      set.seed(i)
      sce <- spatialCluster(sce,
                            q=length(unique(cl)), 
                            platform="Visium", d=np, gamma=2,
                            nrep=1000, burn.in=100)
      
      c1 <- list(cluster = sce$spatial.cluster)
      
      ARI <- c(ARI, mclust::adjustedRandIndex(c1$cluster, cl))
      FMI <- c(FMI, dendextend::FM_index(c1$cluster, cl))
      NMI <- c(NMI, igraph::compare(as.numeric(factor(c1$cluster)),  as.numeric(factor(cl)), method="nmi"))
      Pur <- c(Pur, IntNMF::ClusterPurity(c1$cluster, cl))
      
      df1 = data.frame(
        value = c(ARI, FMI, NMI, Pur),
        metric = rep(c("ARI", "FMI", "NMI", "Pur"), each = 1),
        method = method,
        top_n = top,
        clustering = "BayesSpace"
      )
      
      rm(diet.seurat, sce, seu)
      gc()
      return(df1)
      
    }, mc.cores = length(topList)))
    
    write.csv(df_top, file = here("data/bayesspace/stats/", paste0("bayesspace_", method, "_rep", i, ".csv")))
    
  })
  
}

f = list.files(path=here("data/bayesspace/stats/"))

df_method_bayesspace = do.call(rbind, lapply(f, function(x) {
  
  method = sapply(strsplit(x, "_"), "[[", 2)
  rep = gsub("rep|.csv", "", sapply(strsplit(x, "_"), "[[", 3))
  
  tmp_df = read.table(file = paste0(path, x),
                      sep = ",", header = TRUE, row.names = 1)
  
  return(tmp_df)
  
}))

save(here("data", "res_bayesspace.RData"))
