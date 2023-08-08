###################################
# Run hierarchical and kmeans clustering and evaluate concordance of labels against ground truth
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(kernlab)
library(parallel)
library(patchwork)
library(tidyverse)
library(here)

load(file=here("data", "res_clustering.RData"))

topList = seq(100,2000,100)

df_method = do.call(rbind, lapply(1:length(statsList), function(method_idx) {
  
  print(method_idx)
  np = 20; rep = 15
  
  statsvalue = statsList[[method_idx]]
  pvalue = pvalueList[[method_idx]]
  method = names(statsList)[[method_idx]]
  
  df_top = do.call(rbind, lapply(topList, function(top) {
    
    clList = lapply(1:rep, function(i) {
      set.seed(i)
      s <- sample(colnames(exp.log), size = round(ncol(exp.log)*0.8))
      cl_subset <- cl_full[s]
      return(cl_subset)
    })
    
    ARI <- FMI <- NMI <- Pur <- c()
    for(i in 1:rep) {
      
      set.seed(i)
      cl = clList[[i]]
      
      s <- sample(colnames(exp.log), size = round(ncol(exp.log)*0.8))
      pca <- scater::calculatePCA(exp.log[names(sort(statsvalue, decreasing = T)[1:top]), s], ncomponents = np)
      
      hc = hclust(dist(pca))
      c1 <- list()
      c1$cluster <- cutree(hc, k  = length(unique(cl)))
      
      ARI <- c(ARI, mclust::adjustedRandIndex(c1$cluster, cl))
      FMI <- c(FMI, dendextend::FM_index(c1$cluster, cl))
      NMI <- c(NMI, igraph::compare(as.numeric(factor(c1$cluster)),  as.numeric(factor(cl)), method="nmi"))
      Pur <- c(Pur, IntNMF::ClusterPurity(c1$cluster, cl))
      
    }
    
    df1 = data.frame(
      value = c(ARI, FMI, NMI, Pur),
      metric = rep(c("ARI", "FMI", "NMI", "Pur"), each = rep),
      method = method,
      top_n = top,
      clustering = "pca_kmeans",
      cluster_num = NA
    )
    
    ARI2 <- FMI2 <- NMI2 <- Pur2 <- c()
    for(i in 1:rep) {
      set.seed(i)
      
      cl = clList[[i]]
      
      s <- sample(colnames(exp.log), size = round(ncol(exp.log)*0.8))
      pca <- scater::calculatePCA(exp.log[names(sort(statsvalue, decreasing = T)[1:top]), s], ncomponents = np)
      
      hc = hclust(dist(pca))
      c1 <- list()
      c1$cluster <- cutree(hc, k  = length(unique(cl)))
      
      ARI2 <- c(ARI2, mclust::adjustedRandIndex(c1$cluster, cl))
      FMI2 <- c(FMI2, dendextend::FM_index(c1$cluster, cl))
      NMI2 <- c(NMI2, igraph::compare(as.numeric(factor(c1$cluster)),  as.numeric(factor(cl)), method="nmi"))
      Pur2 <- c(Pur2, IntNMF::ClusterPurity(c1$cluster, cl))
      
    }
    
    df2 = data.frame(
      value = c(ARI2, FMI2, NMI2, Pur2),
      metric = rep(c("ARI", "FMI", "NMI", "Pur"), each = rep),
      method = method,
      top_n = top,
      clustering = "hc",
      cluster_num = NA
    )
    
    df = rbind(df1,df2)
    
    return(df)
    
  }))
  return(df_top)
  
}, mc.cores = 2))

# write.csv(df_method, file=here("data", "res_clustering_results.csv"))



# plot boxplots by top SVGs
df_method = read.csv(here("data", "res_clustering_results.csv"))

is.odd <- function(v) v %% 2 != 0
idx = is.odd(1:length(unique(df_method$top_n)))

col = RColorBrewer::brewer.pal(8, "Spectral")
col = colorRampPalette(col)(length(unique(df_method$top_n)[idx]))


# plot kmeans clustering
p1 = df_method %>%
  filter(top_n %in% unique(df_method$top_n)[idx]) %>%
  filter(metric == "ARI" & clustering == "pca_kmeans") %>%
  ggplot(aes(x = reorder(method, value, mean), y = value, fill = as.factor(top_n))) + 
  geom_boxplot(notch = FALSE) +  ggtitle("ARI") +  scale_fill_manual(values = col) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90))
p2 = df_method %>%
  filter(top_n %in% unique(df_method$top_n)[idx]) %>%
  filter(metric == "NMI" & clustering == "pca_kmeans") %>%
  ggplot(aes(x = reorder(method, value, mean), y = value, fill = as.factor(top_n))) + 
  geom_boxplot(notch = FALSE) + ggtitle("NMI") +  scale_fill_manual(values = col) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90))
p3 = df_method %>%
  filter(top_n %in% unique(df_method$top_n)[idx]) %>%
  filter(metric == "FMI" & clustering == "pca_kmeans") %>%
  ggplot(aes(x = reorder(method, value, mean), y = value, fill = as.factor(top_n))) + 
  geom_boxplot(notch = FALSE) + ggtitle("FMI") +  scale_fill_manual(values = col) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90))
p4 = df_method %>%
  filter(top_n %in% unique(df_method$top_n)[idx]) %>%
  filter(metric == "Pur" & clustering == "pca_kmeans") %>%
  ggplot(aes(x = reorder(method, value, mean), y = value, fill = as.factor(top_n))) + 
  geom_boxplot(notch = FALSE) + ggtitle("Purity") +  scale_fill_manual(values = col) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90))
p = patchwork::wrap_plots(list(p1,p2,p3,p4), ncol = 2) + patchwork::plot_layout(guides = "collect")

fn = here("plots", "kmeans_clustering_concordance.pdf")
ggsave(p, filename = fn, width = 15, height = 10)


# plot hierarchical clustering
p1 = df_method %>%
  filter(top_n %in% unique(df_method$top_n)[idx]) %>%
  filter(metric == "ARI" & clustering == "hc") %>%
  ggplot(aes(x = reorder(method, value, mean), y = value, fill = as.factor(top_n))) + 
  geom_boxplot(notch = FALSE) +  ggtitle("ARI") +  scale_fill_manual(values = col) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90))
p2 = df_method %>%
  filter(top_n %in% unique(df_method$top_n)[idx]) %>%
  filter(metric == "NMI" & clustering == "hc") %>%
  ggplot(aes(x = reorder(method, value, mean), y = value, fill = as.factor(top_n))) + 
  geom_boxplot(notch = FALSE) + ggtitle("NMI") +  scale_fill_manual(values = col) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90))
p3 = df_method %>%
  filter(top_n %in% unique(df_method$top_n)[idx]) %>%
  filter(metric == "FMI" & clustering == "hc") %>%
  ggplot(aes(x = reorder(method, value, mean), y = value, fill = as.factor(top_n))) + 
  geom_boxplot(notch = FALSE) + ggtitle("FMI") +  scale_fill_manual(values = col) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90))
p4 = df_method %>%
  filter(top_n %in% unique(df_method$top_n)[idx]) %>%
  filter(metric == "Pur" & clustering == "hc") %>%
  ggplot(aes(x = reorder(method, value, mean), y = value, fill = as.factor(top_n))) + 
  geom_boxplot(notch = FALSE) + ggtitle("Purity") +  scale_fill_manual(values = col) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90))
p = patchwork::wrap_plots(list(p1,p2,p3,p4), ncol = 2) + patchwork::plot_layout(guides = "collect")

fn = here("plots", "hierarchical_clustering_concordance.pdf")
ggsave(p, filename = fn, width = 15, height = 10)