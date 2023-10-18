###################################
# plot clustering concordance for all methods
# Carissa Chen, updated Oct 2023
###################################

library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(here)

# load bayesspace, spaGCN, SINFONIA, kmeans and hierarchical clustering results
load(file=here("data", "res_bayesspace.RData"))
load(file=here("data", "res_spaGCN.RData"))
load(file=here("data", "res_SINFONIA.RData"))
df_method_old = read.csv(here("data", "res_clustering_results.csv"))

df_method = rbind(df_method_bayesspace,
                  df_method_sinfonia,
                  df_method_spagcn,
                  df_method_old)

is.odd <- function(v) v %% 2 != 0
idx = is.odd(1:length(unique(df_method$top_n)))

df_method$top_n = as.factor(as.numeric(df_method$top_n))
df_method$method = factor(df_method$method)
levels(df_method$method) = c("giotto_kmeans","giotto_rank","meringue","moransi","nnsvg","somde","spark","spatialde")
col = RColorBrewer::brewer.pal(8, "Spectral")
col = colorRampPalette(col)(length(unique(df_method$top_n)[idx]))

df_method$clustering = factor(df_method$clustering, levels=c("BayesSpace", "sinfonia_leiden","sinfonia_louvain","SpaGCN","pca_kmeans", "hc"))

metric_idx = "ARI"
p1 = df_method %>%
  filter(top_n %in% unique(df_method$top_n)[idx] & metric == "ARI") %>%
  ggplot(aes(x = method, y = value, fill = as.factor(top_n))) + 
  geom_boxplot() +
  ggtitle("ARI") +  scale_fill_manual(values = col) + 
  facet_wrap(~clustering, ncol = 3, scales = "free") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90))
p1 |> ggsave(filename = here("plots", paste0("clustering_", metric_idx, ".pdf")), width = 25, height = 10)

metric_idx = "NMI"
p1 = df_method %>%
  filter(top_n %in% unique(df_method$top_n)[idx] & metric == metric_idx) %>%
  ggplot(aes(x = method, y = value, fill = as.factor(top_n))) + 
  geom_boxplot() +
  ggtitle(metric_idx) +  scale_fill_manual(values = col) + 
  facet_wrap(~clustering, ncol = 3, scales = "free") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90))
p1 |> ggsave(filename = here("plots", paste0("clustering_", metric_idx, ".pdf")), width = 25, height = 10)

metric_idx = "FMI"
p1 = df_method %>%
  filter(top_n %in% unique(df_method$top_n)[idx] & metric == metric_idx) %>%
  ggplot(aes(x = method, y = value, fill = as.factor(top_n))) + 
  geom_boxplot() +
  ggtitle(metric_idx) +  scale_fill_manual(values = col) + 
  facet_wrap(~clustering, ncol = 3, scales = "free") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90))
p1 |> ggsave(filename = here("plots", paste0("clustering_", metric_idx, ".pdf")), width = 25, height = 10)

metric_idx = "Pur"
p1 = df_method %>%
  filter(top_n %in% unique(df_method$top_n)[idx] & metric == metric_idx) %>%
  ggplot(aes(x = method, y = value, fill = as.factor(top_n))) + 
  geom_boxplot() +
  ggtitle(metric_idx) +  scale_fill_manual(values = col) + 
  facet_wrap(~clustering, ncol = 3, scales = "free") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90))
p1 |> ggsave(filename = here("plots", paste0("clustering_", metric_idx, ".pdf")), width = 25, height = 10)

