###################################
# Plot boxplots of proportion of ground truth SVGs among significant SVGs
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(ggplot2)
library(here)

load(here("data", "res_simulation_calculate_accuracy.RData"))
load(file=here("data", "res_simulations_prop_sigSVGs.RData"))

methods = c("giotto_kmeans", "giotto_rank", "seurat_moransi", "nnSVG", "somde",  "sparkx", "spatialde", "meringue")
sample_keep = c("Dataset21","Dataset22","Dataset30",
                "Dataset32","Dataset34","Dataset37",
                "Dataset38","DBiTseq_0628cL","DBiTseq_E11-2L")

tmp = goldtruth[dataset_o]
tmp = lapply(tmp, function(x) x$gene[x$category == TRUE])

tmp_df = do.call(rbind, lapply(methods, function(method_name) {
  
  tmp_method = lapply(names(tmp), function(x) {
    
    tmp_subset = tmp_results %>%
      filter(sample == x) %>%
      filter(method == method_name)
    return(unique(as.character(tmp_subset$gene)))
    
  })
  names(tmp_method) = names(tmp)
  
  result_ratio = sapply(names(tmp), function(x) {
    res = sum(tmp_method[[x]] %in% tmp[[x]])/length(tmp_method[[x]])
    if (is.nan(res)) {
      return(0)
    } else {
      return(res)
    }
  })
  df = data.frame(
    dataaset = names(result_ratio),
    results = result_ratio,
    method = method_name,
    sample = names(tmp)
  )
  return(df)
  
}))

tmp_df2 = do.call(rbind, lapply(methods, function(method_name) {
  
  tmp_method = lapply(names(tmp), function(x) {
    
    tmp_subset = tmp_results2 %>%
      filter(sample == x) %>%
      filter(method == method_name)
    return(unique(as.character(tmp_subset$gene)))
    
  })
  names(tmp_method) = names(tmp)
  
  result_ratio = sapply(names(tmp), function(x) {
    res = sum(tmp_method[[x]] %in% tmp[[x]])/length(tmp_method[[x]])
    if (is.nan(res)) {
      return(0)
    } else {
      return(res)
    }
  })
  df = data.frame(
    dataaset = names(result_ratio),
    results = result_ratio,
    method = method_name,
    sample = names(tmp)
  )
  return(df)
  
}))

# plot proportion plots

method_col = paletteer::palettes_d$ggthemes$Jewel_Bright[-3]
names(method_col) = sort(methods)

p1 = ggplot(tmp_df[tmp_df$sample %in% sample_keep,], aes(x = reorder(method, results, median), y = results, fill = "grey")) +
  geom_boxplot() + 
  geom_jitter() + 
  scale_fill_manual(values = method_col) + ggtitle("0.05") + coord_cartesian(ylim = c(0.2, 1)) + 
  ylab("Proportion of ground truth SVG among significant SVGs") + theme_pubclean() 
p2 = ggplot(tmp_df2[tmp_df2$sample %in% sample_keep,], aes(x = reorder(method, results, median), y = results, fill = "grey")) +
  geom_boxplot() + 
  geom_jitter() + 
  scale_fill_manual(values = method_col) + ggtitle("0.01") +  coord_cartesian(ylim = c(0.2, 1)) + 
  ylab("Proportion of ground truth SVG among significant SVGs") + theme_pubclean()

p = patchwork::wrap_plots(list(p1,p2))

fn = here("plots", "simulations_accuracy_prop_boxplots.pdf")
p %>% ggsave(filename = fn, width = 5, height = 5)
