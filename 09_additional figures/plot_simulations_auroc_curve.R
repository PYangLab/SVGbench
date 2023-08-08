###################################
# Calculate AUC curve
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(dplyr)
library(yardstick)
library(ggpubr)
library(ROSE)
library(here)

load(here("data", "res_simulation_calculate_accuracy.RData"))

raw_results = results_method
stats_here = "p-value"
sample_keep = c("Dataset21","Dataset22","Dataset30",
                "Dataset32","Dataset34","Dataset37",
                "Dataset38","DBiTseq_0628cL","DBiTseq_E11-2L")
methods = c("giotto_kmeans", "giotto_rank", "seurat_moransi", "nnSVG", "somde",  "sparkx", "spatialde", "meringue")

tmp_variable = unique(paste(raw_results$sample, 
                            raw_results$method, 
                            sep = "___"))
test = data.frame(
  sample = sapply(strsplit(tmp_variable, "___"), "[[", 1),
  method = sapply(strsplit(tmp_variable, "___"), "[[", 2)
)
test = test[test$sample %in% sample_keep,]

method_col = paletteer::palettes_d$ggthemes$Jewel_Bright[-3]
names(method_col) = sort(methods)

tbl_all <- lapply(1:nrow(test), function(i) {
  
  tmp = raw_results %>%
    filter(sample == test[i, "sample"], 
           method == test[i, "method"])
  
  gt = goldtruth[[test[i, "sample"]]]
  
  all <- tmp$gene
  gt <- gt$gene[gt$category == TRUE]
  truth = ifelse(all %in% gt, "SVG", "Other") %>% factor()
  
  if (stats_here == "p-value") {
    results = -log10(tmp$pvalue)
  } else {
    results = tmp$statistics
  }
  names(results)= tmp$gene
  
  pred =  rank(-results, na.last = TRUE)
  
  
  method = rep(test[i, "method"], length(truth))
  sample = rep(test[i, "sample"], length(truth))
  
  tbl = tibble(truth, pred, sample, method) 
  
  return(tbl)
  
})

tbl_all <- do.call(rbind, tbl_all)

p1 = tbl_all %>% 
  group_by(sample, method) %>%
  yardstick::roc_curve(truth = truth, pred) %>% 
  ggplot(aes(x = 1 - specificity, y = sensitivity,
             colour = method)) +
  geom_path(linewidth = 2, show.legend = TRUE) +
  geom_abline(lty = 1, col="red") +
  scale_color_manual(values = method_col) +
  facet_wrap(~sample, ncol = 3, nrow = 3) + 
  coord_equal() + theme_classic() + 
  theme(
    legend.position = "bottom"
  )

fn = here("plots", "simulations_auc_curve.pdf")
p1 %>% ggsave(filename = fn, width = 20, height = 20)
