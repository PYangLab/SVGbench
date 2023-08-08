###################################
# Calculate accuracy of each method
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(dplyr)
library(yardstick)
library(ggpubr)
library(ROSE)
library(cutpointr)
library(here)

load(here("data", "res_simulation_calculate_accuracy.RData"))

raw_results = results_method
raw_results$value = raw_results$pvalue

tmp_variable = unique(paste(raw_results$sample, 
                            raw_results$method, 
                            sep = "___"))
test = data.frame(
  sample = sapply(strsplit(tmp_variable, "___"), "[[", 1),
  method = sapply(strsplit(tmp_variable, "___"), "[[", 2)
)

sample_keep = c("Dataset21","Dataset22","Dataset30",
                "Dataset32","Dataset34","Dataset37",
                "Dataset38","DBiTseq_0628cL","DBiTseq_E11-2L")
test = test[test$sample %in% sample_keep,]

cutoff=c(1e-100,1e-50,1e-10, 0.01, 0.05, 0.1) 

tbl_cutoff <- lapply(cutoff, function(threshold) {
  
  tbl_all <- lapply(1:nrow(test), function(i) {
    
    tmp = raw_results %>%
      filter(sample == test[i, "sample"], 
             method == test[i, "method"])
    
    gt = goldtruth[[test[i, "sample"]]]
    
    all <- tmp$gene
    gt <- gt$gene[gt$category == TRUE]
    other <- all[!all %in% gt]
    
    truth = ifelse(all %in% gt, "SVG", "Other") %>% factor()
    results = tmp$value
    names(results)= tmp$gene
    
    pred =  na.omit(names(results)[results < threshold])
    pred_other = na.omit(names(results)[results > threshold])
    
    tp = length(pred[pred %in% gt])
    fp = length(pred[!pred %in% gt])
    fn = length(pred_other[pred_other %in% gt])
    tn = length(pred_other[pred_other %in% other])
    
    tpr = cutpointr::tpr(tp, fn)
    fpr = cutpointr::tpr(fp, tn)
    tnr = cutpointr::tpr(fp, tn)
    fnr = cutpointr::tpr(tp, fn)
    accuracy = cutpointr::accuracy(tp, fp, tn, fn)
    fdr = cutpointr::false_discovery_rate(tp, fp, tn, fn)
    f1 = cutpointr::F1_score(tp, fp, tn, fn)
    
    method = test[i, "method"]
    sample = test[i, "sample"]
    
    tbl = tibble(tpr, fpr, tnr, fnr,
                 accuracy, fdr, f1,
                 sample, method, threshold) 
    
    return(tbl)
    
  })
  tbl_all <- do.call(rbind, tbl_all)
  return(tbl_all)
  
})

tbl_cutoff <- do.call(rbind, tbl_cutoff)
dim(tbl_cutoff)

tbl_cutoff2 = tbl_cutoff #%>% 
tbl_cutoff2$threshold_exaggerated = tbl_cutoff2$threshold < 0.01

my_col = c("#f8766d", "#c49a00", "#54b600", "#00c094", "#00b6eb", "#a58aff", "#fe63da")
my_col = c("#05A8AA", "#B8D5B8", "#D7B49E", "#DC602E", "#BC412B")
my_col = colorRampPalette(my_col)(6)
my_col = RColorBrewer::brewer.pal(6, "Spectral")
names(my_col) = cutoff

methods_rank = c("giotto_rank", "giotto_kmeans", "seurat_moransi","meringue",
                 "spatialde", "nnSVG", "somde", "sparkx")
tbl_cutoff2$method = factor(tbl_cutoff2$method, levels = methods_rank)
tbl_cutoff2 = tbl_cutoff2[order(as.numeric(tbl_cutoff2$threshold), decreasing = FALSE),]

p1 <- tbl_cutoff2 %>% 
  ggplot(aes(x = fdr, y = tpr, 
             color = as.factor(threshold))) +
  geom_point(alpha = 0.8, show.legend = TRUE, stroke = NA, size = 4) +
  geom_vline(xintercept = c(0.01, 0.05, 0.1), col = c("red"), lty="dotted") + 
  scale_color_manual(values = my_col) +
  facet_wrap(~threshold_exaggerated+method, ncol = 8, nrow = 2) + 
  coord_equal() +
  theme_classic() + 
  theme(
    legend.position = "bottom"
  )

fn = here("plots", "tpr_fdr_accuracy.pdf")
p1 %>% ggsave(filename = fn, width = 15, height = 7)
