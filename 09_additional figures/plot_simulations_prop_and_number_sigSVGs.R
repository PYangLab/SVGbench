###################################
# Plot proportion and number of SVGs ddetected under different fdr thresholds
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(ggplot2)
library(here)

load(here("data", "res_simulation_calculate_accuracy.RData"))

raw_results = results_method
raw_results$value = raw_results$pvalue

sample_keep = c("Dataset21","Dataset22","Dataset30",
                "Dataset32","Dataset34","Dataset37",
                "Dataset38","DBiTseq_0628cL","DBiTseq_E11-2L")


tbl_all <- do.call(rbind, lapply(sample_keep, function(i) {
  
  tmp = raw_results %>%
    filter(sample == i)
  
  gene_index = unique(tmp$gene)
  svg.num = do.call(rbind, lapply(methods, function(method_name) {
    
    tmp_mat = tmp[tmp$method == method_name, ]
    tmp_mat = tmp_mat[,"value"]
    
    res = c(sum(na.omit(as.numeric(tmp_mat)) == 0),
            sum(na.omit(as.numeric(tmp_mat)) <= 1e-10 &
                  na.omit(as.numeric(tmp_mat)) > 0),
            sum(na.omit(as.numeric(tmp_mat)) <= 0.01 &
                  na.omit(as.numeric(tmp_mat)) > 1e-10),
            sum(na.omit(as.numeric(tmp_mat)) <= 0.05 &
                  na.omit(as.numeric(tmp_mat)) > 0.01),
            sum(na.omit(as.numeric(tmp_mat) > 0.05)))
    
    df= data.frame(
      count = res,
      type = c("zero", "below 1e-10", "below 0.01", "below 0.05", "nonsig"),
      method = method_name,
      sample = i
    )
    return(df)
    
    
  }))
  return(svg.num)
}))

result <- tbl_all

methods_rank = c("giotto_rank", "giotto_kmeans", "seurat_moransi","meringue",
                 "spatialde", "nnSVG", "somde", "sparkx")
result$method = factor(result$method, levels = methods_rank)
result$type= factor(result$type, levels = c("zero", "below 1e-10", "below 0.01", "below 0.05", "nonsig"))

my_col = RColorBrewer::brewer.pal(6, "Spectral")
names(my_col) = cutoff
my_col = my_col[3:5]
my_col = c(my_col, "#EF6F6A", "#6388B4")
names(my_col) = c( "below 1e-10", "below 0.01", "below 0.05", "zero", "nonsig")

g = ggplot(result, aes(x = sample, y = count, fill = type)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = my_col) + 
  facet_wrap(~method, ncol = 8) + 
  theme_classic() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90)) 

fn = here("plots", "simulation_prop_sigSVGs.pdf")
g %>% ggsave(filename = fn, width = 15, height = 5)

g = result %>%
  filter(type != "nonsig") %>%
  ggplot(aes(x = sample, y = count, fill = type)) + 
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(values = my_col) + 
  facet_wrap(~method, ncol = 4) + 
  theme_classic() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90)) 

fn = here("plots", "simulation_number_sigSVGs.pdf")
g %>% ggsave(filename = fn, width = 12, height = 5)
