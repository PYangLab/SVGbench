###################################
# Script to generate gene statistic dependency plots
# Carissa Chen, updated Oct 2023
###################################

# load libraries
library(pheatmap)
library(here)

csv = read.csv(here("data", "dataset_metadata.csv"), head=T)
load(here("data", "res_dataset_QC.RData"))

# load output from SVG detection
svg_lists = readRDS(here("data", "input_SVG_lists_by_dataset.rds"))

methods = c("giotto_kmeans", "giotto_rank", "meringue", "moransi", "nnsvg", "somde", "spark", "spatialde")

index = list(giotto_kmeans= list(method="giotto_kmeans", statistic = "adj.p.value", gene="genes"),
             giotto_rank = list(method="giotto_rank", statistic = "adj.p.value", gene="genes"), 
             meringue =  list(method="meringue", statistic = "observed", gene=NA),
             moransi =  list(method="moransi", statistic = "observed", gene=NA),
             nnsvg =  list(method="nnsvg", statistic = "LR_stat", gene="symbol"),
             somde =  list(method="somde", statistic = "LLR", gene="g"),
             spark = list(method="spark", statistic = "adjustedPval", gene=NA),
             spatialde =  list(method="spatialde", statistic = "LLR", gene="g"))

# create matrix of correlation
m1 <- matrix(NA, 8, length(svg_lists))
m2 <- matrix(NA, 8, length(svg_lists))

# get gene statistics for each method
for (i in 1:length(svg_lists)) {
  gene_stats_list = lapply(1:length(svg_lists[[i]]), function(x) {
    tb = svg_lists[[i]][[x]]
    method = names(svg_lists[[i]])[x]
    statistic_col = index[[method]]$statistic
    
    if(is.na(index[[method]]$gene) == TRUE) {
      g = rownames(tb)
    } else {
      g = tb[, index[[method]]$gene]
    }
    statistic = tb[, statistic_col]
    
    if (statistic_col %in% c("adj.p.value", "adjustedPval")) {
      statistic = -log10(statistic)
    }
    names(statistic) = g
    return(statistic)
  })
  names(gene_stats_list) = names(svg_lists[[i]])
  
  # by average expression
  exprs = avg.exprs.genes[[names(svg_lists)[i]]]
  g = names(exprs)
  
  for(j in 1:length(methods)) {
    cat(j, ".......\n")
    cb <- cbind(exprs[g], gene_stats_list[[j]][g])
    cb.complete <- cb[complete.cases(cb),]
    m1[j, i] <- cor(cb.complete[,1], cb.complete[,2], method = "spearman")
  }
  
  # by percentage of zeros
  pzero = pzero.rows[[names(svg_lists)[i]]]
  g = names(pzero)
  
  for(j in 1:length(methods)) {
    cat(j, ".......\n")
    cb <- cbind(pzero[g], gene_stats_list[[j]][g])
    cb.complete <- cb[complete.cases(cb),]
    m2[j, i] <- cor(cb.complete[,1], cb.complete[,2], method = "spearman")
  }
  
}

rownames(m1) <- rownames(m2) <- methods
colnames(m1) <- colnames(m2) <- names(svg_lists)


# plot heatmaps

# plot dependency by average expression
mean = sort(Matrix::rowMeans(m1), decreasing=T)
fn = here("plots", "avg_expr_stat_boxplot.pdf")
pdf(fn)
boxplot(m1[names(mean),])
dev.off()


# plot dependency by percentage of zeros
mean = sort(Matrix::rowMeans(m2), decreasing=T)
fn = here("plots", "pzero_stat_cor.pdf")
pdf(fn)
pheatmap::pheatmap(m2[names(mean),], cluster_rows = F, cluster_cols = F)
dev.off()

