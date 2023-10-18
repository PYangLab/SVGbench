###################################
# Plot rankings of gene statistics after gene sampling
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(ggrastr)
library(glue)
library(here)

load(file=here("data", "gene_subsampling_stats_SVG_lists.RData"))

cor.list = list()
for (i in 1:length(subsampled_genes)) {
  
  cat(i, "............\n")
  sub = subsampled_genes[[i]]
  full = original_genes[[i]]
  
  cor.list[[i]] = sapply(1:length(all.names), function(j) {
    
    if (all.names[j] %in% names(sub)) {
      
      d = all.names[j]
      d2 = paste(unlist(strsplit(d, "_"))[1], "1", sep="_")
      cat(j, "............\n")
      g = names(sub[[d]])
      cor = cor(sub[[d]][g], full[[d2]][g], method="spearman", use="complete.obs")
      
    } else {
      cor = NA
    }
    
  })
  
}

mat = do.call(rbind, cor.list)
rownames(mat) = names(subsampled_genes)
colnames(mat) = names(subsampled_genes$gkm)

df = reshape2::melt(mat)
df = df[!is.na(df$value), ]

# plot results
pal = paletteer::palettes_d$ggthemes$Jewel_Bright[-3]
names(pal) = c("giotto_kmeans","giotto_rank", "seurat_moransi", "meringue", "nnSVG", "somde", "sparkx", "spatialde")

df$Var1 = factor(df$Var1, levels = c("giotto_kmeans","giotto_rank", "seurat_moransi", "meringue", "nnSVG", "somde", "sparkx", "spatialde"))

ggplot(df, aes(x=reorder(Var1, value), y=value, fill=Var1)) +
  geom_boxplot() +
  scale_fill_manual(values=pal) +
  theme_bw()

fn = here("plots", "gene_sampling_correlation_v2.pdf")
ggsave(filename=fn, width=6, height=4)
