###################################
# Plot rankings of gene statistics after cell sampling
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(ggrastr)
library(glue)
library(here)

load(file=here("data", "spot_subsampling_stats_SVG_lists.RData"))

# plot scatter plots
methods = c("giotto_kmeans", "giotto_rank", "meringue", "moransi", "nnsvg", "somde", "spark", "spatialde")
test_datasets = c("Dataset21", "Dataset41")

pal = paletteer::palettes_d$ggthemes$Jewel_Bright[-3]
names(pal) = methods

tmp = list()
count = 0
for (m in methods) {
  
  count = count + 1
  dat.sub = sub_stats[[grep(m, names(sub_stats))]]
  dat.nonsub = nonsub_stats[[m]]
  
  tmp[[count]] <- lapply(1:length(test_datasets), function(t) {
    sub = dat.sub[grep(test_datasets[t], names(dat.sub))]
    nonsub = dat.nonsub[[test_datasets[t]]]
    
    cor.df = lapply(1:length(sub), function(x) {
      g = names(sub[[x]])
      
      
      # plot scatterplots
      df = data.frame(nonsub = rank(nonsub[g]), sub = rank(sub[[x]][g]))
      
      ggplot(df, aes(x=nonsub, y=sub)) +
        geom_point_rast(color = pal[m], alpha=0.3) +
        theme_base() +
        theme(aspect.ratio = 1)
      
      fn = here("plots", "gene_statistics_spot_subsampling",  glue::glue(paste(m, names(sub)[[x]], sep="_"), ".pdf"))
      ggsave(filename=fn)
       
       
      # calculate correlation
      cor = cor(sub[[x]][g], nonsub[g], method="spearman")
      
    })
    
    names(cor.df) = names(sub)
    cor.df = reshape2::melt(cor.df)
    cor.df$method = m
    return(cor.df)
    
  })
  
  names(tmp[[count]]) = test_datasets
  tmp[[count]] = do.call(rbind, tmp[[count]])
  
}

names(tmp) = methods
plot_df = do.call(rbind,tmp)

# generate barplot of spearmans correlation

plot_df$dataset = sapply(strsplit(plot_df$L1, "_"), "[[", 1)

ggplot(plot_df, aes(x=dataset, y=value, fill=method)) +
geom_bar(position="dodge", stat="identity") +
labs(x=NULL, y="spearman's correlation", fill=NULL, colour=NULL) +
scale_fill_manual(values=pal) +
facet_wrap(vars(dataset), scales="free_x", ncol=5) +
theme_bw() +
theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1), aspect.ratio = 1)

fn = here("plots", "spot_sampling_spearmans_correlation.pdf")
ggsave(filename=fn, width=10, height=10)