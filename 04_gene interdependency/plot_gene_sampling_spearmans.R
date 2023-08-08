###################################
# Plot rankings of gene statistics after gene sampling
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(ggrastr)
library(glue)
library(here)

load(file=here("data", "gene_subsampling_stats_SVG_lists.RData"))

methods = c("giotto_kmeans", "giotto_rank", "meringue", "moransi", "nnsvg", "somde", "spark", "spatialde")
test_datasets = c("Dataset21", "Dataset41")

pal = paletteer::palettes_d$ggthemes$Jewel_Bright[-3]
names(pal) = methods 

each_dataset = list()
count = 0
for (m in methods) {
  count = count + 1
  
  # subset for method
  dat.sub = sub_stats[[grep(m, names(sub_stats))]]
  dat.nonsub = nonsub_stats[[m]]
  
  # for each dataset
  each_dataset[[count]] <- lapply(1:length(test_datasets), function(t) {
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
      
      fn = here("plots", "gene_statistics_gene_subsampling", glue::glue(paste(m, names(sub)[[x]], sep="_"), ".pdf"))
      ggsave(filename=fn)
      
      
      # calculate correlation
      cor = cor(sub[[x]][g], nonsub[g], method="spearman")
    })
    
    names(cor.df) = names(sub)
    cor.df = reshape2::melt(cor.df)
    cor.df$method = m
    return(cor.df)
    
  })
  
  names(each_dataset[[count]]) = test_datasets
  each_dataset[[count]] = do.call(rbind, each_dataset[[count]])
  
}

names(each_dataset) = methods
plot_df = do.call(rbind, each_dataset)

# plot results
pal = paletteer::palettes_d$ggthemes$Jewel_Bright[-3]
names(pal) = methods

plot_df$dataset = sapply(strsplit(plot_df$L1, "_"), "[[", 1)

ggplot(plot_df, aes(x=method, y=value, fill=method)) +
geom_bar(position="dodge", stat="identity") +
labs(x=NULL, y="spearman's correlation", fill=NULL, colour=NULL) +
scale_fill_manual(values=pal) +
facet_wrap(vars(dataset), scales="free_x", ncol=3) +
theme_bw() +
theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1), aspect.ratio = 1)

fn = here("plots", "gene_sampling_correlation.pdf")
ggsave(filename=fn, width=10, height=10)
