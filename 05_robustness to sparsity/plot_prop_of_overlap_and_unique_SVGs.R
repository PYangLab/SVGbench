###################################
# Plot proportion of overlapping sig SVGs after spot subsampling
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(ggplot2)
library(here)

load(file=here("data", "spot_subsampling_pvalue_SVG_lists.RData"))

tmp = list()
count = 0
for (m in methods) {
  
  count = count + 1
  dat.sub = sub_pvalue[[grep(m, names(sub_pvalue))]]
  dat.nonsub = nonsub_pvalue[[m]]
  
  tmp[[count]] <- lapply(1:length(test_datasets), function(t) {
    p.sub = dat.sub[grep(test_datasets[t], names(dat.sub))]
    p.nonsub = dat.nonsub[[test_datasets[t]]]
    
    # calculate overlaps and unique SVGs
    tmp = lapply(p.sub, function(x) {
      sub.sig = names(x)[x < 0.05]
      nonsub.sig = names(p.nonsub)[p.nonsub < 0.05]
      n = intersect(sub.sig, nonsub.sig)
      n2 = setdiff(sub.sig, nonsub.sig)
      df = data.frame(overlap = length(n), unique_sub = length(n2))
      
    })
    
    tmp = do.call(rbind, tmp)
    return(tmp)
    
  })
  
  names(tmp[[count]]) = test_datasets
  tmp[[count]] = do.call(rbind, tmp[[count]])
  
}

names(tmp) = methods
prop.sig = do.call(rbind,tmp)

# plot proportions of unique and overlapping SVGs
prop.sig$rownames = rownames(prop.sig)
df = reshape2::melt(prop.sig)
df$method = sapply(strsplit(df$rownames, "\\."), "[[", 1)
df$dataset = sapply(strsplit(df$rownames, "\\."), "[[", 2)

pal = paletteer::palettes_d$ggthemes$Jewel_Bright[-3]
names(pal) = c("giotto_kmeans", "giotto_rank", "meringue", "moransi", "nnsvg", "somde", "spark", "spatialde")

# order by most overlapping
df[df$variable %in% c("overlap","unique_sub"), ] %>%
  ggplot(aes(y=value, x=method, fill=variable)) + 
  geom_bar(position="fill", stat="identity") +
  facet_wrap(vars(dataset), ncol=5) +
  ggtitle("Proportion of overlapping sig SVG and unique sig SVG") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))

fn = here("plots", "spot_sampling_sigSVGs_prop_plot.pdf")
ggsave(filename=fn, width=15, height=10)

