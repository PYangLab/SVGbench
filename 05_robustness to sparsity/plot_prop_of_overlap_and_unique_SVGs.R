###################################
# Plot proportion of overlapping sig SVGs after spot subsampling
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(ggplot2)
library(here)

load(file=here("data", "spot_subsampling_pvalue_SVG_lists.RData"))

all.names = names(subsampled_cells_pval[[8]])

overlap.list = list()
for (i in 1:length(subsampled_cells_pval)) {
  
  cat(i, "............\n")
  sub = subsampled_cells_pval[[i]]
  full = original_cells_pval[[i]]
  
  overlap.list[[i]] = sapply(1:length(all.names), function(j) {
    
    if (all.names[j] %in% names(sub)) {
      
      d = all.names[j]
      d2 = paste(unlist(strsplit(d, "_"))[1], "1", sep="_")
      
      sub = na.omit(sub[[d]])
      full = na.omit(full[[d2]])
      
      g1 = names(sub)[sub < 0.05]
      g2 = names(full)[full < 0.05]
      
      ifelse(length(setdiff(g1, g2)) > 0, {
        o = length(setdiff(g1, g2))/length(g1)
      }, {o = 0})
      
    } else {
      o = NA
    }
    
  })
}

mat = do.call(rbind, overlap.list)
rownames(mat) = names(subsampled_cells_pval)
colnames(mat) = all.names

df = reshape2::melt(mat)
df = df[!is.na(df$value), ]

pal = paletteer::palettes_d$ggthemes$Jewel_Bright[-3]
methods = c("gkm", "grank", "meringue", "moransi", "nnsvg", "somde", "sparkx", "spatialde")
names(pal) = methods

df$Var1 = factor(df$Var1, levels = c("gkm", "grank", "meringue", "moransi", "nnsvg", "somde", "sparkx", "spatialde"))

ggplot(df, aes(x=reorder(Var1, -value), y=value, fill=Var1)) +
  geom_boxplot() +
  scale_fill_manual(values=pal) +
  labs(y="proportion of unique SVGs") +
  theme_bw()

fn = here("plots", "spot_sampling_sigSVGs_prop_plot.pdf")
ggsave(filename=fn, width=6, height=4)
