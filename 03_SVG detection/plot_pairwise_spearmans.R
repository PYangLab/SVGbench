###################################
# Script to generate pairwise spearmans correlation boxplots
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(ggplot2)
library(ggthemes)
library(paletteer)
library(here)

sp.mat = readRDS(here("data", "res_pairwise_spearmans_corr.rds"))

# colour by method
pal = paletteer::palettes_d$ggthemes$Jewel_Bright[-3]
names(pal) = c("giotto_kmeans", "giotto_rank", "meringue","moransi", "nnsvg", "somde", "spark", "spatialde")

ggplot(sp.mat, aes(x=Var2, y=value, fill=Var2)) +
  geom_boxplot(width=0.7, outlier.shape=NA) +
  geom_jitter(stroke = 0, alpha=0.8, width=0.25) +
  labs(x=NULL, y="spearman's correlation", fill=NULL) +
  scale_fill_manual(values=pal) +
  facet_wrap(~Var1, scales="free_x", ncol=4) +
  ylim(c(-0.3,1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1), aspect.ratio = 1)
ggsave(here("plots", "correlation_boxplot_colour_by_method.pdf"), width=8, height=8)

# colour by technology
ggplot(sp.mat, aes(x=Var2, y=value)) +
  geom_boxplot(width=0.7, color="grey", outlier.shape=NA) +
  geom_jitter(aes(color=platform), stroke = 0, alpha=0.8, width=0.25) +
  labs(x=NULL, y="spearman's correlation", fill=NULL) +
  scale_color_tableau(palette = "Tableau 10") +
  facet_wrap(~Var1, scales="free_x", ncol=4) +
  ylim(c(-0.3,1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1), aspect.ratio = 1)
ggsave(here("plots", "correlation_boxplot_colour_by_platform.pdf"), width=8, height=8)

# colour by binned spot numbers
br = cut(sp.mat$spots, breaks=c(0,500,1000,2000,9000))
sp.mat$spots_bin = br

ggplot(sp.mat, aes(x=Var2, y=value)) +
  geom_boxplot(width=0.7, color="grey", outlier.shape=NA) +
  geom_jitter(aes(color=spots_bin), stroke = 0, alpha=0.8, width=0.25) +
  labs(x=NULL, y="spearman's correlation", fill=NULL) +
  scale_color_brewer(palette = "RdBu", direction = -1) +
  facet_wrap(~Var1, scales="free_x", ncol=4) +
  ylim(c(-0.3,1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1), aspect.ratio = 1)
ggsave(here("plots", "correlation_boxplot_colour_by_spots.pdf"), width=8, height=8)
