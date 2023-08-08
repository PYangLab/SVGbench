###################################
# Plot computational time
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(ggplot2)
library(ggbreak)
library(here)

time = readRDS(here("data", "res_time_calculation.rds"))
csv = read.csv(here("data", "dataset_metadata.csv"), head=T)

time$L1 = factor(time$L1, levels=csv$dataset)
time$method = factor(time$method, levels=c("giotto_kmeans", "giotto_rank","meringue", "seurat_moransi", "nnSVG", "somde", "spark", "spatialde"))

pal = paletteer::palettes_d$ggthemes$Jewel_Bright[-3]

ggplot(time, aes(x=L1, y=value/60, fill=method)) +
  geom_col(position="dodge") +
  labs(y="minutes", x=NULL) +
  scale_fill_manual(values=pal) +
  scale_y_break(c(150,600)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, vjust=1, hjust=1))

fn = here("plots", "time_barplots.pdf")
ggsave(fn, width=6, height=6)

# create bubble plot of rankings

# rank each method per dataset
sub = list()
count = 0
for (d in csv$dataset) {
  count = count + 1
  sub[[count]] = time[which(time$L1 == d), ]
  sub[[count]]$rank = rank(-sub[[count]]$value)
}
rank.df = do.call(rbind, sub)
rank.df$L1 = factor(rank.df$L1, levels=csv$dataset)
rank.df$method = factor(rank.df$method, levels=rev(c("giotto_kmeans", "giotto_rank", "meringue", "seurat_moransi", "nnSVG", "somde", "spark", "spatialde")))

ggplot(rank.df, aes(x=L1, y=reorder(method, rank), fill=as.character(rank), size=rank)) +
  geom_point(pch=21) +
  labs(y="rank", x=NULL) +
  scale_fill_brewer(palette = "RdYlBu", direction=-1) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, vjust=1, hjust=1))

fn = here("plots", "time_bubblerank.pdf")
ggsave(fn, width=8, height=6)


