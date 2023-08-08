###################################
# Plot computational memory
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(ggplot2)
library(here)

memory = readRDS(here("data", "res_memory_calculation.rds"))
csv = read.csv(here("data", "dataset_metadata.csv"), head=T)

memory$L1 = factor(memory$L1, levels=csv$dataset)
memory$method = factor(memory$method, levels=c("giotto_kmeans", "giotto_rank","meringue", "seurat_moransi", "nnSVG", "somde", "spark", "spatialde"))

pal = paletteer::palettes_d$ggthemes$Jewel_Bright[-3]

ggplot(memory.df, aes(x=L1, y=value/1000, fill=method)) +
  geom_col(position="dodge") +
  labs(y="memory usage GiB", x=NULL) +
  scale_fill_manual(values=pal) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, vjust=1, hjust=1))

fn = here("plots", "memory_barplots.pdf")
ggsave(fn, width=6, height=6)

# create bubble plot of rankings

# rank each method per dataset
sub = list()
count = 0
for (d in csv$dataset) {
  count = count + 1
  sub[[count]] = memory[which(memory$L1 == d), ]
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

fn = here("plots", "memory_bubblerank.pdf")
ggsave(filename=fn, width=8, height=6)