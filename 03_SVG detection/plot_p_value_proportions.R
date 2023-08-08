###################################
# Script to plot p-value proportional plot
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(ggplot2)
library(patchwork)
library(paletteer)
library(dplyr)

result = readRDS(here("data", "res_p_value_thres.rds"))

# order each method by count of < 0.05 genes

# giotto kmeans
o = result %>%
  group_by(type, dataset) %>%
  filter(method == "giotto_kmeans", type == 0.05) %>%
  arrange(desc(count)) %>%
  ungroup() %>%
  pull(dataset)

result$dataset = factor(result$dataset, levels=o)  

g1 = result %>%
  filter(method == "giotto_kmeans") %>%
  ggplot(aes(x = dataset, y = count, fill = type)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_tableau(palette="Superfishel Stone", direction=-1) +
  theme_classic() + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1))

# giotto rank
o = result %>%
  group_by(type, dataset) %>%
  filter(method == "giotto_rank", type == 0.05) %>%
  arrange(desc(count)) %>%
  ungroup() %>%
  pull(dataset)

result$dataset = factor(result$dataset, levels=o)  

g2 = result %>%
  filter(method == "giotto_rank") %>%
  ggplot(aes(x = dataset, y = count, fill = type)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_tableau(palette="Superfishel Stone", direction=-1) +
  theme_classic() + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1))

# meringue
o = result %>%
  group_by(type, dataset) %>%
  filter(method == "meringue", type == 0.05) %>%
  arrange(desc(count)) %>%
  ungroup() %>%
  pull(dataset)

result$dataset = factor(result$dataset, levels=o)  

g3 = result %>%
  filter(method == "meringue") %>%
  ggplot(aes(x = dataset, y = count, fill = type)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_tableau(palette="Superfishel Stone", direction=-1) +
  theme_classic() + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1))

# moransi
o = result %>%
  group_by(type, dataset) %>%
  filter(method == "moransi", type == 0.05) %>%
  arrange(desc(count)) %>%
  ungroup() %>%
  pull(dataset)

result$dataset = factor(result$dataset, levels=o)  

g4 = result %>%
  filter(method == "moransi") %>%
  ggplot(aes(x = dataset, y = count, fill = type)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_tableau(palette="Superfishel Stone", direction=-1) +
  theme_classic() + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1))

# nnsvg
o = result %>%
  group_by(type, dataset) %>%
  filter(method == "nnsvg", type == 0.05) %>%
  arrange(desc(count)) %>%
  ungroup() %>%
  pull(dataset)

result$dataset = factor(result$dataset, levels=o)  

g5 = result %>%
  filter(method == "nnsvg") %>%
  ggplot(aes(x = dataset, y = count, fill = type)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_tableau(palette="Superfishel Stone", direction=-1) +
  theme_classic() + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1))

# somde
o = result %>%
  group_by(type, dataset) %>%
  filter(method == "somde", type == 0.05) %>%
  arrange(desc(count)) %>%
  ungroup() %>%
  pull(dataset)

result$dataset = factor(result$dataset, levels=o)  

g6 = result %>%
  filter(method == "somde") %>%
  ggplot(aes(x = dataset, y = count, fill = type)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_tableau(palette="Superfishel Stone", direction=-1) +
  theme_classic() + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1))

# sparkx
o = result %>%
  group_by(type, dataset) %>%
  filter(method == "spark", type == 0.05) %>%
  arrange(desc(count)) %>%
  ungroup() %>%
  pull(dataset)

result$dataset = factor(result$dataset, levels=o)  

g7 = result %>%
  filter(method == "spark") %>%
  ggplot(aes(x = dataset, y = count, fill = type)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_tableau(palette="Superfishel Stone", direction=-1) +
  theme_classic() + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1))

# spatialde
o = result %>%
  group_by(type, dataset) %>%
  filter(method == "spatialde", type == 0.05) %>%
  arrange(desc(count)) %>%
  ungroup() %>%
  pull(dataset)

result$dataset = factor(result$dataset, levels=o) 

g8 = result %>%
  filter(method == "spatialde") %>%
  ggplot(aes(x = dataset, y = count, fill = type)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_tableau(palette="Superfishel Stone", direction=-1) +
  theme_classic() + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle=45, vjust=1, hjust=1))

glist = list(g1,g2,g3,g4,g5,g6,g7,g8)

patch = patchwork::wrap_plots(glist, ncol=4)

fn = here("plots", "p_value_prop_plot_ordered.pdf")
patch %>% ggsave(filename=fn, width=20, height=15)
