#################
# Plot masks of datasets used for simulation
# Last updated: Carissa Chen Jul 2023
#################

# load libraries
library(ggplot2)
library(tidyverse)
library(here)

path="../Simulation/test_scDesign3/"
f = list.files(path)

f_tmp = gsub(".RData", "", f)
f = f[f_tmp%in% c("Dataset21","Dataset22","Dataset30",
                  "Dataset32","Dataset34","Dataset37",
                  "Dataset38","DBiTseq_0628cL","DBiTseq_E11-2L")]

lapply(1:length(f), function(x) {
  
  print(x)
  load(paste0(path, f[[x]]), verbose = TRUE)
  print(dim(simu_sce))
  
  loc = colData(simu_sce)[,c("spatial1","spatial2")]
  loc$nGene = MatrixGenerics::colSums2(counts(simu_sce)>0)
  loc$nUMI = MatrixGenerics::colSums2(counts(simu_sce))
  long = as.data.frame(loc)
  
  long$x = loc[,1]
  long$y = loc[,2]
  
  p0 = as_tibble(long, rownames = "Cell") %>% 
    ggplot(aes(x = x, y = y, color = nGene)) +
    ggrastr::geom_point_rast(size = 0.1) +
    scale_colour_gradientn(colors = viridis_pal(option = "magma")(10)) + 
    coord_fixed(ratio = 1) + theme(axis.text.x = element_text(angle = 45))
  
  p1 = as_tibble(long, rownames = "Cell") %>% 
    ggplot(aes(x = x, y = y, color = nGene)) +
    ggrastr::geom_point_rast(size = 0.2) +
    scale_colour_gradientn(colors = viridis_pal(option = "magma")(10)) + 
    coord_fixed(ratio = 1) + theme(axis.text.x = element_text(angle = 45))
  
  p2 = as_tibble(long, rownames = "Cell") %>% 
    ggplot(aes(x = x, y = y, color = nGene)) +
    ggrastr::geom_point_rast(size = 0.4) +
    scale_colour_gradientn(colors = viridis_pal(option = "magma")(10)) + 
    coord_fixed(ratio = 1) + theme(axis.text.x = element_text(angle = 45))
  
  p3 = as_tibble(long, rownames = "Cell") %>% 
    ggplot(aes(x = x, y = y, color = nGene)) +
    ggrastr::geom_point_rast(size = 1) +
    scale_colour_gradientn(colors = viridis_pal(option = "magma")(10)) + 
    coord_fixed(ratio = 1) + theme(axis.text.x = element_text(angle = 45))
  
  p4 = as_tibble(long, rownames = "Cell") %>% 
    ggplot(aes(x = x, y = y, color = nGene)) +
    ggrastr::geom_point_rast(size = 2) +
    scale_colour_gradientn(colors = viridis_pal(option = "magma")(10)) + 
    coord_fixed(ratio = 1) + theme(axis.text.x = element_text(angle = 45))
  
  p5 = as_tibble(long, rownames = "Cell") %>% 
    ggplot(aes(x = x, y = y, color = nGene)) +
    ggrastr::geom_point_rast(size = 5) +
    scale_colour_gradientn(colors = viridis_pal(option = "magma")(10)) + 
    coord_fixed(ratio = 1) + theme(axis.text.x = element_text(angle = 45))
  
  p = patchwork::wrap_plots(list(p0,p1,p2,p3,p4,p5), ncol = 3) + patchwork::plot_layout(guides = "collect")
  
  fn = here("plots", "simulated_masks", gsub(".RData", "", f[[x]]), ".pdf")
  ggsave(p, filename = fn, width = 10, height = 5)
  
})
