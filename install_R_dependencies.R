install_if_missing <- function(packages, installer = install.packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    installer(missing)
  }
}

cran_packages <- c(
  "corrplot",
  "cutpointr",
  "dplyr",
  "ggbreak",
  "ggplot2",
  "ggpubr",
  "ggrastr",
  "ggthemes",
  "glue",
  "here",
  "jaccard",
  "kernlab",
  "paletteer",
  "patchwork",
  "pheatmap",
  "RColorBrewer",
  "reshape2",
  "ROSE",
  "scales",
  "tidyr",
  "tidyverse",
  "viridis",
  "yardstick"
)

install_if_missing(cran_packages)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "BiocParallel",
  "DESeq2",
  "S4Vectors",
  "scran",
  "scuttle",
  "SingleCellExperiment",
  "SpatialExperiment",
  "tradeSeq",
  "zellkonverter"
)

install_if_missing(bioc_packages, function(pkgs) {
  BiocManager::install(pkgs, ask = FALSE, update = FALSE)
})

message("Core R dependency installation complete.")
message("Some method-specific packages, including Giotto, MERINGUE, SPARK, ")
message("scDesign3 and PseudotimeDE, may require upstream installation steps.")
