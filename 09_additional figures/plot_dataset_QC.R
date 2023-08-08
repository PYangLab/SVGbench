###################################
# Plot gene expression QC across datasets
# Carissa Chen, updated Jul 2023
###################################

# load libraries
library(here)

load(here("additional_files", "dataset_QC_normCounts.RData"))

pzero.rows = lapply(normCounts, function(x) {
  pzero.rows= rowSums(x == 0)/ncol(x)
  names(pzero.rows) = rownames(x)
  return(pzero.rows)
})

avg.exprs.genes = lapply(normCounts, function(x) {
  avg.exprs = rowMeans2(x)
  names(avg.exprs) = rownames(x)
  return(avg.exprs)
})

# plot

pdf(here("plots", "dataset_QC.pdf"), height=10, width=6)
par(mfrow=(c(2,1)))
boxplot(pzero.rows, las=2, main="pzero_rows")
boxplot(avg.exprs.genes, las=2, main="avg expr of genes")
dev.off()

save(pzero.rows, avg.exprs.genes, file=here("additional_files", "dataset_QC_res.RData"))