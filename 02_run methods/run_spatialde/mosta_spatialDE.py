###################################
# call SVGs for E9.5 dataset using SpatialDE python
# Carissa Chen, updated Jul 2023
###################################

import pandas as pd
import NaiveDE
import SpatialDE
import numpy as np
import glob
import re
import random
import anndata

# create vector of filenames
inpath = '../../../SpatialData/MOSTA/'
datasets = [glob.glob(inpath + 'E9.5.MOSTA.h5ad')]

data_idx = re.sub('.MOSTA.h5ad', '', datasets)
data_idx = re.sub(inpath, '', datasets)
print(str('..................' + data_idx))
    
# read data
adata = anndata.read(name)

# filter
counts = adata.to_df(layer='count')
counts = counts.transpose()

# filter any genes in fewer than 30 cells
# filter any cells where more than 50% of counts are contributed by top 50 genes
counts = counts.loc[adata.var.n_cells > 30, adata.obs.pct_counts_in_top_50_genes < 50]

dup = (counts.index).duplicated()
counts = counts.loc[dup == False, ]

colzero = counts.sum(axis=0) != 0
rowzero = counts.sum(axis=1) != 0
counts = counts.loc[rowzero, colzero]

# get spatial coordinates
spatial_loc = pd.DataFrame(adata.obsm['spatial'], columns=['x_coord', 'y_coord'], index=adata.obs_names)
spatial_loc = spatial_loc.loc[counts.columns,]
spatial_loc['total_counts'] = counts.sum(axis=0)

# run spatialde

norm_expr = NaiveDE.stabilize(counts)
resid_expr = NaiveDE.regress_out(spatial_loc, norm_expr, 'np.log(total_counts)').T
result = SpatialDE.run(spatial_loc, resid_expr)

# write csv of results

outpath = '../SVG/MOSTA/spatialde/' + data_idx + '_svg_spatialde.csv'
result.to_csv(outpath)
