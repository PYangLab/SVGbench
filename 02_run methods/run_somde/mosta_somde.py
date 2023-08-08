###################################
# call SVGs for E9.5 dataset using SOMDE
# Carissa Chen, updated Jul 2023
###################################

import pandas as pd
import numpy as np
import string

import tqdm
from somde import SomNode
from somde.util import *
import anndata
import glob
import re
import random

inpath = '../../../SpatialData/MOSTA/""
outpath = '../SVG/MOSTA/somde/'
datasets = str(inpath + 'E9.5.MOSTA.h5ad')

data_idx = re.sub('.MOSTA.h5ad', '', datasets)
data_idx = re.sub(inpath, '', datasets)
print(str('..................' + data_idx))

# read data
adata = anndata.read(filename=datasets)

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
coord = pd.DataFrame(adata.obsm['spatial'], columns=['x_coord', 'y_coord'], index=adata.obs_names)
coord = coord.loc[counts.columns,]
spatial_loc=coord[['x_coord','y_coord']].values.astype(np.float32)

# run somde

som = SomNode(spatial_loc,10)
ndf,ninfo = som.mtx(counts)
nres = som.norm()
result, SVnum = som.run()

# write csv of results
result.to_csv(outpath + data_idx + '_svg_somde.csv')
