###################################
# Script to subsample spots and call SVGs using SOMDE
# Carissa Chen, updated Oct 2023
###################################

import pandas as pd
import numpy as np
import string

import tqdm
import tracemalloc
from somde import SomNode
from somde.util import *
import anndata
import glob
import re
import time
import resource
import random

# create vector of filenames
inpath = '../../../SpatialData/SpatialBenchmark/DataUploadSubset/'

count_filenames = glob.glob(inpath + '*/Spatial_count.csv')
meta_filenames = glob.glob(inpath + '*/Locations.csv')

cp = 0.8

random.seed(2022)
for i in range(0,2):

    name = matched_count_files[i]
    data_idx = re.sub(inpath, '', name)
    print(str(i) + '..................' + data_idx)
  
  # read data
    counts = pd.read_csv(count_files[i], sep =",", index_col = 0).T
    spatial_loc = pd.read_csv(meta_files[i], sep =",", index_col = 0)
  
  # filter
    dup = (counts.index).duplicated()
    counts = counts.loc[dup == False, ]
  
    colzero = counts.sum(axis=0) != 0
    rowzero = counts.sum(axis=1) != 0
    counts = counts.loc[rowzero, colzero]
  
    spatial_loc = spatial_loc.loc[counts.columns, ]
  
  # subsample percentage of spots
    counts = counts.sample(frac=cp,axis=1, replace=False)
    spatial_loc = spatial_loc.loc[counts.columns, ].values.astype(np.float32)
  
  # run somde
  
    som = SomNode(spatial_loc,10)
    ndf,ninfo = som.mtx(counts)
    nres = som.norm()
    result, SVnum = som.run()
  
  # write csv of results
  
    outpath = '../SVG/subsample_spots/somde/' + data_idx + '_svg_somde_subsample_spots_' + str(cp) + '.csv'
    result.to_csv(outpath)
