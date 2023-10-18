###################################
# Script to spots spots and call SVGs using SpatialDE python
# Carissa Chen, updated Jul 2023
###################################

import pandas as pd
import NaiveDE
import SpatialDE
import numpy as np
import tracemalloc
import glob
import re
import time
import resource
import random
import csv

# create vector of filenames
inpath = '../../../SpatialData/SpatialBenchmark/DataUploadSubset/'
count_filenames = glob.glob(inpath + '*/Spatial_count.csv')
meta_filenames = glob.glob(inpath + '*/Locations.csv')
        
cp = [0.8]

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
    
# spots percentage of spots
    counts = counts.sample(frac=cp,axis=1, replace=False)
    spatial_loc = spatial_loc.loc[counts.columns, ]
    spatial_loc['total_counts'] = counts.sum(axis=0)
    
# run spatialde

    norm_expr = NaiveDE.stabilize(counts)
    resid_expr = NaiveDE.regress_out(spatial_loc, norm_expr, 'np.log(total_counts)').T
    result = SpatialDE.run(spatial_loc, resid_expr)

# write csv of results

    outpath = '../SVG/subsample_spots/spatialde/' + data_idx + '_svg_spatialde_subsample_spots_' + str(cp) + '.csv'
    result.to_csv(outpath)
