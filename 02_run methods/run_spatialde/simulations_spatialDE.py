###################################
# call SVGs for simulation datasets using SpatialDE python
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
import os

# create vector of filenames
inpath = '../../../SpatialData/Simulations'
outpath = '../SVG/simulations/call_variable_genes/spatialde/'

count_filenames = glob.glob(inpath + '*/Spatial_count.csv')
meta_filenames = glob.glob(inpath + '*/Locations.csv')

random.seed(1)
for i in range(0,10):

    name = count_filenames[i]
    data_idx = re.sub(inpath, '', name)
    print(str(i) + '..................' + data_idx)
    
# read data
    counts = pd.read_csv(count_filenames[i], sep =",", index_col = 0).T
    spatial_loc = pd.read_csv(meta_filenames[i], sep =",", index_col = 0)
    
# filter
    dup = (counts.index).duplicated()
    counts = counts.loc[dup == False, ]
    
    colzero = counts.sum(axis=0) != 0
    rowzero = counts.sum(axis=1) != 0
    counts = counts.loc[rowzero, colzero]
    
    spatial_loc = spatial_loc.loc[counts.columns, ]
    spatial_loc['total_counts'] = counts.sum(axis=0)
    
# run spatialde

    time_start = time.perf_counter()
    tracemalloc.start()
    norm_expr = NaiveDE.stabilize(counts) # columns are samples, rows are genes
    resid_expr = NaiveDE.regress_out(spatial_loc, norm_expr, 'np.log(total_counts)').T
    result = SpatialDE.run(spatial_loc, resid_expr)
 
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    time_elapsed = (time.perf_counter() - time_start)
    tracemalloc.clear_traces()
    
    peak = peak/1024.0/1024.0

# write csv of computational times
    csv_path = str(outpath + data_idx + '_time_spatialde.csv')

    with open(csv_path, 'w+') as output:
        columns = ['time', 'n_cell', 'n_gene', 'dataset', 'Peak_RAM_Used_MiB', 'method']
        output.write(','.join(columns) + '\n')
        n_cell = len(counts.columns)
        n_gene = len(counts)
        dataset = data_idx
        method = 'spatialde'
          
        comp = [str(time_elapsed), str(n_cell), str(n_gene), dataset, str(peak), method]
        output.write(','.join(comp))
    
        output.write('\n')
    
# write csv of results

    outpath1 = '../SVG/simulations/call_variable_genes/spatialde/' + data_idx + '_svg_spatialde.csv'
    result.to_csv(outpath1)
