###################################
# call SVGs for benchmarking datasets using SOMDE
# Carissa Chen, updated Jul 2023
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
outpath = '../SVG/simulations/call_variable_genes/somde/'

count_filenames = glob.glob(inpath + '*/Spatial_count.csv')
meta_filenames = glob.glob(inpath + '*/Locations.csv')

random.seed(1)
for i in range(0,22):

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
    
    spatial_loc = spatial_loc.loc[counts.columns, ].values.astype(np.float32)
    
# run somde

    time_start = time.perf_counter()
    tracemalloc.start()
    
    som = SomNode(spatial_loc,10)
    ndf,ninfo = som.mtx(counts)
    nres = som.norm()
    result, SVnum = som.run()
    
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    time_elapsed = (time.perf_counter() - time_start)
    tracemalloc.clear_traces()
    
    peak = peak/1024.0/1024.0

# write csv of computational times
    csv_path = str(outpath + data_idx + '_time_somde.csv')

    with open(csv_path, 'w+') as output:
        columns = ['time', 'n_cell', 'n_gene', 'dataset', 'Peak_RAM_Used_MiB', 'method']
        output.write(','.join(columns) + '\n')
        n_cell = len(counts.columns)
        n_gene = len(counts)
        dataset = data_idx
        method = 'somde'
          
        comp = [str(time_elapsed), str(n_cell), str(n_gene), dataset, str(peak), method]
        output.write(','.join(comp))
    
        output.write('\n')
    
# write csv of results
    outpath1 = outpath + data_idx + '_svg_somde.csv'
    result.to_csv(outpath1)
