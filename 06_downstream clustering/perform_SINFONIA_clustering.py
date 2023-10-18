import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import sinfonia
import warnings
import re

adata=sc.read("../../../SpatialData/MOSTA/E9.5_E1S1.MOSTA.h5ad")

r_seed=t_seed=n_seed=1000

#Normalize and take log for UMI
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

path_stats = "../../../SINFONIA/stats/"
f_celltypes = [fn for fn in os.listdir(path_stats)
              if any(fn.endswith(ext) for ext in ".csv")]
f_tissue = []
for i in range(0, len(f_celltypes)):
    f_tissue.append(f_celltypes[i].split("_")[1])

methods = []
for i in range(len(f_tissue)):
    tmp = f_tissue[i]
    tmp = re.sub('.csv', '', tmp)
    methods.append(tmp)

reps = range(0,10)
folder_path_tosave = "../../../SINFONIA/sinfonia_clustering_results/"
tops = range(100, 2000+1, 100)

for method_idx in methods:
    
    print(method_idx)
    
    statsvalue = pd.read_csv(path_stats + "stats_" + method_idx + ".csv", index_col = 0)
    method = method_idx
    
    for top in tops:
        
        print(top)
        genes_to_keep = statsvalue.iloc[:top,0]

        for rep_idx in reps:
            
            print(rep_idx)
            np.random.seed(rep_idx)
            sampled_columns = adata.obs_names.to_list()
            sample_size = int(round(adata.n_obs * 0.8))
            s = np.random.choice(adata.obs_names, size=sample_size, replace=False)
            
            adata_subset = adata[s,genes_to_keep]
            cl_subset = adata.obs['annotation']
            
            sc.pp.pca(adata_subset)
            sc.pp.neighbors(adata_subset)
            sc.tl.umap(adata_subset)
            
            sc.tl.louvain(adata_subset, key_added="default_louvain")
            sc.tl.leiden(adata_subset, key_added="default_leiden")
            
            #Set seed
            import random
            random.seed(r_seed)
            np.random.seed(n_seed)
            adata_subset = sinfonia.get_N_clusters(adata_subset, n_cluster=adata_subset.obs['annotation'].nunique(), cluster_method='louvain')
            adata_subset = sinfonia.get_N_clusters(adata_subset, n_cluster=adata_subset.obs['annotation'].nunique(), cluster_method='leiden')
            
            data = {
                'cells' : adata_subset.obs_names,
                'label' : adata_subset.obs['annotation'],
                'louvain_label' : adata_subset.obs['louvain'],
                'leiden_label' : adata_subset.obs['leiden']
            }
            pd.DataFrame(data).to_csv(folder_path_tosave + '{}_{}_{}.csv'.format(method_idx, top, rep_idx))
