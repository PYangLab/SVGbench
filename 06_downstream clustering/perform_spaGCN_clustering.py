import os,csv,re
import pandas as pd
import numpy as np
import math
from scipy.sparse import issparse
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import cv2
import scanpy as sc
import SpaGCN as spg
import random, torch
from scipy.stats import rankdata
import re

adata=sc.read("../../../SpatialData/MOSTA/E9.5_E1S1.MOSTA.h5ad")

r_seed=t_seed=n_seed=1000

#Normalize and take log for UMI
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

#Set coordinates
adata.obs["spatial"]=adata.obsm['spatial'][:,0]
adata.obs["spatial"]=adata.obsm['spatial'][:,1]
adata.obs["x_pixel"]=adata.obsm['spatial'][:,0]
adata.obs["y_pixel"]=adata.obsm['spatial'][:,1]
adata.obs["x_array"]=adata.obsm['spatial'][:,0]
adata.obs["y_array"]=adata.obsm['spatial'][:,1]

path_stats = "../../../SpaGCN/stats/"
f_celltypes = [fn for fn in os.listdir(path_stats)
              if any(fn.endswith(ext) for ext in ".csv")]
f_tissue = []
for i in range(0, len(f_celltypes)):
    f_tissue.append(f_celltypes[i].split("_")[1])
    
def unique(list1):
    x = np.array(list1)
    print(np.unique(x))
    return(np.unique(x))

methods = []
for i in range(len(f_tissue)):
    tmp = f_tissue[i]
    tmp = re.sub('.csv', '', tmp)
    methods.append(tmp)
    
methods

folder_path_tosave = "../../../SpaGCN/SpaGCN_results/"

reps = range(0,10)
tops = range(100, 2000+1, 100)

for method_idx in methods:
    
    print(method_idx)
    
    np = 20
    rep = 10
    
    statsvalue = pd.read_csv(path_stats + "stats_" + method_idx + ".csv", index_col = 0)
    method = method_idx
    
    for top in tops:
        
        print(top)
        genes_to_keep = statsvalue.iloc[:top,0]

        for rep_idx in reps:
            
            print(rep_idx)
            import numpy as np
            np.random.seed(rep_idx)
            sampled_columns = adata.obs_names.to_list()
            sample_size = int(round(adata.n_obs * 0.8))
            s = np.random.choice(adata.obs_names, size=sample_size, replace=False)
            
            adata_subset = adata[s,genes_to_keep]
            cl_subset = adata.obs['annotation']
            
            x_array=adata_subset.obs["x_array"].tolist()
            y_array=adata_subset.obs["y_array"].tolist()
            x_pixel=adata_subset.obs["x_pixel"].tolist()
            y_pixel=adata_subset.obs["y_pixel"].tolist()
            
            adj_no_img=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)
            p=0.5
            l=spg.search_l(p, adj_no_img, start=0.01, end=1000, tol=0.01, max_run=100)
            
            n_clusters = len(np.unique(adata_subset.obs["annotation"]))
            res=spg.search_res(adata_subset,
                               adj_no_img,
                               l,
                               n_clusters,
                               start=0.7,
                               step=0.1,
                               tol=5e-3,
                               lr=0.05,
                               max_epochs=20,
                               r_seed=r_seed,
                               t_seed=t_seed,
                               n_seed=n_seed)
            
            clf=spg.SpaGCN()
            clf.set_l(l)
            #Set seed
            import random
            random.seed(r_seed)
            torch.manual_seed(t_seed)
            np.random.seed(n_seed)
                        
            #Run
            clf.train(
                adata_subset,
                adj_no_img,
                init_spa=True,
                init="louvain",
                res=res,
                tol=5e-3,
                lr=0.05,
                max_epochs=200
            )
            
            y_pred, prob=clf.predict()
            adata_subset.obs["pred"]= y_pred
            adata_subset.obs["pred"]= adata_subset.obs["pred"].astype('category')
            
            #Do cluster refinement(optional)
            #shape="hexagon" for Visium data, "square" for ST data.
            adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
            
            refined_pred=spg.refine(sample_id=adata_subset.obs.index.tolist(), pred=adata_subset.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
            adata_subset.obs["refined_pred"]=refined_pred
            adata_subset.obs["refined_pred"]=adata_subset.obs["refined_pred"].astype('category')
            
            data = {
                'cells' : adata_subset.obs_names,
                'label' : adata_subset.obs["annotation"],
                'new_label' : adata_subset.obs["refined_pred"]
            }
            pd.DataFrame(data).to_csv(folder_path_tosave + '{}_{}_{}.csv'.format(method_idx, top, rep_idx))

