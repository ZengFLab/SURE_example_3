import subprocess 
import glob
import pandas as pd
import scanpy as sc

import SEACells as sea 

import datatable as dt


cmd='rm *_B*metacells.h5ad'
subprocess.call(cmd, shell=True)

# load h5ad
adata_ = sc.read_h5ad('b9171f05-8112-4a55-95f2-4cf8a57df8a2.h5ad')

# process each batch
files = glob.glob('*_B*_soft_assignments.txt.gz')
n_batches = len(files)

for i in range(n_batches):
    print(f'{i+1} / {n_batches}')

    fl = files[i]

    A=dt.fread(fl, header=True).to_pandas().set_index('index')
    A.columns = [f'human_brain_B{i+1}_{j}' for j in A.columns]
    metacells = A.idxmax(axis=1)

    adata = adata_[metacells.index.tolist()].copy()
    adata.obs['metacells'] = metacells[adata.obs_names]

    adata_metacells = sea.core.summarize_by_SEACell(adata, SEACells_label='metacells', celltype_label='cell_type', summarize_layer='X')

    adata_metacells.obs['cell_type'] = adata.obs.groupby('metacells').apply(
        lambda x: x['supercluster_term'].mode().iloc[0]
        ).loc[adata_metacells.obs_names]

    adata_metacells.obs['tissue'] = adata.obs.groupby('metacells').apply(
        lambda x: x['tissue'].mode().iloc[0]
        ).loc[adata_metacells.obs_names]

    adata_metacells.obs['Batch'] = f'B{i+1}'
    
    adata_metacells.write_h5ad(f'human_brain_B{i+1}_metacells.h5ad')
    
    
# combine all
files = glob.glob('*_B*_metacells.h5ad')

adata_list = []

for fl in files:
    tmp = sc.read_h5ad(fl)
    adata_list.append(tmp)
    
adata = sc.concat(adata_list)
adata.write_h5ad('human_brain_primary_metacells.h5ad')