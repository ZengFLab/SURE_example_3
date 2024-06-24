import subprocess 
import glob
import numpy as np
import pandas as pd


from SURE import SURE
from utils.scdata_cached import setup_data_loader, SingleCellCached

import torch
from torch.utils.data import DataLoader

import datatable as dt


comp_ratio = 1.0 / 20 # 20 cells to 1 metacell
comp_lb = 200 # lower bound of the metacell size
comp_ub = 900 # upper bound of the metacell size

mtxFiles = glob.glob('/home/oem/SURE_example_3/human_brain_*_counts.txt.gz')
for i in np.arange(len(mtxFiles)):
    meta = pd.read_csv(f'/home/oem/SURE_example_3/human_brain_B{i+1}_cell.csv.gz')
    num_cells = meta.shape[0]
    
    num_metacells = int(num_cells * comp_ratio)
    num_metacells = 10 * int(num_metacells / 10)
    num_metacells = num_metacells if num_metacells > comp_lb else comp_lb
    num_metacells = num_metacells if num_metacells < comp_ub else comp_ub
    
    print('>>>')
    print(f'#{i+1} / {len(mtxFiles)}: /home/oem/SURE_example_3/human_brain_B{i+1}_counts.txt.gz')
    print(f'#cells: {num_cells}')
    print(f'#metacells: {num_metacells}')
    
    cmd = f'CUDA_VISIBLE_DEVICES=0  python SURE.py --data-file "/home/oem/SURE_example_3/human_brain_B{i+1}_counts.txt.gz" \
                        --undesired-factor-file "/home/oem/SURE_example_3/human_brain_B{i+1}_uwv.txt.gz" \
                        --jit \
                        -rt \
                        -lr 0.0001 \
                        -n 300 \
                        -de 50 \
                        -bs 1000 \
                        --seed 0 \
                        --cuda \
                        -cd {num_metacells} \
                        -zi exact \
                        -dirichlet \
                        -likeli negbinomial \
                        --save-model "human_brain_B{i+1}.pth" 2>&1 | tee "human_brain_B{i+1}.log"'
                        
    subprocess.call(cmd, shell=True)
    
    ModelPath = f'human_brain_B{i+1}.pth'
    Omic1Path=f'/home/oem/SURE_example_3/human_brain_B{i+1}_counts.txt.gz'
    ConditionPath2=None
    
    # load model
    model = torch.load(ModelPath)
    batch_size = 10000

    use_float64 = False
    use_cuda = True
    log_trans = False

    # load data
    data_cached = SingleCellCached(Omic1Path, ConditionPath2, log_trans=log_trans, use_cuda=False, use_float64 = use_float64)
    data_loader = DataLoader(data_cached, batch_size = batch_size, shuffle = False)
    
    cells = pd.read_csv(f'/home/oem/SURE_example_3/human_brain_B{i+1}_cell.csv.gz', header=0)
    genes = pd.read_csv(f'/home/oem/SURE_example_3/human_brain_B{i+1}_gene.csv.gz', header=0)
    
    assigns = []

    for xs1,ks2 in data_loader:
        if use_cuda:
            xs1 = xs1.cuda()

        zs = model.soft_assignments(xs1)
    
        if use_cuda:
            zs = zs.cpu().detach().numpy()
        else:
            zs = zs.detach().numpy()

        assigns.append(zs)

    assigns = np.concatenate(assigns, axis=0)
    df=pd.DataFrame(assigns, index=cells.iloc[:,0].values, columns=['MC_B{}_{}'.format(i+1,j+1) for j in np.arange(assigns.shape[1])])
    dt.Frame(df.reset_index()).to_csv(f'/home/oem/SURE_example_3/human_brain_B{i+1}_soft_assignments.txt.gz')

