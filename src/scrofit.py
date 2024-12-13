from itertools import combinations
import os
import math
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix


import src.preprocessing as pp
from src.alignment import align_slides
import src.plotting as pl
from src.util import print_msg
from src.mapping import mapping_MCMF, mapping_1NN, adjusting_position


class SCRoFit(object):

    def __init__(self, adata_dict=None, out_dir='./', sample='sample'):
        self.adata_dict = adata_dict
        self.out_dir = out_dir
        self.sample = sample
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        for key, adata in adata_dict.items():
            if key == 'ST':
                pp.preprocess(adata, min_genes=100, min_cells=5)
            if key.startswith('SM'):
                pp.preprocess(adata, min_genes=100, min_cells=5)

    def check_slides(self, s=1):
        pl.plot_CCS(self, s=s)

    def flip_slides(self, key, direction='hv'):
        if 'v' in direction:
            max = self.adata_dict[key].obsm['spatial'][:, 1].max()
            self.adata_dict[key].obsm['spatial'][:, 1] = max - self.adata_dict[key].obsm['spatial'][:, 1]
        if 'v' in direction:
            max = self.adata_dict[key].obsm['spatial'][:, 0].max()
            self.adata_dict[key].obsm['spatial'][:, 0] = max - self.adata_dict[key].obsm['spatial'][:, 0]
        pl.plot_CCS(self)

    def align_slides(self, anchor_key='ST', method='rir', figsize=(10, 8),
                     echo=True, plot=True):
        if echo:
            print(f'COSPA aligning slides to {anchor_key} with {method}...')

        if 'spatial_' + method not in self.adata_dict[anchor_key].obsm.keys():
            align_slides(self.adata_dict, anchor_key=anchor_key, method=method)

        if plot:
            fig_fn = f'{self.out_dir}/{self.sample}_alignment.png'
            pl.plot_CCS(self, figsize=figsize, fig_fn=fig_fn)

    def mapping(self, source_key='SM', target_key='ST', mapping_method='MCMF',
                ccs_type='spatial_align', distance_method='euclidean', 
                n_thread=4, alpha=1, beta=1, k=3, n_batch=1000,
                adata_layer='log1p', verbose=True):
        X_adata = self.adata_dict[source_key]
        Y_adata = self.adata_dict[target_key]
        if mapping_method == '1NN':
            mapping_1NN(X_adata, Y_adata, ccs_type=ccs_type, distance_method=distance_method)        
        if mapping_method == 'MCMF':
            mapping_MCMF(X_adata, Y_adata, ccs_type=ccs_type, k=k,
                    n_thread=n_thread, n_batch=n_batch, alpha=alpha, beta=beta,
                    adata_layer=adata_layer, verbose=verbose)
            
        adjusting_position(X_adata, Y_adata, mapping_method=mapping_method)



