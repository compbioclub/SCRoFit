import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap, to_rgba
from scipy import stats

from src.util import get_array


def get_value_boundary(v):
    if type(v) == pd.DataFrame:
        v = v.to_numpy()
    if abs(v.max()) > abs(v.min()):
        return abs(v.max())
    else:
        return abs(v.min())

def plt_util(title):
    plt.xticks([])
    plt.yticks([])
    plt.title(title)
    plt.colorbar()
    


def customize_cmap(cmap, position='center', color='black'):
    if cmap is None:
        base_cmap = customize_exp_cmap()
    else:
        base_cmap = plt.get_cmap(cmap)
    # Define a custom colormap with the value 0 set to black
    colors = [base_cmap(i) for i in range(base_cmap.N)]

    if position == 'center':
        colors[base_cmap.N // 2] = to_rgba(color)  # Set the middle color
    if position == 'min':
        colors[0] = to_rgba(color)
    if position == 'max':
        colors[-1] = to_rgba(color)
    custom_cmap = LinearSegmentedColormap.from_list('custom', colors, N=base_cmap.N)
    return custom_cmap


def customize_exp_cmap():
    cmap_name = 'viridis_r'
    base_cmap = plt.get_cmap(cmap_name)

    truncate_ratio = 0.7
    colors = base_cmap(np.linspace(0, truncate_ratio, 256))[::-1]
    custom_cmap = LinearSegmentedColormap.from_list('custom', colors)
    return custom_cmap


def plot_CCS(obj, s=1, ncol=3, nrow=None, figsize=(10, 10), fig_fn=None):

    adata_dict = obj.adata_dict
    adata = list(adata_dict.values())[1]

    fig = plt.figure(figsize=figsize)
    nsub = len(adata.obsm.keys()) * len(adata_dict)
    if nrow is None:
        nrow = int(np.ceil(nsub / ncol))

    i = 1
    for key, adata in adata_dict.items():
        plt.subplot(nrow, ncol, i)
        points = adata.obsm['spatial']
        plt.scatter(points[:, 0], points[:, 1], s=s)
        plt.title(f'OCS: {key}\n Size: {points.shape[0]}')
        i += 1

    for ccs_type in adata.obsm.keys():
        if not ccs_type.startswith('spatial_'):
            continue
        if ccs_type.endswith('_map'):
            continue        
        plt.subplot(nrow, ncol, i)
        for key, adata in adata_dict.items():
            if ccs_type in adata.obsm:
                points = adata.obsm[ccs_type]
                plt.scatter(points[:, 0], points[:, 1], label=key, s=s)
                plt.legend()
        plt.title(f'CCS: {ccs_type}')
        i += 1

    for ccs_type in adata.obsm.keys():
        if not ccs_type.startswith('spatial_'):
            continue
        if not ccs_type.endswith('_map'):
            continue        
        plt.subplot(nrow, ncol, i)
        for key, adata in adata_dict.items():
            if ccs_type in adata.obsm:
                points = adata.obsm[ccs_type]
                plt.scatter(points[:, 0], points[:, 1], label=key, s=s)
                plt.legend()
        plt.title(f'OCS: {ccs_type}\n Size: {points.shape[0]}')
        i += 1

    #plt.tight_layout()
    plt.suptitle(obj.sample)
    if fig_fn is not None:
        fig.savefig(fig_fn)


def _spatial_scatter(adata, color, ccs_type='spatial',
                    s=5, cmap='viridis'):

    coords = adata.obsm[ccs_type]
    plt.scatter(coords[:,0], coords[:,1], c=adata[:,color].X.toarray(), s=s, cmap=cmap)
    plt_util(color)