# -*- coding: utf-8 -*-
import os 
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math

import hotspot
import matplotlib.pyplot as plt
import matplotlib.colors
import seaborn as sns
import mplscience

from scipy.io import mmread
from scipy.sparse import issparse, csr_matrix

import pickle
import warnings; warnings.simplefilter('ignore')

def create_adata_sub(
        exp_path,
        adata_path,
        gender
    ):
        """Creat the pseudobulk annadta form group difference expression data and compute the PCA latent(loadings) 
        and compute the cell PCA score use the groupDE PCA loadings 

        The resulting adata_sub containing the PCA latet of the pseudobulk group DE gene  

        Parameters
        ----------
        exp_path: the path of the group_DE_adjust_expression data (reg:subtype/0.pseudobulk/rawcount/DE/one_vs_one)
            Count matrix (shape is genes by sample_id )
        adata_path: the .h5ad file (reg: Analysis/h5adFiles)
        gender: M or F

        """
        df = pd.read_table(exp_path,
                   sep = '\t', index_col=0)
        adata_pb = sc.AnnData(pd.DataFrame(
                        df.T,  
                        index=df.columns,  
                        columns=df.index,  
                        ))
        metadata = pd.read_table('/data/yangyu/Project/skinAtlas/70samples/sample_info.txt',sep='\t',index_col=0) 
        cell_names = adata_pb.obs_names
        adata_pb.obs = metadata.loc[cell_names]
        #sc.pp.pca(adata_pb,n_comps=50)
        n_features = adata_pb.shape[1]  
        n_samples = adata_pb.shape[0]   
        adata_pb.obs['celltype_Granular'] = [exp_path.split('_DE')[0].split('/')[-1]] * n_samples
        n_components = min(n_samples, n_features)  
        print(n_components)
        if n_components >= 50 :
            sc.pp.pca(adata_pb,n_comps=50)
        else:
            n_components -= 1
            sc.pp.pca(adata_pb,n_comps=n_components)
            print(n_components)

        adata = sc.read_h5ad(adata_path)
        
        adata_sub = adata[:,adata_pb.var_names]
        if gender == 'M':
             adata_sub = adata_sub[adata_sub.obs['sample_group'].isin(['NE','NT','PS'])]
        else:
             adata_sub = adata_sub[adata_sub.obs['sample_group'].isin(['NE','NF','NP','NT'])]

       
        sc.tl.pca(adata_sub)

        
        nor_counts = adata_sub.X #the normalize count of the adata 
        if not issparse(nor_counts):
              nor_counts = np.asarray(nor_counts)
        
        
        loadings = adata_pb.varm['PCs'] # pseudobulk loading 
        pca_scores = np.dot(nor_counts, loadings)
        pca_scores = np.array(pca_scores)
        sc.pl.pca(adata_sub,color='sample_group')

        adata_sub.obsm['X_pca'] = pca_scores  #replace the raw pca scores
        sc.pp.neighbors(adata_sub, n_pcs=n_components)
        sc.tl.umap(adata_sub)
        sc.pl.umap(adata_sub, color=['sample_group'])
        with open('adata_pb.pkl', 'wb') as f:
            pickle.dump(adata_pb, f)
        with open('adata_sub.pkl', 'wb') as f:
            pickle.dump(adata_sub, f)

        return adata_pb,adata_sub


def run_gene_module(adata_sub):
    hs = hotspot.Hotspot(
            adata_sub,
            model='normal',
            latent_obsm_key="X_pca",
            umi_counts_obs_key="total_counts"
        )
    hs.create_knn_graph(weighted_graph=False, n_neighbors=30)
    hs_results = hs.compute_autocorrelations(jobs=96)

    with open('hs_results.pkl', 'wb') as f:
        pickle.dump(hs_results, f) 
    hs_results.to_csv('hs_results.csv')   

    hs_genes = hs_results.loc[hs_results.FDR < 0.05].index # Select genes

    local_correlations = hs.compute_local_correlations(hs_genes, jobs=96) 
    with open('local_correlations.pkl', 'wb') as f:
        pickle.dump(local_correlations, f)

    modules = hs.create_modules(
        min_gene_threshold=8, core_only=True, fdr_threshold=0.05
    ) 
    with open('modules.pkl', 'wb') as f:
        pickle.dump(modules, f)
    with open('hs.pkl', 'wb') as f:
        pickle.dump(hs, f)

# for KC example:
adata_pb,adata_sub = create_adata_sub('/home/yangyu/Project/skinAtlas/70samples/Result/1.celltype_pseudobulkDE/Female_KC_DE_adjusted_expression.txt','/home/yangyu/Project/skinAtlas/70samples/Analysis/results/23.KC_VIMKC/3.kc_label_final.h5ad','F')
with open('adata_sub.pkl', 'wb') as f:
    pickle.dump(adata_sub, f)
with open('adata_pb.pkl', 'wb') as f:
    pickle.dump(adata_pb, f)
run_gene_module(adata_sub)

hs = pickle.load(open('hs.pkl', 'rb'))
modules = pickle.load(open('modules.pkl', 'rb'))
local_correlations = pickle.load(open('local_correlations.pkl', 'rb'))
hs_results = pickle.load(open('hs_results.pkl', 'rb'))

module_scores = hs.calculate_module_scores()

module_scores.to_csv('module_scores.csv')

hs_genes = hs_results.loc[hs_results.FDR < 0.05].index # Select genes
hs_genes

sc.pl.umap(adata_sub, color=['kc_cell_type_scVI1','sample_group'], frameon=False,ncols=2,save='samplegroup_umap.pdf')

hs.plot_local_correlations(mod_cmap="tab20",vmin=-10,vmax=10)
plt.savefig('figures/local_correlations.pdf')

module_cols = []
for c in module_scores.columns:
    key = f"Module {c}"
    adata_sub.obs[key] = module_scores[c]
    module_cols.append(key)

with mplscience.style_context():
    sc.pl.umap(adata_sub, color=module_cols, frameon=False,cmap='Reds',save='module_umap.pdf')

pap_df = modules.to_frame()
pap_df
pap_df.to_csv( "Hotspot_gene_modules.csv")

# get all modules in a dataframe
all_mods = hs.results.join(hs.modules)
all_mods.dropna( inplace=True)
# drop non assigned genes set to -1.0
all_mods = all_mods[all_mods['Module']!=-1.0]
all_mods.shape

all_mods.Module = all_mods.Module.astype('int64')
all_mods.to_csv('all-cells_DEGs_hotspot-gene-modules.csv')

def member_test( A, B):
    b = set(B)
    return( [x in b for x in A])

# dataframe to hold cluster stage module expression
csme = pd.DataFrame([])
# loop through groups
for group in ['NF', 'NT', 'NE', 'NP']:
        #print(group)
    group_mk = adata_sub.obs['sample_group']==group
        # only want to look at group with over 10 cells in a cluster
    if( sum( group_mk)<10):
        continue
    adata_itr1 = adata_sub[group_mk]
    ind_nm = f"{group}"
        # loop through modules
    for mod_itr in np.sort(all_mods['Module'].unique()):
            #print(mod_itr)
        if( int(mod_itr)==-1):
            continue
        mod_gene_mk = all_mods['Module']==mod_itr
        mod_genes = all_mods.index.values[mod_gene_mk]
        adata_gene_mk = np.array( member_test(  adata_sub.var_names, mod_genes))
        gene_csr = adata_itr1[:,adata_gene_mk].X
        if( gene_csr.sum(0).sum()==0.0):
            print( f"nothing in {ind_nm}")
            csme.loc[ind_nm,f"Module-{str(int(mod_itr))}"] = 0.0
        else:
            gene_mean = gene_csr.mean(0)
            csme.loc[ind_nm,f"Module-{str(int(mod_itr))}"] = np.mean( gene_mean)

cmap = sns.color_palette( "ch:start=.2,rot=-.3", n_colors=1000)
sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white'})
g = sns.clustermap( csme.T, col_cluster=False, #row_cluster=True,
                   standard_scale=0, cmap=cmap, yticklabels=csme.columns.values,
                   linewidths=0.025, cbar_pos=(0.25, 0.92, 0.325, 0.05), 
                    cbar_kws={"orientation": "horizontal"})
g.ax_row_dendrogram.set_visible(False)
g.cax.set_visible(True)
plt.setp( g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0);
plt.savefig( "figures/hotspot_deg-samplegroup-vs-modules.pdf",dpi=1200, bbox_inches='tight')

