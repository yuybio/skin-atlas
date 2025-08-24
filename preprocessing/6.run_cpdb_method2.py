# -*- coding: utf-8 -*-
import scanpy as sc
import anndata as ad
import pandas as pd
import ktplotspy as kpy
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import pickle

cpdb_file_path = '/data01/home/yangyu/Pipeline/cellphonedb/v5.0.0/cellphonedb.zip'
meta_file_path = 'data/metadata.tsv'
counts_file_path = '/data01/home/yangyu/Project/skinAtlas/70samples/Analysis/h5adFiles.bak/4.atlas_Granular.h5ad'
microenvs_file_path = 'data/microenvironment.txt'
out_path = 'results/method2_withScore'

os.chdir('/data01/home/yangyu/Project/skinAtlas/70samples/Analysis/results/4.scVIintegrated/3.cpdb/celltype')

#adata = sc.read_h5ad(counts_file_path)
#adata.obs['group_type']=adata.obs['sample_group'].astype('str')+'_'+adata.obs['cell_type'].astype('str')
#adata.obs['group_type']=adata.obs['group_type'].astype('category')
#adata.write('/data01/home/yangyu/Project/skinAtlas/70samples/Analysis/h5adFiles.bak/4.atlas_Granular.h5ad')


from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

# cell type
cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path = cpdb_file_path,                 # mandatory: CellphoneDB database zip file.
    meta_file_path = meta_file_path,                 # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = counts_file_path,             # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
    counts_data = 'hgnc_symbol',                     # defines the gene annotation in counts matrix.
    #active_tfs_file_path = active_tf_path,           # optional: defines cell types and their active TFs.
    microenvs_file_path = microenvs_file_path,       # optional (default: None): defines cells per microenvironment.
    score_interactions = True,                       # optional: whether to score interactions or not. 
    iterations = 1000,                               # denotes the number of shufflings performed in the analysis.
    threshold = 0.1,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
    threads = 50,                                     # number of threads to use in the analysis.
    debug_seed = 222,                                 # debug randome seed. To disable >=0.
    result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
    pvalue = 0.05,                                   # P-value threshold to employ for significance.
    subsampling = False,                             # To enable subsampling the data (geometri sketching).
    subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
    subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
    subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
    separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
    output_path = out_path,                          # Path to save results.
    output_suffix = 'all'                             # Replaces the timestamp in the output files by a user defined string in the  (default: None).
    )
with open('cpdb_results_method2_celltype.pkl', 'wb') as f:
        pickle.dump(cpdb_results, f)


# sub type
os.chdir('/data01/home/yangyu/Project/skinAtlas/70samples/Analysis/results/4.scVIintegrated/3.cpdb/sub_celltype')
#adata = sc.read_h5ad(counts_file_path)
#adata.obs['group_type']=adata.obs['sample_group'].astype('str')+'_'+adata.obs['cell_type'].astype('str')
#adata.obs['group_type']=adata.obs['group_type'].astype('category')
#adata.write('/data01/home/yangyu/Project/skinAtlas/70samples/Analysis/h5adFiles.bak/4.atlas_Granular.h5ad')
cpdb_file_path = '/data01/home/yangyu/Pipeline/cellphonedb/v5.0.0/cellphonedb.zip'
meta_file_path = 'data/metadata.tsv'
counts_file_path = '/data01/home/yangyu/Project/skinAtlas/70samples/Analysis/h5adFiles.bak/4.atlas_Granular.h5ad'
microenvs_file_path = 'data/microenvironment.txt'
out_path = 'results/method2_withScore'

cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path = cpdb_file_path,                 # mandatory: CellphoneDB database zip file.
    meta_file_path = meta_file_path,                 # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = counts_file_path,             # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
    counts_data = 'hgnc_symbol',                     # defines the gene annotation in counts matrix.
    #active_tfs_file_path = active_tf_path,           # optional: defines cell types and their active TFs.
    microenvs_file_path = microenvs_file_path,       # optional (default: None): defines cells per microenvironment.
    score_interactions = True,                       # optional: whether to score interactions or not. 
    iterations = 1000,                               # denotes the number of shufflings performed in the analysis.
    threshold = 0.1,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
    threads = 50,                                     # number of threads to use in the analysis.
    debug_seed = 222,                                 # debug randome seed. To disable >=0.
    result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
    pvalue = 0.05,                                   # P-value threshold to employ for significance.
    subsampling = False,                             # To enable subsampling the data (geometri sketching).
    subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
    subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
    subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
    separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
    output_path = out_path,                          # Path to save results.
    output_suffix = 'all'                             # Replaces the timestamp in the output files by a user defined string in the  (default: None).
    )
with open('cpdb_results_method2_subcelltype.pkl', 'wb') as f:
        pickle.dump(cpdb_results, f)







