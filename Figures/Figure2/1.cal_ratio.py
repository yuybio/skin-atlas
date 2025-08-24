# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# =============== read metadata  ===============
metadata = pd.read_csv('/home/yangyu/Project/skinAtlas/70samples/Fig/Figure1/seu.metadata.csv', index_col=0)
subtype_order = [
    'Basal.Ker.COL17A1',
    'Basal.Ker.ADGRL3',
    'Proliferating.Ker',
    'Spinous.Ker.ADGRL3',
    'Spinous.Ker',
    'Granular.Ker','Transitional.Ker','Hillock-Club.Epi',
    'Channel-ATPase_1', 'Channel-ATPase_2', 'Channel-ATPase_3', 'Channel-gap_1', 'Channel-gap_2', 'Sweat_gland_duct_cell',
    'IRS', 'ORS','HF-MSCs','HF-PCs', 'HF-MFCs',
    'Clear.cell.CA2', 'Clear.cell.WFDC2', 'Dark.cell.MUCL1','Sebocytes','MEC',
    'Progenitor.Mel', 'Mel_1', 'Mel_2','KC.Mel',
    'Merkel', 'SCs_1', 'SCs_2',
    'Papillary.Fb',
'Mesenchymal.Fb',
'Mesenchymal.Fb.FGFBP2',
'Mesenchymal.Fb.INHBA',
'Mesenchymal.Fb.COCH',
'Pro-inflammatory.Fb',
'arteriole.ECs', 
'capillary.ECs', 
'post-capillary-venule.ECs',
'venule.ECs',
'EMT.ECs',
'LEC1','LEC2',
 'Pc1',
 'Pc2',
 'Pc3',
 'vSMC1',
 'vSMC2',
    'M1', 'M2','Monocytes','cDC1', 'cDC2' ,'MigDC',
    'LC1', 'LC2','MigLC',   
    'Mast cells','helper T cells',
    'Regulatory T cells', 
    'cytotoxic T cells', 
    'NKT cells',
    'ILC3', 
    'Memory B cells',
    'Plasma cells'
]
celltype_order = ['KC','KC_Channel','HFC','SGC', 'Sebocyte','MEC','MEL', 'Schwann','FB', 'VEC','LEC','Pc-vSMC',
                       'Mac-DC','LC','Mast', 'Lymphocyte']


for col in ["cell_type", "celltype_Granular", "sample_id", "sample_group"]:
    if col in metadata.columns:
        metadata[col] = metadata[col].astype("category")

female_metadata = metadata.loc[metadata["sex"] == "F"].copy()
male_metadata   = metadata.loc[metadata["sex"] == "M"].copy()


if "celltype_order" in globals():
    for df in (female_metadata, male_metadata):
        present = [x for x in celltype_order if x in df["cell_type"].unique().tolist()]
        if present:
            df["cell_type"] = pd.Categorical(df["cell_type"], categories=present, ordered=True)

def build_global_expected(df, sample_col='sample_id', group_col='sample_group'):
    
    counts = (
        df.groupby(sample_col).size()
          .reset_index(name='sample_expected')
          .merge(df[[sample_col, group_col]].drop_duplicates(), on=sample_col, how='left')
    )
    total = len(df)
    counts['Total_expected'] = total
    counts['Percentage_expected'] = counts['sample_expected'] / counts['Total_expected'] * 100.0
    return counts

# =============== Female Broad ===============
final_df   = pd.DataFrame()
results_df = pd.DataFrame()
results    = []

cell_types_f = (female_metadata["cell_type"].cat.categories
                if hasattr(female_metadata["cell_type"], "cat") else female_metadata["cell_type"].unique())
for cell_type in cell_types_f:
    print(cell_type)
    cell_type_metadata = female_metadata[female_metadata['cell_type'] == cell_type]

    #
    counts_per_sample = (
        cell_type_metadata.groupby(['sample_id', 'cell_type'])
        .size().reset_index(name='sample_observed')
        .merge(cell_type_metadata[['sample_id', 'sample_group']].drop_duplicates(), on='sample_id')
    )
    total_observed = (
        cell_type_metadata.groupby('cell_type').size()
        .reset_index(name='Total_observed')
        .merge(cell_type_metadata[['sample_id','cell_type']].drop_duplicates(), on='cell_type')
    )
    actual_counts = pd.merge(counts_per_sample, total_observed, on=['cell_type','sample_id'])
    actual_counts['Percentage_observed'] = (actual_counts['sample_observed'] / actual_counts['Total_observed']) * 100

    # 
    theoretical_counts = build_global_expected(female_metadata)

    # 
    combined_df = pd.merge(actual_counts, theoretical_counts, on=['sample_id','sample_group'])
    combined_df['Percentage_Ratio'] = combined_df['Percentage_observed'] / combined_df['Percentage_expected']

    # 
    group_stats = combined_df.groupby(['sample_group', 'cell_type'])['Percentage_Ratio'].agg(mean_raw='mean', std_raw='std').reset_index()
    combined_df = combined_df.merge(group_stats, on=['sample_group', 'cell_type'])
    combined_df_filtered = combined_df[np.abs(combined_df['Percentage_Ratio'] - combined_df['mean_raw']) <= 3 * combined_df['std_raw']]

    # 
    final_group_stats = (
        combined_df_filtered.groupby(['sample_group', 'cell_type'])['Percentage_Ratio']
        .mean().reset_index().rename(columns={'Percentage_Ratio': 'Percentage_Ratio_Ave'})
    )
    combined_df_final = combined_df_filtered.merge(final_group_stats, on=['sample_group', 'cell_type'])
    combined_df_final['Log_Percentage_Ratio_Ave'] = np.log2(combined_df_final['Percentage_Ratio_Ave'].replace(0, np.nan))

    # 
    combined_df.to_csv(f'Female_{cell_type}.Percentage_Ratio_raw.csv',index=False)
    combined_df_final.to_csv(f'Female_{cell_type}.Percentage_Ratio_ave.csv',index=False)

    # 
    for ct in combined_df_final['cell_type'].unique():
        combined_df_subset = combined_df_final[combined_df_final['cell_type'] == ct]
        for target_group in combined_df_final['sample_group'].unique():
            target_data = combined_df_subset[combined_df_subset['sample_group'] == target_group]['Percentage_Ratio']
            rest_data   = combined_df_subset[combined_df_subset['sample_group'] != target_group]['Percentage_Ratio']
            if not target_data.empty and not rest_data.empty:
                _, p = mannwhitneyu(target_data, rest_data, alternative='two-sided')
                results.append({
                    'cell_type': cell_type,
                    'sample_group': target_group,
                    'Percentage_Ratio_Ave': combined_df_subset[combined_df_subset['sample_group'] == target_group]['Percentage_Ratio_Ave'].values[0],
                    'Log_Percentage_Ratio_Ave': combined_df_subset[combined_df_subset['sample_group'] == target_group]['Log_Percentage_Ratio_Ave'].values[0],
                    'p_value': p
                })

results_df = pd.DataFrame(results)
if len(results_df):
    results_df['corrected_p_value'] = multipletests(results_df['p_value'], alpha=0.05, method='fdr_bh')[1]
final_df = pd.concat([final_df, results_df], ignore_index=True)
final_df.to_csv('Female_ALL_cell.Percentage_Ratio_ave.csv',index=False)

#
pivot_table_df = final_df.pivot(index="cell_type", columns="sample_group", values="Log_Percentage_Ratio_Ave")
pivot_sig_df   = final_df.pivot(index="cell_type", columns="sample_group", values="p_value")  
pivot_table_df.loc['HFC','NP'] = np.nan
pivot_sig_df.loc['HFC','NP'] = np.nan
# 
if "celltype_order" in globals():
    current_order = [item for item in celltype_order if item in pivot_table_df.index]
    if len(current_order):
        pivot_table_df = pivot_table_df.loc[current_order]
        pivot_sig_df   = pivot_sig_df.loc[current_order]

#
want_cols = ['NF','NP','NE','NT']
present_cols = [c for c in want_cols if c in pivot_table_df.columns]
if len(present_cols):
    pivot_table_df = pivot_table_df.reindex(columns=present_cols)
    pivot_sig_df   = pivot_sig_df.reindex(columns=present_cols)

def significance_mark(p):
    return "*" if (pd.notna(p) and p <= 0.05) else ""

plt.figure(figsize=(3,7))
sns.heatmap(
    pivot_table_df,
    annot=np.vectorize(significance_mark)(pivot_sig_df), fmt='',
    annot_kws={"size": 10, "color": "black"},
    cmap='RdBu_r', center=0,
    linewidths=0.5, linecolor='white',
    vmin=-2, vmax=2,
    cbar_kws={'label': 'Log2(Ro/e)', 'shrink': 0.8}
)
plt.xticks(rotation=0)
plt.yticks(rotation=0)
plt.savefig('Female_all_Ratio.pdf',dpi=1200,bbox_inches='tight')  # 后面 granular 会覆盖
plt.show()

# =============== Male Broad ===============
final_df   = pd.DataFrame()
results_df = pd.DataFrame()
results    = []

cell_types_m = (male_metadata["cell_type"].cat.categories
                if hasattr(male_metadata["cell_type"], "cat") else male_metadata["cell_type"].unique())
for cell_type in cell_types_m:
    print(cell_type)
    cell_type_metadata = male_metadata[male_metadata['cell_type'] == cell_type]

    counts_per_sample = (
        cell_type_metadata.groupby(['sample_id', 'cell_type'])
        .size().reset_index(name='sample_observed')
        .merge(cell_type_metadata[['sample_id', 'sample_group']].drop_duplicates(), on='sample_id')
    )
    total_observed = (
        cell_type_metadata.groupby('cell_type').size()
        .reset_index(name='Total_observed')
        .merge(cell_type_metadata[['sample_id','cell_type']].drop_duplicates(), on='cell_type')
    )
    actual_counts = pd.merge(counts_per_sample, total_observed, on=['cell_type','sample_id'])
    actual_counts['Percentage_observed'] = (actual_counts['sample_observed'] / actual_counts['Total_observed']) * 100

    theoretical_counts = build_global_expected(male_metadata)

    combined_df = pd.merge(actual_counts, theoretical_counts, on=['sample_id','sample_group'])
    combined_df['Percentage_Ratio'] = combined_df['Percentage_observed'] / combined_df['Percentage_expected']

    group_stats = combined_df.groupby(['sample_group', 'cell_type'])['Percentage_Ratio'].agg(mean_raw='mean', std_raw='std').reset_index()
    combined_df = combined_df.merge(group_stats, on=['sample_group', 'cell_type'])
    combined_df_filtered = combined_df[np.abs(combined_df['Percentage_Ratio'] - combined_df['mean_raw']) <= 3 * combined_df['std_raw']]

    final_group_stats = (
        combined_df_filtered.groupby(['sample_group', 'cell_type'])['Percentage_Ratio']
        .mean().reset_index().rename(columns={'Percentage_Ratio': 'Percentage_Ratio_Ave'})
    )
    combined_df_final = combined_df_filtered.merge(final_group_stats, on=['sample_group', 'cell_type'])
    combined_df_final['Log_Percentage_Ratio_Ave'] = np.log2(combined_df_final['Percentage_Ratio_Ave'].replace(0, np.nan))

    combined_df.to_csv(f'Male_{cell_type}.Percentage_Ratio_raw.csv',index=False)
    combined_df_final.to_csv(f'Male_{cell_type}.Percentage_Ratio_ave.csv',index=False)

    for ct in combined_df_final['cell_type'].unique():
        combined_df_subset = combined_df_final[combined_df_final['cell_type'] == ct]
        for target_group in combined_df_final['sample_group'].unique():
            target_data = combined_df_subset[combined_df_subset['sample_group'] == target_group]['Percentage_Ratio']
            rest_data   = combined_df_subset[combined_df_subset['sample_group'] != target_group]['Percentage_Ratio']
            if not target_data.empty and not rest_data.empty:
                _, p = mannwhitneyu(target_data, rest_data, alternative='two-sided')
                results.append({
                    'cell_type': cell_type,
                    'sample_group': target_group,
                    'Percentage_Ratio_Ave': combined_df_subset[combined_df_subset['sample_group'] == target_group]['Percentage_Ratio_Ave'].values[0],
                    'Log_Percentage_Ratio_Ave': combined_df_subset[combined_df_subset['sample_group'] == target_group]['Log_Percentage_Ratio_Ave'].values[0],
                    'p_value': p
                })

results_df = pd.DataFrame(results)
if len(results_df):
    results_df['corrected_p_value'] = multipletests(results_df['p_value'], alpha=0.05, method='fdr_bh')[1]
final_df = pd.concat([final_df, results_df], ignore_index=True)

pivot_table_df_m = final_df.pivot(index="cell_type", columns="sample_group", values="Log_Percentage_Ratio_Ave")
pivot_sig_df_m   = final_df.pivot(index="cell_type", columns="sample_group", values="p_value") 

if "celltype_order" in globals():
    current_order_m = [item for item in celltype_order if item in pivot_table_df_m.index]
    if len(current_order_m):
        pivot_table_df_m = pivot_table_df_m.loc[current_order_m]
        pivot_sig_df_m   = pivot_sig_df_m.loc[current_order_m]

want_cols_m = ['PS','NE','NT']
present_cols_m = [c for c in want_cols_m if c in pivot_table_df_m.columns]
if len(present_cols_m):
    pivot_table_df_m = pivot_table_df_m.reindex(columns=present_cols_m)
    pivot_sig_df_m   = pivot_sig_df_m.reindex(columns=present_cols_m)

def _star(p):
    return "*" if (pd.notna(p) and p <= 0.05) else ""

plt.figure(figsize=(2.3,7))
sns.heatmap(
    pivot_table_df_m,
    annot=np.vectorize(_star)(pivot_sig_df_m), fmt='',
    annot_kws={"size": 10, "color": "black"},
    cmap='RdBu_r', center=0,
    linewidths=0.5, linecolor='white',
    vmin=-2, vmax=2,
    cbar_kws={'label': 'Log2(Ro/e)', 'shrink': 0.8}
)
plt.xticks(rotation=0)
plt.yticks(rotation=0)
plt.savefig('Male_all_Ratio.pdf', dpi=1200, bbox_inches='tight')  
plt.show()

# =============== Female subtype ===============
final_df   = pd.DataFrame()
results_df = pd.DataFrame()
results    = []

pivot_table = {}
pivot_sig   = {}

cell_types_f = (female_metadata["cell_type"].cat.categories
                if hasattr(female_metadata["cell_type"], "cat") else female_metadata["cell_type"].unique())
for cell_type in cell_types_f:
    print(cell_type)

    if cell_type not in ['MEC','Sebocyte','Mast']:
        #
        cell_type_metadata = female_metadata[female_metadata['cell_type'] == cell_type]
        counts_per_sample = (
            cell_type_metadata.groupby(['sample_id', 'celltype_Granular'])
            .size().reset_index(name='sample_observed')
            .merge(cell_type_metadata[['sample_id', 'sample_group']].drop_duplicates(), on='sample_id')
        )
        total_observed = (
            cell_type_metadata.groupby('celltype_Granular').size()
            .reset_index(name='Total_observed')
            .merge(cell_type_metadata[['sample_id','celltype_Granular','cell_type']].drop_duplicates(), on='celltype_Granular')
        )
        actual_counts = pd.merge(counts_per_sample, total_observed, on=['celltype_Granular','sample_id'])
        actual_counts['Percentage_observed'] = (actual_counts['sample_observed'] / actual_counts['Total_observed']) * 100

        # 
        conut_per_sample_theoretical_counts = (
            cell_type_metadata.groupby(['sample_id', 'cell_type'])
            .size().reset_index(name='sample_expected')
        )
        conut_per_sample_theoretical_counts = conut_per_sample_theoretical_counts[conut_per_sample_theoretical_counts.cell_type==cell_type]
        conut_per_sample_theoretical_counts = conut_per_sample_theoretical_counts.merge(
            cell_type_metadata[['sample_id','sample_group','cell_type']].drop_duplicates(),
            on=['sample_id','cell_type']
        )
        total_expected = (
            cell_type_metadata.groupby('cell_type').size()
            .reset_index(name='Total_expected')
        )
        total_expected = total_expected[total_expected.cell_type==cell_type]
        total_expected = total_expected.merge(
            cell_type_metadata[['sample_group','cell_type']].drop_duplicates(),
            on=['cell_type']
        )
        theoretical_counts = pd.merge(conut_per_sample_theoretical_counts, total_expected, on=['sample_group','cell_type'])
        theoretical_counts['Percentage_expected'] = (theoretical_counts['sample_expected'] / theoretical_counts['Total_expected']) * 100

        combined_df = pd.merge(actual_counts, theoretical_counts, on=['sample_id','sample_group','cell_type'])
        combined_df['Percentage_Ratio'] = combined_df['Percentage_observed'] / combined_df['Percentage_expected']
        #combined_df = combined_df[combined_df['sample_observed'] >= 5]
    else:
        # 
        cell_type_metadata = female_metadata[female_metadata['cell_type'] == cell_type]
        counts_per_sample = (
            cell_type_metadata.groupby(['sample_id', 'celltype_Granular'])
            .size().reset_index(name='sample_observed')
            .merge(cell_type_metadata[['sample_id', 'sample_group']].drop_duplicates(), on='sample_id')
        )
        total_observed = (
            cell_type_metadata.groupby('celltype_Granular').size()
            .reset_index(name='Total_observed')
            .merge(cell_type_metadata[['sample_id','celltype_Granular','cell_type']].drop_duplicates(), on='celltype_Granular')
        )
        actual_counts = pd.merge(counts_per_sample, total_observed, on=['celltype_Granular','sample_id'])
        actual_counts['Percentage_observed'] = (actual_counts['sample_observed'] / actual_counts['Total_observed']) * 100

        # 
        theoretical_counts = build_global_expected(female_metadata)

        combined_df = pd.merge(actual_counts, theoretical_counts, on=['sample_id','sample_group'])
        combined_df['Percentage_Ratio'] = combined_df['Percentage_observed'] / combined_df['Percentage_expected']
        #combined_df = combined_df[combined_df['sample_observed'] >= 5]

    # 
    group_stats = combined_df.groupby(['sample_group', 'celltype_Granular'])['Percentage_Ratio'].agg(mean_raw='mean', std_raw='std').reset_index()
    combined_df = combined_df.merge(group_stats, on=['sample_group', 'celltype_Granular'])
    combined_df_filtered = combined_df[np.abs(combined_df['Percentage_Ratio'] - combined_df['mean_raw']) <= 3 * combined_df['std_raw']]

    final_group_stats = (
        combined_df_filtered.groupby(['sample_group', 'celltype_Granular'])['Percentage_Ratio']
        .mean().reset_index().rename(columns={'Percentage_Ratio': 'Percentage_Ratio_Ave'})
    )
    combined_df_final = combined_df_filtered.merge(final_group_stats, on=['sample_group', 'celltype_Granular'])
    combined_df_final['Log_Percentage_Ratio_Ave'] = np.log2(combined_df_final['Percentage_Ratio_Ave'].replace(0, np.nan))

    combined_df.to_csv(f'Female_{cell_type}.Percentage_Ratio_raw.csv',index=False)
    combined_df_final.to_csv(f'Female_{cell_type}.Percentage_Ratio_ave.csv',index=False)

    # 
    for celltype_granular in combined_df_final['celltype_Granular'].unique():
        combined_df_subset = combined_df_final[combined_df_final['celltype_Granular'] == celltype_granular]
        for target_group in combined_df_final['sample_group'].unique():
            target_data = combined_df_subset[combined_df_subset['sample_group'] == target_group]['Percentage_Ratio']
            rest_data   = combined_df_subset[combined_df_subset['sample_group'] != target_group]['Percentage_Ratio']
            if not target_data.empty and not rest_data.empty:
                _, p = mannwhitneyu(target_data, rest_data, alternative='two-sided')
                results.append({
                    'cell_type': cell_type,
                    'celltype_Granular': celltype_granular,
                    'sample_group': target_group,
                    'Percentage_Ratio_Ave': combined_df_subset[combined_df_subset['sample_group'] == target_group]['Percentage_Ratio_Ave'].values[0],
                    'Log_Percentage_Ratio_Ave': combined_df_subset[combined_df_subset['sample_group'] == target_group]['Log_Percentage_Ratio_Ave'].values[0],
                    'p_value': p
                })

results_df = pd.DataFrame(results)
if len(results_df):
    results_df['corrected_p_value'] = multipletests(results_df['p_value'], alpha=0.05, method='fdr_bh')[1]
final_df = pd.concat([final_df, results_df], ignore_index=True)

# 
pivot_table = {}
pivot_sig   = {}

default_cols = ['NF','NT','NE','NP']
short_cols   = ['NF','NT','NE']  # Sebocyte, HFC

for cell_type in final_df['cell_type'].dropna().unique():
    df_sub = final_df[final_df['cell_type'] == cell_type]
    pvt = df_sub.pivot(index="celltype_Granular", columns="sample_group", values="Log_Percentage_Ratio_Ave")
    sig = df_sub.pivot(index="celltype_Granular", columns="sample_group", values="p_value")  # 原始 p 值

    if pvt.empty:
        continue

    want_cols = short_cols if cell_type in ['Sebocyte','HFC'] else default_cols
    present_cols = [c for c in want_cols if c in pvt.columns]
    if len(present_cols) == 0:
        continue

    pvt = pvt.reindex(columns=present_cols)
    sig = sig.reindex(columns=present_cols)

    if 'subtype_order' in globals():
        row_present = [r for r in subtype_order if r in pvt.index]
        if len(row_present):
            pvt = pvt.loc[row_present]
            sig = sig.loc[row_present]

    pivot_table[cell_type] = pvt
    pivot_sig[cell_type]   = sig

# 
if 'celltype_order' in globals():
    keys_in_order = [k for k in celltype_order if k in pivot_table]
else:
    keys_in_order = list(pivot_table.keys())

# 
cell_height = 0.22
cell_width  = 0.4
total_width = 0.0
total_height = 0.0
for ct in keys_in_order:
    pvt = pivot_table[ct]
    if pvt.shape[0] == 0 or pvt.shape[1] == 0:
        continue
    nrows, ncols = pvt.shape
    total_height += nrows * cell_height
    total_width  = max(total_width, ncols * cell_width)

extra_space = 0.05
total_height += extra_space * len(keys_in_order)

def star(p):
    return "*" if (pd.notna(p) and p <= 0.05) else ""

fig = plt.figure(figsize=(max(total_width, 3), max(total_height, 3)))
current_bottom = 0.0
for i, ct in enumerate(keys_in_order):
    pvt = pivot_table[ct]
    sig = pivot_sig[ct]
    if pvt.shape[0] == 0 or pvt.shape[1] == 0:
        continue
    sig_anno = np.vectorize(star)(sig)

    nrows, ncols = pvt.shape
    height_ratio = (nrows * cell_height) / total_height if total_height > 0 else 1.0
    width_ratio  = (ncols * cell_width)  / total_width  if total_width  > 0 else 1.0

    ax = fig.add_axes([0, current_bottom, width_ratio, height_ratio])
    sns.heatmap(
        pvt, annot=sig_anno, fmt='',
        annot_kws={"size": 10, "color": "black"},
        cmap='RdBu_r', center=0,
        linewidths=0.5, linecolor='white',
        vmin=-2, vmax=2,
        cbar_kws={'label': 'Log2(Ro/e)', 'shrink': 0.5},
        ax=ax
    )
    if i == 0:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    else:
        ax.set_xticks([])
    ax.set_ylabel(''); ax.set_xlabel('')
    ax.axhline(y=pvt.shape[0], color='white', lw=0.1)

    current_bottom += height_ratio + (extra_space / total_height if total_height > 0 else 0.0)

plt.savefig('Female_all_Ratio_subtype.pdf', dpi=1200, bbox_inches='tight')  
plt.show()

# =============== Male subtype ===============
final_df   = pd.DataFrame()
results_df = pd.DataFrame()
results    = []

cell_types_m = (male_metadata["cell_type"].cat.categories
                if hasattr(male_metadata["cell_type"], "cat") else male_metadata["cell_type"].unique())
for cell_type in cell_types_m:
    print(cell_type)

    if cell_type not in ['MEC','Sebocyte','Mast']:
        cell_type_metadata = male_metadata[male_metadata['cell_type'] == cell_type]
        counts_per_sample = (
            cell_type_metadata.groupby(['sample_id', 'celltype_Granular'])
            .size().reset_index(name='sample_observed')
            .merge(cell_type_metadata[['sample_id', 'sample_group']].drop_duplicates(), on='sample_id')
        )
        total_observed = (
            cell_type_metadata.groupby('celltype_Granular').size()
            .reset_index(name='Total_observed')
            .merge(cell_type_metadata[['sample_id','celltype_Granular','cell_type']].drop_duplicates(), on='celltype_Granular')
        )
        actual_counts = pd.merge(counts_per_sample, total_observed, on=['celltype_Granular','sample_id'])
        actual_counts['Percentage_observed'] = (actual_counts['sample_observed'] / actual_counts['Total_observed']) * 100

        # 
        conut_per_sample_theoretical_counts = (
            cell_type_metadata.groupby(['sample_id', 'cell_type'])
            .size().reset_index(name='sample_expected')
        )
        conut_per_sample_theoretical_counts = conut_per_sample_theoretical_counts[conut_per_sample_theoretical_counts.cell_type==cell_type]
        conut_per_sample_theoretical_counts = conut_per_sample_theoretical_counts.merge(
            cell_type_metadata[['sample_id','sample_group','cell_type']].drop_duplicates(), on=['sample_id','cell_type']
        )
        total_expected = (
            cell_type_metadata.groupby('cell_type').size()
            .reset_index(name='Total_expected')
        )
        total_expected = total_expected[total_expected.cell_type==cell_type]
        total_expected = total_expected.merge(
            cell_type_metadata[['sample_group','cell_type']].drop_duplicates(), on=['cell_type']
        )
        theoretical_counts = pd.merge(conut_per_sample_theoretical_counts, total_expected, on=['sample_group','cell_type'])
        theoretical_counts['Percentage_expected'] = (theoretical_counts['sample_expected'] / theoretical_counts['Total_expected']) * 100

        combined_df = pd.merge(actual_counts, theoretical_counts, on=['sample_id','sample_group','cell_type'])
        combined_df['Percentage_Ratio'] = combined_df['Percentage_observed'] / combined_df['Percentage_expected']
        #combined_df = combined_df[combined_df['sample_observed'] >= 5]
    else:
        cell_type_metadata = male_metadata[male_metadata['cell_type'] == cell_type]
        counts_per_sample = (
            cell_type_metadata.groupby(['sample_id', 'celltype_Granular'])
            .size().reset_index(name='sample_observed')
            .merge(cell_type_metadata[['sample_id', 'sample_group']].drop_duplicates(), on='sample_id')
        )
        total_observed = (
            cell_type_metadata.groupby('celltype_Granular').size()
            .reset_index(name='Total_observed')
            .merge(cell_type_metadata[['sample_id','celltype_Granular','cell_type']].drop_duplicates(), on='celltype_Granular')
        )
        actual_counts = pd.merge(counts_per_sample, total_observed, on=['celltype_Granular','sample_id'])
        actual_counts['Percentage_observed'] = (actual_counts['sample_observed'] / actual_counts['Total_observed']) * 100

        # 
        theoretical_counts = build_global_expected(male_metadata)

        combined_df = pd.merge(actual_counts, theoretical_counts, on=['sample_id','sample_group'])
        combined_df['Percentage_Ratio'] = combined_df['Percentage_observed'] / combined_df['Percentage_expected']
        #combined_df = combined_df[combined_df['sample_observed'] >= 5]

    # 
    group_stats = combined_df.groupby(['sample_group', 'celltype_Granular'])['Percentage_Ratio'].agg(mean_raw='mean', std_raw='std').reset_index()
    combined_df = combined_df.merge(group_stats, on=['sample_group', 'celltype_Granular'])
    combined_df_filtered = combined_df[np.abs(combined_df['Percentage_Ratio'] - combined_df['mean_raw']) <= 3 * combined_df['std_raw']]

    final_group_stats = (
        combined_df_filtered.groupby(['sample_group', 'celltype_Granular'])['Percentage_Ratio']
        .mean().reset_index().rename(columns={'Percentage_Ratio': 'Percentage_Ratio_Ave'})
    )
    combined_df_final = combined_df_filtered.merge(final_group_stats, on=['sample_group', 'celltype_Granular'])
    combined_df_final['Log_Percentage_Ratio_Ave'] = np.log2(combined_df_final['Percentage_Ratio_Ave'].replace(0, np.nan))

    combined_df.to_csv(f'Male_{cell_type}.Percentage_Ratio_raw.csv',index=False)
    combined_df_final.to_csv(f'Male_{cell_type}.Percentage_Ratio_ave.csv',index=False)

    for celltype_granular in combined_df_final['celltype_Granular'].unique():
        combined_df_subset = combined_df_final[combined_df_final['celltype_Granular'] == celltype_granular]
        for target_group in combined_df_final['sample_group'].unique():
            target_data = combined_df_subset[combined_df_subset['sample_group'] == target_group]['Percentage_Ratio']
            rest_data   = combined_df_subset[combined_df_subset['sample_group'] != target_group]['Percentage_Ratio']
            if not target_data.empty and not rest_data.empty:
                _, p = mannwhitneyu(target_data, rest_data, alternative='two-sided')
                results.append({
                    'cell_type': cell_type,
                    'celltype_Granular': celltype_granular,
                    'sample_group': target_group,
                    'Percentage_Ratio_Ave': combined_df_subset[combined_df_subset['sample_group'] == target_group]['Percentage_Ratio_Ave'].values[0],
                    'Log_Percentage_Ratio_Ave': combined_df_subset[combined_df_subset['sample_group'] == target_group]['Log_Percentage_Ratio_Ave'].values[0],
                    'p_value': p
                })

results_df = pd.DataFrame(results)
if len(results_df):
    results_df['corrected_p_value'] = multipletests(results_df['p_value'], alpha=0.05, method='fdr_bh')[1]
final_df = pd.concat([final_df, results_df], ignore_index=True)

pivot_table = {}
pivot_sig   = {}

target_cols = ['PS','NE','NT']  

for cell_type in final_df['cell_type'].dropna().unique():
    df = final_df[final_df['cell_type'] == cell_type]

    pvt = df.pivot(index="celltype_Granular", columns="sample_group",
                   values="Log_Percentage_Ratio_Ave")
    sig = df.pivot(index="celltype_Granular", columns="sample_group",
                   values="p_value")

    if pvt.empty:
        continue

    present = [c for c in target_cols if c in pvt.columns]
    if not present:
        continue
    pvt = pvt.reindex(columns=present)
    sig = sig.reindex(columns=present)

    if 'subtype_order' in globals():
        row_present = [r for r in subtype_order if r in pvt.index]
        if row_present:
            pvt = pvt.loc[row_present]
            sig = sig.loc[row_present]

    pivot_table[cell_type] = pvt
    pivot_sig[cell_type]   = sig



if 'subtype_order' in globals():
    keys_in_order = [k for k in subtype_order if k in pivot_table]
else:
    keys_in_order = list(pivot_table.keys())

cell_height = 0.22
cell_width  = 0.4
total_width = 0.0
total_height = 0.0
for ct in keys_in_order:
    p = pivot_table[ct]
    if p.shape[0] == 0 or p.shape[1] == 0:
        continue
    nrows, ncols = p.shape
    total_height += nrows * cell_height
    total_width  = max(total_width, ncols * cell_width)

extra_space = 0.05
total_height += extra_space * max(1, len(keys_in_order))

def star(p):
    return "*" if (pd.notna(p) and p <= 0.05) else ""

fig = plt.figure(figsize=(max(total_width, 3), max(total_height, 3)))
current_bottom = 0.0

for i, cell_type in enumerate(keys_in_order):
    pvt = pivot_table[cell_type]
    sig = pivot_sig[cell_type]
    if pvt.shape[0] == 0 or pvt.shape[1] == 0:
        continue

    sig_anno = np.vectorize(star)(sig)
    nrows, ncols = pvt.shape
    height_ratio = (nrows * cell_height) / total_height if total_height > 0 else 1.0
    width_ratio  = (ncols * cell_width)  / total_width  if total_width  > 0 else 1.0

    ax = fig.add_axes([0, current_bottom, width_ratio, height_ratio])
    sns.heatmap(
        pvt, annot=sig_anno, fmt='',
        annot_kws={"size": 10, "color": "black"},
        cmap='RdBu_r', center=0,
        linewidths=0.5, linecolor='white',
        vmin=-2, vmax=2,
        cbar_kws={'label': 'Log2(Ro/e)', 'shrink': 0.5},
        ax=ax
    )
    if i == 0:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    else:
        ax.set_xticks([])

    ax.set_ylabel(''); ax.set_xlabel('')
    ax.axhline(y=pvt.shape[0], color='white', lw=0.1)

    current_bottom += height_ratio + (extra_space / total_height if total_height > 0 else 0.0)

plt.savefig('Male_all_Ratio_subtype.pdf', dpi=1200, bbox_inches='tight')
plt.show()