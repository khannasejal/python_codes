#%%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt 
import numpy as np
from scipy.stats import zscore
import os
# %%
gene_exp_df = pd.read_csv('TRANSCRIPT_TAMR_FASR_MCF7_expression_data.tsv', sep='\t')
#%%
cols_to_transform = ['Rep1.TAMR.TPM','Rep2.TAMR.TPM','Rep3.TAMR.TPM','Rep1.FASR.TPM','Rep2.FASR.TPM','Rep3.FASR.TPM','Rep1.MCF7.TPM','Rep2.MCF7.TPM','Rep3.MCF7.TPM']
#do log transformation
gene_exp_df[cols_to_transform] = gene_exp_df[cols_to_transform].apply(lambda x: np.log2(x+1))
#do z normalisation across the rows
gene_exp_df[cols_to_transform] = gene_exp_df[cols_to_transform].apply(lambda x: zscore(x,axis=0),axis=1)

# %%
root_dir = "path_to_genesets"
genesets = ['epi_corr','mes','partial','emt','lum','basal','er']
#function to read the genesets which start with a prefix 

def get_genesets(geneset):
    for file_name in os.listdir(root_dir):
        if file_name.startswith(geneset):
            file_path = os.path.join(root_dir,file_name)
            break
    genes = pd.read_csv(file_path, sep='\t',header = None)
    return genes.values.flatten().tolist()
#%%
def get_heatmap_gene_exp_for_a_geneset(geneset):
    gene_list = get_genesets(geneset)
    gene_exp_df_filt = gene_exp_df[gene_exp_df['gene.symbol'].isin(gene_list)]
    gene_exp_df_filt.set_index('transcript.name', inplace=True)
    col = ['Rep1.MCF7.TPM','Rep2.MCF7.TPM','Rep3.MCF7.TPM','Rep1.FASR.TPM','Rep2.FASR.TPM','Rep3.FASR.TPM','Rep1.TAMR.TPM','Rep2.TAMR.TPM','Rep3.TAMR.TPM']
    gene_exp_df_filt_plot = gene_exp_df_filt[col]
    gene_exp_df_filt_plot = gene_exp_df_filt_plot[gene_exp_df_filt_plot[col].sum(axis=1) != 0]
    map = sns.clustermap(gene_exp_df_filt_plot, row_cluster=True, col_cluster=False,cmap='RdBu_r', figsize=(8,20))
    map.figure.suptitle(f'Clustermap for {geneset}')
    plt.show()

#%%
def get_boxplot_gene_exp_for_a_geneset(geneset):
    gene_list = get_genesets(geneset)
    gene_exp_df_filt = gene_exp_df[gene_exp_df['gene.symbol'].isin(gene_list)]
    gene_exp_df_filt.set_index('transcript.name', inplace=True)
    col = ['Rep1.MCF7.TPM', 'Rep2.MCF7.TPM', 'Rep3.MCF7.TPM', 
           'Rep1.FASR.TPM', 'Rep2.FASR.TPM', 'Rep3.FASR.TPM', 
           'Rep1.TAMR.TPM', 'Rep2.TAMR.TPM', 'Rep3.TAMR.TPM']
    gene_exp_df_filt_plot = gene_exp_df_filt[col]
    gene_exp_df_filt_plot = gene_exp_df_filt_plot[gene_exp_df_filt_plot[col].sum(axis=1) != 0]
    gene_exp_long = gene_exp_df_filt_plot.reset_index().melt(id_vars='transcript.name', 
                                                             value_vars=col, 
                                                             var_name='Replicate', 
                                                             value_name='TPM')
    palette = {'Rep1.MCF7.TPM': 'lightyellow', 'Rep2.MCF7.TPM': 'lightyellow', 'Rep3.MCF7.TPM': 'lightyellow',
        'Rep1.FASR.TPM': '#F7A541', 'Rep2.FASR.TPM': '#F7A541', 'Rep3.FASR.TPM': '#F7A541',
        'Rep1.TAMR.TPM': '#E9564A', 'Rep2.TAMR.TPM': '#E9564A', 'Rep3.TAMR.TPM': '#E9564A'}
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=gene_exp_long, x='Replicate', y='TPM', palette=palette)
    plt.title(f'Box Plot for {geneset} Gene Expression')
    plt.xticks(rotation=45)
    plt.show()
# %%
for i in genesets:
    get_boxplot_gene_exp_for_a_geneset(i)
