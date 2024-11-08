#%%
import pandas as pd
import gseapy as gp
import numpy as np

# %%
#read the gene expression (RPKM values) file from CCLE
ccle_breast_cancer_rnaseq = pd.read_csv("CCLE_RNAseq_genes_rpkm_20180929.gct", sep='\t', skiprows=2)
first_two_columns = ccle_breast_cancer_rnaseq.iloc[:, :2]
#Filter only breast cancer cell lines 
filtered_columns = [col for col in ccle_breast_cancer_rnaseq.columns[2:] if col.split('_')[1] == 'BREAST']
filtered_df = ccle_breast_cancer_rnaseq[filtered_columns]
result_df = pd.concat([first_two_columns, filtered_df], axis=1)
result_df = result_df.drop('Name', axis=1)
data_frame_expression_values = pd.DataFrame(result_df)  #rows are genes and columns are breast cancer cell lines

# %%
#read the geneset to be used for ssgsea
luminal = pd.read_excel('luminal_basal_geneset.xlsx', header =0)["Luminal"]
luminal = luminal.dropna()
luminal_geneset = {'luminal_set': luminal}

#function to do ssgsea (calculates enrichment score for a given geneset in a sample)
def do_ssgsea_with_geneset(geneset):
    ss = gp.ssgsea(data =data_frame_expression_values, gene_sets = geneset, outdir=None, sample_norm_method='rank', no_plot=True, min_size=10, max_size =2000)
    nes = ss.res2d.pivot(index='Term', columns='Name', values='NES')
    return nes

#%%
ss_luminal = do_ssgsea_with_geneset(luminal_geneset)