

import numpy as np
#import snakemake
import pandas as pd
from sklearn.decomposition import PCA


def load_screens(gene_effect, sample_info, olfactory_genes):
    # Load screens
    screens = pd.read_csv(gene_effect, index_col=0).T
    screens.index = screens.index.str.split(' ').str.get(0)
    # Map Broad ID to CCLE name
    cell_lines = pd.read_csv(sample_info, 
                             index_col='Broad_ID',
                             usecols=['Broad_ID', 'CCLE_name'], squeeze=True)
    screens.columns = screens.columns.reindex(cell_lines)[0]
    # Bias-correct using "molecular function:olfactory receptor activity" genes
    olfactory_genes = pd.read_csv(
        olfactory_genes, header=None, squeeze=True)
    olfactory_data = screens.reindex(olfactory_genes).dropna()
    transformation = PCA(n_components=4)
    transformation.fit(olfactory_data)
    top_PC_effects = transformation.inverse_transform(
        transformation.transform(screens))
    screens -= top_PC_effects
    screens = screens.iloc[:, :-4]
    return screens


# Load batch-corrected screens
screens = load_screens(snakemake.input['gene_effect'],
                       snakemake.input['sample_info'],
                       snakemake.input['olfactory_genes'],
                      )

# Remove cell lines with any missing genes
# (not required for DepMap 18Q3, but is for more recent releases)
# You can use other strategies to remove NaNs instead, like imputing,
# removing genes with any missing cell lines

screens.dropna(axis=1, inplace=True)

# Warp screen data and intercept based on covariance of screens

screens_df = pd.DataFrame(screens, index = screens.index)

# # map to ENSG

mapping_table = pd.read_csv('resources/hugo_symbol_to_ENSG_mapping.tsv', sep = '\t')
mapping_table = mapping_table[['Ensembl_gene_ID', 'any_symbol']].drop_duplicates().rename(columns = {"Ensembl_gene_ID": "gene_id"})
screens_df = screens_df.merge(mapping_table, left_index=True, right_on = 'any_symbol', how = 'left')
screens_df


screens_df = screens_df.drop_duplicates('gene_id').drop(columns=['any_symbol']).set_index('gene_id')

screens_df.to_csv(snakemake.output[0], sep = '\t')
