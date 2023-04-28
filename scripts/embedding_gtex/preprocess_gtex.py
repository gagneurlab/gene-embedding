# import packages
import pandas as pd
import polars as pl
import numpy as np
from sklearn.preprocessing import StandardScaler

# file paths
gtex_tpm_path = snakemake.input['gtex_data_path']
selected_genes_path = snakemake.input['which_genes']
gtex_subset_path = snakemake.output['gtex_subset_path']

gtex_tpm_table = pd.read_csv(gtex_tpm_path,
                             sep='\t', 
                             skiprows=2)
# reshape and prepare
del gtex_tpm_table["Description"]

# remove transcript notation as it is not used by other data
gtex_tpm_table["Name"] = gtex_tpm_table["Name"].apply(lambda x: x.split('.')[0])
gtex_tpm_table = gtex_tpm_table.drop_duplicates(subset=['Name'])
gtex_tpm_table.reset_index(inplace=True, drop=True)

# check which genes to use
selected_gene_set = pl.read_csv(selected_genes_path, sep='\t')
selected_gene_set = selected_gene_set.filter(
    pl.col("Transcript type") == "protein_coding"
).filter(
    pl.col("Chromosome/scaffold name").is_in([str(i) for i in range(1,23)] + ["X", "Y", "MT"])
).select("Gene stable ID").unique().to_pandas()

selected_gene_set = selected_gene_set.dropna().rename(columns={'Gene stable ID':'gene_id'}).drop_duplicates()
gtex_tpm_table = selected_gene_set[['gene_id']].merge(gtex_tpm_table, left_on='gene_id', right_on='Name',how='inner')
del gtex_tpm_table['Name']
gtex_tpm_table.dropna(inplace=True)
gtex_tpm_table = gtex_tpm_table.drop_duplicates(subset='gene_id')
gtex_tpm_table = gtex_tpm_table.set_index('gene_id')


gtex_tpm_table = np.log10((gtex_tpm_table + 1)) # log10 transform to get closer to normal distribution

# save to file
gtex_tpm_table.to_csv(gtex_subset_path, sep='\t')

