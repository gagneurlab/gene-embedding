# script to precompute stratified kfold

# packages
import pandas as pd
import polars as pl
from sklearn.model_selection import StratifiedKFold

# paths from snakemake
coding_genes_path = snakemake.input['coding_genes_path']
disease_gene_path = snakemake.input['disease_gene_path']
fold_splits_path = snakemake.output["fold_splits_path"]
label_table_path = snakemake.input['label_table_path']


# load tables
# coding_genes_df = pd.read_csv(coding_genes_path, sep='\t').dropna().rename(columns={'Gene stable ID':'gene_id'}).drop_duplicates().reset_index()
# check which genes to use
coding_genes_df = pl.read_csv(coding_genes_path, sep='\t')
coding_genes_df = coding_genes_df.filter(
    pl.col("Transcript type") == "protein_coding"
).filter(
    pl.col("Chromosome/scaffold name").is_in([str(i) for i in range(1,23)] + ["X", "Y", "MT"])
).select("Gene stable ID").unique().to_pandas().rename(columns={'Gene stable ID':'gene_id'})
disease_gene_table = pd.read_csv(disease_gene_path, sep='\t')

skf = StratifiedKFold(shuffle=True, random_state=1234)

split_table = []
for disease_name, disease_group in disease_gene_table.groupby('disease'):
    print(f'creating splits for {disease_name}')
    # get targets
    y = coding_genes_df['gene_id'].isin(disease_group['gene_id']).values
    
    for fold, (train_index, test_index) in enumerate(skf.split(y, y)):
        split_table.append(pd.DataFrame({
            'gene_id' : coding_genes_df.gene_id,
            'target': y,
            'train' : coding_genes_df.gene_id.isin(coding_genes_df.loc[train_index].gene_id),
            'test' : coding_genes_df.gene_id.isin(coding_genes_df.loc[test_index].gene_id),
            'fold' : fold + 1,
            'disease' : disease_name,
            'min_auprc': y.sum() / len(coding_genes_df)
        }))

        
split_table = pd.concat(split_table)

# add short names:
short_label_table = pd.read_csv(label_table_path, sep='\t')
split_table = split_table.merge(short_label_table, on='disease')


split_table.to_csv(fold_splits_path, sep='\t', index = False)
