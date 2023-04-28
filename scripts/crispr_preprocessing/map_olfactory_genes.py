
import pandas as pd

olfactory_genes = pd.read_csv(
        snakemake.input['olfactory_genes'], 
        header=None, names = ['olfactory_genes'])
    

mapping_table = pd.read_csv('resources/hugo_symbol_to_ENSG_mapping.tsv', sep = '\t')
mapping_table = mapping_table[['Ensembl_gene_ID', 'any_symbol']].drop_duplicates().rename(columns = {"Ensembl_gene_ID": "gene_id"})
olfactory_genes = olfactory_genes.merge(mapping_table, left_on = 'olfactory_genes', right_on = 'any_symbol', how = 'left')
olfactory_genes


olfactory_genes = olfactory_genes.drop_duplicates('gene_id').drop(columns=['any_symbol'])

olfactory_genes.to_csv(snakemake.output[0], sep = '\t', index = False)
