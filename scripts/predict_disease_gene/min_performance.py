import pandas as pd

disease_gene_path = snakemake.input['disease_gene_path']
coding_genes_path = snakemake.input['coding_genes_path']
label_table_path = snakemake.input['label_table_path']

min_performance_path = snakemake.output['min_performance_path']

disease_gene_table = pd.read_csv(disease_gene_path, sep='\t')

coding_gene_table = pd.read_csv(coding_genes_path, sep='\t').dropna().drop_duplicates()

short_label_table = pd.read_csv(label_table_path, sep='\t')

min_perf_df = pd.DataFrame(disease_gene_table.disease.value_counts()).reset_index()
min_perf_df.columns = ['disease', 'count']
# compute the minimaly attainable area under the precision recall curve
min_perf_df['min_auprc'] = min_perf_df['count'] / len(coding_gene_table)

min_perf_df = min_perf_df.merge(short_label_table, on='disease')

# create label for plotting
min_perf_df['label'] = min_perf_df[['short_name', 'count']].apply(lambda x: f'{x[0]} ({x[1]})', axis=1)

min_perf_df.to_csv(min_performance_path, sep='\t', index=False)
