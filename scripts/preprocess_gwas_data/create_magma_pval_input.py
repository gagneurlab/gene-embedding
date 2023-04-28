# import packages
import pandas as pd

print(snakemake)
# summary statistics
# select wanted data
finemap_table = pd.read_csv(snakemake.input.ukbb_subset, sep='\t', header=None)
rsid_mapping_table = finemap_table[finemap_table[11] == snakemake.params.trait][[3, 4]]
rsid_mapping_table.columns = ['variant', 'rsid']

nealelab_table = pd.read_csv(snakemake.input.raw_sum_stats, sep='\t')
nealelab_table['variant'] = nealelab_table['variant'].transform(lambda x: "chr" + x)

# merge rsid_mapping to nealelab in order to have the correct IDs with the Pvals
merged_table = rsid_mapping_table.merge(nealelab_table, how='left', on = 'variant')
#rows_to_drop = merged_table[not merged_table['rsis'].startswith('rs', end=3)].index
#merged_table.drop(rows_to_drop, inplace=True)
# put columns needed for magma to file
merged_table[merged_table['rsid'].str.contains('rs')].to_csv(snakemake.output.trait_sum_stats,
                                                             sep='\t', 
                                                             columns=['rsid', 'pval'], 
                                                             header=['SNP', 'P'],
                                                             index=False)

# meta data
ukbb_trait_table = pd.read_csv(snakemake.input.trait_annotation, sep='\t')
selected_trait = ukbb_trait_table[ukbb_trait_table["trait"] == snakemake.params.trait]
# create table with wanted format
meta_table = pd.DataFrame({
    "key":["trait", "number_of_individuals"],
    "value": [snakemake.params.trait, selected_trait["n"].values[0]]
})
#save to file
meta_table.to_csv(snakemake.output.meta_data_file, index=False, sep=':')
