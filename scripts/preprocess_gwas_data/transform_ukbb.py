# import packages
import pandas as pd
from scipy.stats import chi2

print(snakemake)
# summary statistics
# select wanted data
ukbb_table = pd.read_csv(snakemake.input.ukbb_subset, sep='\t', header=None)
subset_table = ukbb_table[ukbb_table[11] == snakemake.params.trait]
sum_stats_table = subset_table.iloc[:,[4, 16]]
# rename columns to wanted format
sum_stats_table.columns = ["SNP", "P"]
sum_stats_table["P"] = 1 - chi2.cdf(sum_stats_table["P"], 1)
# save summary statistics
sum_stats_table.to_csv(snakemake.output.trait_sum_stats, index=False, sep='\t', float_format='%.6f')

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
