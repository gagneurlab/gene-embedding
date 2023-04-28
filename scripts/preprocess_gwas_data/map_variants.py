import pandas as pd
import numpy as np 


# load raw files
# variant annotation
variant_annotation_table = pd.read_csv(snakemake.input.snp_annotation_file, sep='\t',
                                       usecols=['chrom', 'pos','ref', 'alt', 'rsid'],
                                       dtype= str)
variant_annotation_table = variant_annotation_table.astype(str)

# meta information
nsamples_table = pd.read_csv(snakemake.input.sum_stats_sources, sep='\t', usecols=['TRAIT', 'NSAMPLES'])

# summary statistics
gwas_snp_table = None
# decide if we have to create the variant column
# for that read first line of table and look for pval_EUR
columns = pd.read_csv(snakemake.input.raw_sum_stats, sep='\t', nrows=1).columns
if 'pval_EUR' in columns:
    # we are dealing with LDLC or HDL
    gwas_snp_table = pd.read_csv(snakemake.input.raw_sum_stats, sep='\t',
                                 usecols=['chr', 'pos', 'ref', 'alt', 'pval_EUR'],
                                 dtype= str)
    gwas_snp_table = gwas_snp_table.astype(str)
    # build variant id to merge on
    gwas_snp_table['variant'] = gwas_snp_table['chr'] + ':' + gwas_snp_table['pos'] + ':' + gwas_snp_table['ref'] + ':' + gwas_snp_table['alt']
    gwas_snp_table.rename(columns={"pval_EUR": "pval"}, inplace=True)
    # remove variants not having a pval for EU samples
    gwas_snp_table = gwas_snp_table[~gwas_snp_table['pval'].str.contains('nan')]
else:
    #we are dealing with MCH or RBC
    gwas_snp_table = pd.read_csv(snakemake.input.raw_sum_stats, sep='\t')
    
# set to number of samples from phenotype manifest (n_cases_EUR)
gwas_snp_table['n_complete_samples'] = nsamples_table.query('TRAIT == @snakemake.params.trait')['NSAMPLES'].item()

# build the variant id used in the summary statistics tables
variant_annotation_table['variant'] = variant_annotation_table['chrom'] + ':' + variant_annotation_table['pos'] + ':' + variant_annotation_table['ref'] + ':' + variant_annotation_table['alt']
# merge tables
snp_merged_annotation = variant_annotation_table.merge(gwas_snp_table, on='variant', how='left')
snp_merged_annotation.dropna(inplace=True)

# remove rows which don't have rsids in the rsid column
rows_to_drop = snp_merged_annotation[~snp_merged_annotation['rsid'].str.contains('rs')].index
snp_merged_annotation.drop(rows_to_drop, inplace=True)


import numpy as np 
snp_merged_annotation['pval'] = np.exp(snp_merged_annotation['pval'].astype('double'))


# save to file
snp_merged_annotation[['rsid', 'pval', 'n_complete_samples']].to_csv(snakemake.output.trait_sum_stats, sep='\t', header=['SNP', 'P', 'N'], index=False)
