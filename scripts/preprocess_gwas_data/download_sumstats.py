import pandas as pd
import requests
 
print('Download Starting...')

# parameters
sources_table_path = snakemake.input['sum_stats_sources']
trait = snakemake.params['trait']
out_path = snakemake.output['raw_sum_stats']


# read table
sources_table = pd.read_csv(sources_table_path, sep='\t')
# get url

print(f"TRAIT: {trait}")
print(sources_table)
url = sources_table.query('TRAIT == @trait')['LINK'].item()
r = requests.get(url)
 
with open(out_path,'wb') as output_file:
    output_file.write(r.content)

print('saved file at ' + out_path)