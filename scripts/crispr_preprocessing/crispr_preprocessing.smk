

DEP_MAP_DATA = config['datadir'] + '/dep_map_data'


download_paths = {
    "DepMap18Q3": {"gene_effect": "https://ndownloader.figshare.com/files/12704099",
                   "sample_info": "https://ndownloader.figshare.com/files/12704612"},
    "DepMap22Q2": {"gene_effect": "https://ndownloader.figshare.com/files/34990036",
                   "sample_info": "https://ndownloader.figshare.com/files/35020903"}
}




rule download_dep_map_data:
    output: 
        gene_effect = protected(DEP_MAP_DATA + '/{release}_gene_effect.csv'),
        sample_info = protected(DEP_MAP_DATA + '/{release}_sample_info.csv')
    shell: 'wget -c https://ndownloader.figshare.com/files/12704099 -O {output.gene_effect} \n \
            wget -c https://ndownloader.figshare.com/files/12704612 -O {output.sample_info}'
    

rule download_olfactory_genes:
    output: protected(DEP_MAP_DATA + '/olfactory_genes.txt')
    shell: 'wget -c https://raw.githubusercontent.com/kundajelab/coessentiality/master/olfactory_genes.txt -O {output}'
        


rule normalize_and_map_ids:
    input: 
        gene_effect = DEP_MAP_DATA + '/{release}_gene_effect.csv',
        sample_info = DEP_MAP_DATA + '/{release}_sample_info.csv',
        olfactory_genes = DEP_MAP_DATA + '/olfactory_genes.txt'
    output:
        mapped_scores = DEP_MAP_DATA + '/{release}_gene_effect_mapped.tsv'
    script: 'normalize_and_map_ids.py'
    
rule map_olfactory_genes:
    input:
        olfactory_genes = DEP_MAP_DATA + '/olfactory_genes.txt'
    output: 
        olfactory_genes = DEP_MAP_DATA + '/olfactory_genes_ensg.tsv'
    script: 'map_olfactory_genes.py'
