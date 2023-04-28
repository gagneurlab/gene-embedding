

rule download_gtex:
    output:
        gtex_data_path = config['datadir'] + '/input_data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz'
    shell: "wget -c https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz -O {output.gtex_data_path}"



rule preprocess_gtex:
    input:
        gtex_data_path = config['datadir'] + '/input_data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz',
        which_genes = config['datadir'] + '/input_data/gene_annotation/mart_export.tsv'
    output:
        gtex_subset_path = config['datadir'] + '/processed_data/gtex/log10_tpm_subset.gz',
    script: "preprocess_gtex.py"
