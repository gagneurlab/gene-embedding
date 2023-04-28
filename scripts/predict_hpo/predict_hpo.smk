

list_to_compare = ['combined_STRING', 'combined_STRING_EXP', 'dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding']


rule hpo_preprocessing:
    input:
        gene2pheno = config['datadir'] + '/input_data/hpo/genes_to_phenotype.txt',
        gene_mapping  = config['datadir'] + '/input_data/hpo/hugo_gene_mapping.txt'
    output:
        gene2pheno = config['datadir'] + '/processed_data/hpo/genes_to_phenotype_ensg.tsv'
    script: 'hpo_preprocessing.R'
    

rule predict_hpo_terms:
    input:
        gene2pheno = config['datadir'] + '/processed_data/hpo/genes_to_phenotype_ensg.tsv',
        emb = lambda wildcards: config['datadir'] + '/embedding/combination/' + wildcards.emb + '.tsv'
    output:
        hpo_pred = config['datadir'] + '/processed_data/hpo/term_pred_{emb}.tsv'
    script: "predict_hpo_terms.py"
 

rule compute_pr_curves:
    input:
        hpo_pred = config['datadir'] + '/processed_data/hpo/term_pred_{emb}.tsv'
    output:
        pr_curves = config['datadir'] + '/processed_data/hpo/pr_curves/term_pred_{emb}.tsv'
    script: 'compute_hpo_pr_curve.R'
    
    
rule plot_pr_curves:
    input:
        hpo_pred = config['datadir'] + '/processed_data/hpo/term_pred_{emb}.tsv'
    output:
        plotpath = config['datadir'] + '/figures/hpo/{emb}_pr_curve.png',
        fold_figure = config['datadir'] + '/figures/hpo/{emb}_fold_figure.png'
    script: 'plot_hpo_pr_curve_new.R'


rule compute_auc:
    input:
        hpo_pred = config['datadir'] + '/processed_data/hpo/term_pred_{emb}.tsv'
    output:
        pr_curves = config['datadir'] + '/processed_data/hpo/auc/term_pred_{emb}.tsv'
    script: 'compute_hpo_auc.R'

    
rule plot_emb_comparisons:
    input:
        pr_curves = expand(config['datadir'] + '/processed_data/hpo/pr_curves/term_pred_{emb}.tsv', 
                           emb = list_to_compare
                           )
    output:
        hpo_comparison = config['datadir'] + '/figures/hpo/comparison_figure.png'
    run:
        touch(output)
        

rule plot_emb_comparisons_auc:
    input:
        pr_curves = expand(config['datadir'] + '/processed_data/hpo/auc/term_pred_{emb}.tsv', 
                           emb = list_to_compare
                           )
    output:
        hpo_comparison = config['datadir'] + '/figures/hpo/comparison_figure_auc.png'
    run:
        touch(output)