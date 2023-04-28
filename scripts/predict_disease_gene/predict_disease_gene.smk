# Snakemake pipeline to compute metrics and plot comparison of optimized gradient boosting trees for the trait associated gene prediction


# do we need the next line?

LABEL = "disease_pred"

SELECTTRAITS = ['mito', 'Neurology', 'neuromuscular', 'Ophthalmology', 'IEI']

PLOT_DIR = config['datadir'] + '/figures/trait_gene_pred/'

SELECTTRAITS = ['epilepsy', 'mito', 'Neurology', 'neuromuscular', 'Ophthalmology', 'IEI']


# move label dict to general Snakefile, so names are consistent across tasks.

LABELDICT = {
    "combined_STRING": "STRING",
    "combined_STRING_EXP": "STRING Exp.",
    "dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding": "Omics",
    "pops_mat": "PoPS",
    "pops_exp_mat": "PoPS Exp",
    #"dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_STRING_embedding": "Omics + STRING",
    #"dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_STRING_EXP_embedding": "Omics + STRING Exp."
    }

COLOR_DICT = {
    "combined_STRING": "#1f78b4",
    "combined_STRING_EXP": "#a6cee3",
    "dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding":  '#ff7f00',
    "pops_mat": "#33a02c",
    "pops_exp_mat": '#b2df8a',
    #"dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_STRING_embedding": '#6a3d9a',
    #"dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_STRING_EXP_embedding": '#cab2d6'
}


# define paths
DATA_DIR = config['datadir'] + "/"
INPUT_DIR = "data/input_data/disease_gene_eval/"
#EMBED_DIR = DATA_DIR + "AE_Embeddings/combinations/"
EVAL_DIR = DATA_DIR + "Evaluation/"
OPT_DIR = DATA_DIR + "Evaluation/hyperopt/"
DIST_DIR = EVAL_DIR + "distance_values/"
PRED_DIR = OPT_DIR+ "predictions/"
PREREC_DIR = OPT_DIR + "prerec_values/"
PLOT_DIR = DATA_DIR + "figures/trait_gene_pred/"

# automatic input modifications
order = all_embeddings
#order.insert(0, "pca-gtex_pca-crispr_norm-stringnolit_norm-stringnolitverse_norm-proteinsmall")
# print(order)

    
rule compute_folds:
    input:
        coding_genes_path = config['datadir'] + "/input_data/gene_annotation/mart_export.tsv", # list of protein coding genes exportet from biomart
        disease_gene_path = INPUT_DIR + "disease_gene_table.tsv",
        label_table_path = INPUT_DIR + "short_trait_labels.tsv",
    output:
        fold_splits_path = config['datadir'] + "/modeling_splits/disease_gene_splits.tsv"
    script:
        "compute_folds.py"
        
rule exp_predict_disease_gene_pops_target:
    input:
        prediction_results_path = expand(PRED_DIR + "{disease}/pops_exp_mat.tsv", disease = SELECTTRAITS)
        
rule predict_disease_gene_pops:
    input:
        fold_splits_path = config['datadir'] + "/modeling_splits/disease_gene_splits.tsv",
        embedding_path = "data/input_data/PoPS/PoPS.features.txt.gz",
    output: 
        prediction_results_path = PRED_DIR + "{disease}/pops_exp_mat.tsv"
    script: "predict_disease_gene_pops.py"        

        
rule predict_disease_gene:
    input:
        fold_splits_path = config['datadir'] + "/modeling_splits/disease_gene_splits.tsv",
        embedding_path = lambda wildcards: config['datadir'] + '/embedding/combination/' + wildcards.emb + '.tsv'
    output: 
        prediction_results_path = PRED_DIR + "{disease}/{emb}.tsv"
    script: "predict_disease_gene.py"
    

rule compute_average_precision:
    input:
        prediction_results_path = PRED_DIR + "{disease}/{emb}.tsv"
    output:
        eval_table_path = PREREC_DIR + "{disease}/{emb}_stratify_num_pub_{stratify_num_pub}_eval.tsv"
    script: "compute_performance_metrics.py"
    
     
rule plot_disease_gene:
    input:
        predictions = expand(PREREC_DIR + "{disease}/{emb}_stratify_num_pub_False_eval.tsv", disease=SELECTTRAITS, emb=LABELDICT.keys(), stratify_num_pub = [True, False]),
    params:
        select_traits = SELECTTRAITS,
        emb_labels = LABELDICT,
        color_dict = COLOR_DICT,
        #box_order = order
    output:
        plot_path = config['datadir'] + "/figures/trait_gene_pred/trait_gene_pred_individaul.png",
    script: "plot_disease_gene_individual.R"
    # "plot_disease_gene_recall.R" 
    #
    
    
