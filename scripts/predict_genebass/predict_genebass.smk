# Snakemake pipeline to compute metrics and plot comparison of optimized gradient boosting trees for the trait associated gene prediction

LABEL = "disease_pred_genebass"
PLOT_DIR = config['datadir'] + '/figures/trait_gene_pred/'


# GENEBASS_SELECT_TRAITS 
SELECTTRAITS = ['Alkaline_phosphatase', 'Apolipoprotein_A', 'Basal_metabolic_rate', 'Corneal_resistance_factor_right_custom', 'Creatinine', 'Cystatin_C', 'Eosinophill_count', 'Frequency_of_memory_loss_due_to_drinking_alcohol_in_last_year', 'Glycated_haemoglobin_(HbA1c)', 'HDL_cholesterol', 'IGF-1', 'Mean_corpuscular_haemoglobin', 'Monocyte_count', 'Neutrophill_count', 'Platelet_count', 'Recent_changes_in_speed/amount_of_moving_or_speaking', 'Red_blood_cell_(erythrocyte)_count', 'Red_blood_cell_(erythrocyte)_distribution_width', 'Reticulocyte_count', 'SHBG', 'Standing_height', 'Total_bilirubin', 'Triglycerides', 'Trunk_fat-free_mass', 'Whole_body_fat_free_mass_adjusted_for_body_mass_index_(WBFFMadjBMI)', 'Whole_body_water_mass']

LABELDICT = {
    "dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding": "Embedding",
    "combined_STRING": "STRING (Node2Vec + VERSE)",
    "dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_STRING_embedding": "Embedding + STRING",
    "combined_STRING_EXP":"STRING Experimental (Node2Vec + VERSE)",
    "dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_STRING_EXP_embedding": "Embedding + STRING_EXP",
}

COLOR_DICT = {
    "combined_STRING_EXP": "#a6cee3",
    "dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding": "#1f78b4",
}


# define paths
DATA_DIR = config['datadir'] + "/"
INPUT_DIR = "data/input_data/disease_gene_eval/"
#EMBED_DIR = DATA_DIR + "AE_Embeddings/combinations/"
EVAL_DIR = DATA_DIR + "Evaluation/"
OPT_DIR = DATA_DIR + "Evaluation/hyperopt/"
DIST_DIR = EVAL_DIR + "distance_values/"
PRED_DIR = OPT_DIR+ "predictions/"
PREREC_DIR = OPT_DIR + "prerec_values_genebass500k_burden_or_skato/"  # !!!!!!!!!!!!!! MODIFIED !!!!!!!!!!!!!!!!!!
PLOT_DIR = DATA_DIR + "figures/trait_gene_pred/"


    
rule compute_folds:
    input:
        coding_genes_path = config['datadir'] + "/input_data/gene_annotation/mart_export.tsv", # list of protein coding genes exportet from biomart
        disease_gene_path = INPUT_DIR + "disease_gene_table_genebass500k_burden_or_skato.tsv",
        # label_table_path = INPUT_DIR + "short_trait_labels_genebass500k.tsv",
    output:
        fold_splits_path = config['datadir'] + "/modeling_splits/disease_gene_splits_genebass500k_burden_or_skato.tsv"
    script:
        "compute_folds.py"
        
        
rule predict_disease_gene:
    input:
        fold_splits_path = config['datadir'] + "/modeling_splits/disease_gene_splits_genebass500k_burden_or_skato.tsv",
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
    
     
