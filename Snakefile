# Snakefile


### the following embeddings will be concatenated to form the baseline combinations
list_individual_embeddings = [
    #'STRING_NO_LIT',
    'STRING_NO_LIT2',
    'STRING_NFC',
    'STRING_EXP',
    #'STRING_NO_BIAS',
    #'crispr_emb',
    'DepMap18Q3_crispr_emb',
    #'DepMap22Q2_crispr_emb',
    'prot_t5_128_embedding',
    'gtex_emb'
]

# list_individual_embeddings = [
#     'STRING'
# ]



import glob

configfile: "config/config.yaml"

# print(config)
# import ipdb; ipdb.set_trace()

include: 'scripts/graph_preprocessing/graph_preprocessing.smk'
include: 'scripts/graph_embedding/embedding.smk'
include: 'scripts/crispr_preprocessing/crispr_preprocessing.smk'
include: 'scripts/preprocess_gwas_data/preprocess_gwas_data.smk'
include: 'scripts/create_joint_embeddings/create_joint_embeddings.smk'
include: 'scripts/embedding_gtex/embedding_gtex.smk'
include: 'scripts/misc_figures/misc_figures.smk'
include: 'scripts/embedding_vtf/embedding_vtf.smk'


#include: 'scripts/other_preprocessing/disease_gene_preprocessing.smk'



### list additional embeddings to be included in the figures.
list_additional_embeddings = [
#    "combined_STRING",
    "ZEROS",
    "dtf_depMap_gtex_happy-oath-195",
    "vtf_most_lucky-galaxy-20_embedding_200_epoch",
   # "dtf_model_depMap_480_480", "dtf_model_depMap_64_64",
   # "dtf_depMap_gtex_dashing-glade-21_2516141_embedding",
   # "dtf_depMap_gtex_expert-firefly-74_2591521_embedding",
   # "dtf_network_sleek-sound-16_2610004_embedding",
    # "dtf_depMap_gtex_perturbSeq_happy-vortex-75_2591534_embedding",
    # "dtf_most_firm-glade-998_2675299_embedding",
    # "dtf_most_TS_crimson-donkey-1040_2758658_embedding",
    # "vtf_depMap_gtex_perturbSeq_lemon-glade-15_2901568_embedding",
    # "vtf_most-distinctive-cloud-29_3015554_embedding_250_epoch",
    "vtf_most_ts_test_lucky-violet-24_2982432_embedding_50_epoch",
    "vtf_most_ts_test_laced-snowball-53_3099270_embedding_50_epoch",
    "vtf_gtex_depMap_lyric-brook-22_3222375_embedding",
    "vtf_most_laced-universe-11_3222391_embedding",
    "vtf_gtex_depMap_stellar-frog-91_3223981_embedding",
    "vtf_string_avid-serenity-84_3249168_embedding",
    "dtf_depMap_dirichlet_helpful-shadow-98_embedding",
    "vtf_depMap_dirichlet_cerulean-glitter-115_3266661_embedding",
    "vtf_string_dirichlet_pleasant-shadow-183_3269363_embedding",
    "vtf_gtex_dirichlet_zesty-valley-60_3274786_embedding",
    "vtf_string_dirichlet_ethereal-firefly-53_3277825_embedding",
    "vtf_string_lucky-hill-30_3277825_embedding",
    "vtf_gtex_dirichlet_comfy-water-12_3278013_embedding_150_epoch",
    "vtf_gtex_TS_endothelial_honest-cosmos-1169_3313125_embedding",
    "vtf_gtex_ts_many_hopeful-pond-1166_3313082_embedding_50_epoch",
    "vtf_gtex_ts_organs_splendid-butterfly-17_3341650_embedding_50_epoch",
    "vtf_gtex_humap_prott5_2ts_peach-galaxy-201_3381543_embedding",
    "vtf_gtex_depMap_protT5_humap_mild-bird-236_3381517_embedding",
    "vtf_gtex_depMap_huMap_prott5_ts_organs_different-resonance-134_3381568_embedding",
    "vtf_gtex_depMap_prott5_perturbSeq_glad-feather-351_3382939_embedding",
    "vtf_gtex_depMap_prott5_perturbSeq_rose-mountain-430_3383575_embedding",
    "vtf_gtex_depMap_prott5_perturbSeq_polished-sun-512_3384348_embedding",
    "vtf_gtex_depMap_perturbSeq_protT5_huMap_string_ts_organs_fancy-dream-547_3384392_embedding_100_epoch",
    "vtf_gtex_depMap_perturbSeq_protT5_EMB_DIM_512_apricot-elevator-579_3404825_embedding",
    "vtf_gtex_depMap_perturbSeq_protT5_huMap_string_ts_organs_elated-snowflake-641_3415579_embedding_100_epoch",
    "vtf_gtex_depMap_perturbSeq_protT5_huMap_string_ts_organs_elated-snowflake-641_3415579_embedding_200_epoch",
    "vtf_gtex_depMap_perturbSeq_protT5_comfy-cherry-606_3404825_embedding",
    "vtf_gtex_depMap_perturbSeq_protT5_huMap_string_ts_organs_trim-snowball-629_3415577_embedding",
    "vtf_gtex_depMap_perturbSeq_protT5_royal-plasma-721_3507179_embedding",
    "vtf_gtex_depMap_perturbSeq_protT5_playful-donkey-729_3507180_embedding",
    "vtf_gtex_depMap_perturbSeq_protT5_absurd-river-752_3510204_embedding",
    "dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding",
    "dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_humap_embedding",
    "dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_STRING_embedding",
    "dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_Co_evolution_embedding",
    "vtf_gtex_depMap_protT5_huMap_worldly-spaceship-5_4456762_embedding",
    "dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_STRING_EXP_embedding",
    "vtf_gtex_huMap_depMap_coEvo_protT5_youthful-music-439_4520579_embedding",
    "vtf_gtex_depMap_protT5_huMap_sleek-durian-664_4520403_embedding",
    #"pops_mat",
    #"pops_exp_mat",
]

all_embeddings = all_combination_list + list_additional_embeddings 
print("# ----------------- List of all embeddings ------------------------ ")
for emb in all_embeddings:
    print(f"# {emb}")
print("#  ----------------------------------------------------------------- ")


include: 'scripts/eval_embedding/eval_emb.smk'
include: 'scripts/predict_hpo/predict_hpo.smk'
include: 'scripts/predict_gwas/predict_gwas.smk'
include: 'scripts/predict_disease_gene/predict_disease_gene.smk'
#include: 'scripts/predict_cancer_gene/predict_cancer_gene.smk'




GENE_PREDICTION = config['datadir'] + '/gene_prediction'
GRAPH_DATA = config['datadir'] + '/graphdata'
SIMULATION = config['datadir'] + '/simulation3'

# print(all_combination)

import numpy as np
dims = np.concatenate([2**np.arange(3,9), np.array([24])])
dims = [24, 128]


rule all:
    input:  
        config['datadir'] + "/figures/trait_gene_pred/trait_gene_pred_individaul.png",
        #config['datadir'] + "/figures/trait_gene_pred/trait_gene_pred_all_emb.png",
        config['datadir'] + '/figures/hpo/comparison_figure_auc.png',
        # FINAL FIGURE PATHS
        
        # ## DIM SEARCH
        # expand(config['datadir'] + '/figures/{data}/dim_search.png',
        #        data = ['gtex', 'crispr_DepMap18Q3', 'crispr_DepMap22Q2']
        #        ),
        
        ## TRAIT GENE PREDICTION
        config['datadir'] + "/figures/trait_gene_pred/trait_gene_pred_individaul.png",
        # config['datadir'] + "/figures/trait_gene_pred/trait_gene_pred_stratified.png",
        # config['datadir'] + "/figures/trait_gene_pred/trait_gene_pred_individaul_crispr.png",
        # config['datadir'] + "/figures/trait_gene_pred/trait_gene_pred_individaul_clip.png",
        # config['datadir'] + "/figures/trait_gene_pred/trait_gene_pred_individaul_no_string_emb.png",
        
        ## UMAP
        # expand(PLOT_DIR + "umap/{emb}/trait_gene_umap.png", emb = all_embeddings), 
        
        ## CANCER
        expand(config['datadir'] + "/figures/cancer/{emb}/pr_curve.png",
               emb = all_embeddings
               ),
        
        ## HPO
        config['datadir'] + '/figures/hpo/comparison_figure.png',
        
        
#         ## STRING eval dist
#         expand(config['datadir'] + '/figures/eval_dist/{emb}_boxplot.pdf',
#                emb = all_embeddings
#                ),
        
#         config['datadir'] + '/paper_numbers.html',
        
        
#         # GWAS Figure
#         # expand(config['datadir'] + '/figures/gwas_prediction/{emb}_pops_embedding_comparison.png',
#         #        emb = all_embeddings
#         #        ),
           
            
#         # # GWAS Intermediates.
#         # expand("data/PoPS/pops_results/{study}/{study}.{chr}.results", 
#         #        study=STUDIES, 
#         #        chr=CHROMOSOMES
#                # ),
#         # expand(config['datadir'] + "/gwas_prediction/{emb}/{study}/{study}_pred.tsv",
#         #        study=STUDIES,
#         #        emb = ['STRING_NO_LIT',
#         #               'STRING_NO_LIT.crispr_emb.gtex_emb.prot_t5_128_embedding',
#         #               'crispr_emb.gtex_emb.prot_t5_128_embedding',
#         #               "vtf_gtex_depMap_perturbSeq_protT5_EMB_DIM_512_apricot-elevator-579_3404825_embedding",
#         #               'vtf_no_assay_var_pert3_dim_256_vtf_model_depMap.gtex.protT5_run_1_layers_3_embedding',
#         #               'vtf_no_assay_var_pert3_dim_256_vtf_model_depMap.gtex.perturbSeq3.protT5_run_1_layers_3_embedding',
#         #               'combined_gtex_emb',
#         #               ]
#         #        ),
        
#         # expand(config['datadir'] + "/gwas_prediction/{emb}/{study}/{study}_pred_group_cv.tsv",
#         #        study=STUDIES,
#         #        emb = ['STRING_NO_LIT',
#         #               'STRING_NO_LIT.crispr_emb.gtex_emb.prot_t5_128_embedding',
#         #               'crispr_emb.gtex_emb.prot_t5_128_embedding', 
#         #               "vtf_gtex_depMap_perturbSeq_protT5_EMB_DIM_512_apricot-elevator-579_3404825_embedding",
#         #               'vtf_no_assay_var_pert3_dim_256_vtf_model_depMap.gtex.protT5_run_1_layers_3_embedding',
#         #               'vtf_no_assay_var_pert3_dim_256_vtf_model_depMap.gtex.perturbSeq3.protT5_run_1_layers_3_embedding',
#         #               'dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding',
#         #               'combined_gtex_emb',
#         #               'ZEROS',
#         #              ]
#         #        ),      
#         # # Predict emb using POPS matrix
#         # expand(config['datadir'] + "/gwas_prediction/POPSMat/{study}/{study}_pred_group_cv.tsv",
#         #        study=STUDIES
#         #        ),
           
            
# #         expand(config['datadir'] + "/gwas_prediction/{emb}/ElasticNet/{study}/{study}_pred_group_cv.tsv",
# #                study=STUDIES,
# #                emb = ['STRING_NO_LIT',
# #                       'STRING_NO_LIT.crispr_emb.gtex_emb.prot_t5_128_embedding',
# #                       'dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding',
# #                       'ZEROS',
# #                      ]
# #                ),      
           
#         # DEP MAP Intermediates
#         DEP_MAP_DATA + '/olfactory_genes_ensg.tsv',    
            
        
            
#         #OTHER
#         #expand(PLOT_DIR + "umap/{emb}/trait_gene_umap.png", emb = all_combination),   
        
#         #config['datadir'] + '/processed_data/gtex/find_dim_res_median.tsv'
#         DEP_MAP_DATA + '/pca_residuals.tsv',
        
#         # config['datadir'] + '/embedding/individual/HUMAP_HURI_PROPER.tsv',
        


rule plot_cancer:
    input:
        expand(config['datadir'] + "/figures/cancer/{emb}/pr_curve.png",
               emb = all_combination
               ),


rule plot_hpo:
    input:
        expand(config['datadir'] + '/figures/hpo/{emb}_pr_curve.png',
               emb = all_embeddings
               ),
        expand(config['datadir'] + '/figures/hpo/{emb}_fold_figure.png',
               emb = all_embeddings
               ),
        
rule plot_umap_target:
    input:
        expand(PLOT_DIR + "umap/{emb}/trait_gene_umap.png", 
               emb = [
                    # "vtf_depMap_gtex_perturbSeq_lemon-glade-15_2901568_embedding",
                    # "vtf_most-distinctive-cloud-29_3015554_embedding_250_epoch",
                    # "dtf_model_depMap_128_128", "dtf_model_depMap_16_16", 
                    # "vtf_most_ts_test_laced-snowball-53_3099270_embedding_50_epoch",
                    # "dtf_model_depMap_256_256", "dtf_model_depMap_32_32", 
                    # "dtf_depMap_swept-wildflower-15_2493617", 
                    # "dtf_depMap_gtex_happy-oath-195",
                    # "dtf_most_firm-glade-998_2675299_embedding",
                    # "dtf_depMap_peachy-brook-10_2492927",
                    # "dtf_depMap_gtex_swept-butterfly-10_2498062",
                    # "dtf_most_TS_crimson-donkey-1040_2758658_embedding",
                    # "dtf_network_sleek-sound-16_2610004_embedding",
                    # "dtf_depMap_lyric-wind-9_2492568", 
                    # "dtf_networks_jolly-spaceship-5_2524221_embedding_100_epoch",
                    # "dtf_depMap_gtex_summer-energy-20_2516150_embedding_99_epoch",
                    # "vtf_depMap_gtex_summer-galaxy-10_2846270_embedding",
                    # "dtf_most_gentle-puddle-1035_2732599_embedding",
                    # "dtf_depMap_gtex_unique-silence-2_2486065", 
                    # "dtf_depMap_gtex_dashing-glade-21_2516141_embedding", 
                    # "dtf_depMap_gtex_perturbSeq_wild-durian-25_2521985_embedding_100_epoch", "dtf_depMap_gtex_fragrant-grass-17_2513316_embedding",
                    # "vtf_most_ts_test_lucky-violet-24_2982432_embedding_50_epoch",
                    # "dtf_all_spring-vortex-89_2624879_embedding",
                    # "vtf_most_lucky-galaxy-20_embedding_200_epoch",          
                    # "dtf_depMap_gtex_expert-firefly-74_2591521_embedding", 
                    # "dtf_depMap_gtex_perturbSeq_happy-vortex-75_2591534_embedding",
                   "vtf_gtex_depMap_lyric-brook-22_3222375_embedding",
               ]
        ), 
        
rule plot_disease_gene_target:
    input:
        config['datadir'] + "/figures/trait_gene_pred/trait_gene_pred_individaul.png",
        config['datadir'] + "/figures/trait_gene_pred/trait_gene_pred_stratified.png",
