# MAGMADIR = config['datadir'] + '/magma/'
MAGMADIR = 'data/magma/'
EMBED_DIR = config['datadir'] + '/embedding/combination/'

## ----- NEW ----- ##
# rule all:
#    input:
#        config['datadir'] + "/gwas_prediction/POPSMat/{study}/{study}_pred_group_cv.tsv"

rule predictMagmaEmbCovLDsklearnCV:
    input:
        emb = EMBED_DIR + "{emb}.tsv",
        magma_output = MAGMADIR + "{study}/{study}.genes.out",
        magma_gene_raw = MAGMADIR + "{study}/{study}.genes.raw"
    output:
        pred = config['datadir'] + "/gwas_prediction/{emb}/{study}/{study}_pred.tsv"
    script:
        "predict_magma_emb_cov_LD_sklearnCV.py"


rule predictMagmaEmbCovLDsklearnCVGroupCV:
    input:
        emb = EMBED_DIR + "{emb}.tsv",
        magma_output = MAGMADIR + "{study}/{study}.genes.out",
        magma_gene_raw = MAGMADIR + "{study}/{study}.genes.raw"
    output:
        pred = config['datadir'] + "/gwas_prediction/{emb}/{study}/{study}_pred_group_cv.tsv"
    script:
        "predict_magma_emb_cov_LD_sklearnCV_GroupSplit.py"


# rule predict_magma_emb_cov_LD_sklearnCV_GroupSplitPoPSMat:
#     input:
#         emb = "data/input_data/PoPS/PoPS.features.txt.gz",
#         magma_output = MAGMADIR + "{study}/{study}.genes.out",
#         magma_gene_raw = MAGMADIR + "{study}/{study}.genes.raw"
#     output:
#         pred = config['datadir'] + "/gwas_prediction/POPSMat/{study}/{study}_pred_group_cv.tsv"
#     script:
#         "predict_magma_emb_cov_LD_sklearnCV_GroupSplitPoPSMat.py"
        
        
rule predict_magma_emb_cov_LD_sklearnCV_GroupSplitPoPSMat_CHR:
    input:
        emb = "data/input_data/PoPS/PoPS.features.txt.gz",
        magma_output = MAGMADIR + "{study}/{study}.genes.out",
        magma_gene_raw = MAGMADIR + "{study}/{study}.genes.raw"
    output:
        pred = config['datadir'] + "/gwas_prediction/POPSMat/{study}/{study}_pred_{chrom}_group_cv.tsv"
    script:
        "predict_magma_emb_cov_LD_sklearnCV_GroupSplitPoPSMat_perCHR.py"
        
rule aggregate_pops_mat_predictions:
    input:
        expand(config['datadir'] + "/gwas_prediction/POPSMat/{{study}}/{{study}}_pred_{chrom}_group_cv.tsv", chrom = range(1,23))
    output:
         pred = config['datadir'] + "/gwas_prediction/POPSMat/{study}/{study}_pred_group_cv.tsv"
    script: 'combine_pops_pred.R'

# POPS_COMP_EMB = ['STRING_NO_LIT.crispr_emb.gtex_emb.prot_t5_128_embedding']
POPS_COMP_EMB = ['crispr_emb.gtex_emb.prot_t5_128_embedding']
        
rule plot_PoPS_emb_comparison:
    input:
        magma = expand(MAGMADIR + "{study}/{study}.genes.out", study = STUDIES),
        pops = expand(config['datadir'] + "/PoPS/pops_results/{study}/{study}", study = STUDIES),
        emp_pred = expand(config['datadir'] + "/gwas_prediction/{emb}/{study}/{study}_pred.tsv",
                          emb = POPS_COMP_EMB,
                          study = STUDIES
                          ) 
    output:
        config['datadir'] + '/figures/gwas_prediction/{emb}_pops_embedding_comparison.png'
    script: 'pops_comparison.R'


###
# ------- Elastic Net ---------
###

rule predict_magma_emb_cov_LD_sklearnCV_GroupSplitPoPSMat_CHR_ElasticNet:
    input:
        embedding_path = lambda wildcards: config['datadir'] + '/embedding/combination/' + wildcards.emb + '.tsv',
        # emb = "data/input_data/PoPS/PoPS.features.txt.gz",
        magma_output = MAGMADIR + "{study}/{study}.genes.out",
        magma_gene_raw = MAGMADIR + "{study}/{study}.genes.raw"
    output:
        pred = config['datadir'] + "/gwas_prediction/POPSMat_ElasticNet/{study}/{study}_pred_{chrom}_group_cv.tsv"
    script:
        "predict_magma_emb_cov_LD_sklearnCV_GroupSplitPoPSMat_perCHR_ElasticNet.py"
        
rule aggregate_pops_mat_predictions_ElasticNet:
    input:
        expand(config['datadir'] + "/gwas_prediction/POPSMat_ElasticNet/{{study}}/{{study}}_pred_{chrom}_group_cv.tsv", chrom = range(1,23))
    output:
         pred = config['datadir'] + "/gwas_prediction/POPSMat_ElasticNet/{study}/{study}_pred_group_cv.tsv"
    script: 'combine_pops_pred.R'