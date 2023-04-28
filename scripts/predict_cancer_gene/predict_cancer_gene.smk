cancer



rule preprocess_emogi:
    input: config['datadir'] + "/input_data/cancer_eval/emogi_predictions.tsv"
    output: config['datadir'] + "/cancer_pred/emogi_preprocessed.tsv"
    script: "preprocess_emogi.R"
    
    
rule preprocess_onco_var:
    input: config['datadir'] + "/input_data/cancer_eval/TCGA.PanCancer.all.genes.OncoVar.tsv.gz"
    output: config['datadir'] + "/cancer_pred/onco_var_preprocessed.tsv"
    script: "preprocess_onco_var.R"
     

rule cancer_pred:
    input: 
        model = config['datadir'] + "/cancer_pred/{model}_preprocessed.tsv",
        emb = config['datadir'] + "/embedding/combination/{emb}.tsv"
    output: 
        pred = config['datadir'] + "/cancer_pred/{emb}/{model}_prediction.tsv"
    script: "cancer_pred.py"
    
rule plot_cancer_res:
    input:
        emogi = config['datadir'] + "/cancer_pred/{emb}/emogi_prediction.tsv",
        onco_var = config['datadir'] + "/cancer_pred/{emb}/onco_var_prediction.tsv"
    output:
        config['datadir'] + "/figures/cancer/{emb}/pr_curve.png"
    script: "plot_cancer_pr_curve.R"