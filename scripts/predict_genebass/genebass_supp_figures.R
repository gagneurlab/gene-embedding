library(data.table)
library(cowplot)
library(ggplot2)
library(dplyr)
library(yaml)

config <- read_yaml('config/config.yaml')
DATADIR <- config$datadir[1]

trait_genes_dir <- file.path(DATADIR,"/input_data/disease_gene_eval/disease_gene_table_genebass500k_burden_or_skato.tsv")
trait_genes = fread(trait_genes_dir)
gene_counts <- raw_data %>% count(disease)
# gene_counts[order(-n)]

# SELECT_TRAITS = c("Total_bilirubin", "Mean_corpuscular_haemoglobin", "Mean_sphered_cell_volume", "Red_blood_cell_(erythrocyte)_distribution_width", "Glycated_haemoglobin_(HbA1c)", "Red_blood_cell_(erythrocyte)_count", "Reticulocyte_count", "Platelet_distribution_width", "Mean_platelet_(thrombocyte)_volume", "Platelet_count", "Lymphocyte_count", "Monocyte_count", "Eosinophill_count", "Lymphocyte_percentage", "Neutrophill_count", "Sitting_height", "Standing_height", "Forced_expiratory_volume_in_1-second_(FEV1),_predicted", "Whole_body_fat-free_mass", "Forced_vital_capacity_(FVC)", "Frequency_of_inability_to_cease_drinking_in_last_year", "Frequency_of_memory_loss_due_to_drinking_alcohol_in_last_year", "Frequency_of_failure_to_fulfil_normal_expectations_due_to_drinking_alcohol_in_last_year", "Recent_changes_in_speed/amount_of_moving_or_speaking", "Cholesterol", "HDL_cholesterol", "Triglycerides", "Alkaline_phosphatase", "Gamma_glutamyltransferase", "Creatinine", "Cystatin_C", "Urate", "SHBG", "IGF-1")
SELECT_TRAITS <- config$SELECT_TRAITS$GENEBASS

gene_counts[,selected:=disease %in% SELECT_TRAITS]
gene_counts <- gene_counts[order(n)]
gene_counts$disease = factor(gene_counts$disease, level=gene_counts$disease)


# Number of genes associated with traits
num_genes <- ggplot(gene_counts, aes(x=disease, y=n)) +
  geom_bar(stat="identity", fill = "#1561A9") +
  coord_flip() +
  xlab("Trait") +
  ylab("Number of associated genes") +
  theme_cowplot(font_size = 11)
num_genes



#  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

SELECT_TRAITS = c("Total_bilirubin", "Mean_corpuscular_haemoglobin", "Mean_sphered_cell_volume", "Red_blood_cell_(erythrocyte)_distribution_width", "Glycated_haemoglobin_(HbA1c)", "Red_blood_cell_(erythrocyte)_count", "Reticulocyte_count", "Platelet_distribution_width", "Mean_platelet_(thrombocyte)_volume", "Platelet_count", "Lymphocyte_count", "Monocyte_count", "Eosinophill_count", "Lymphocyte_percentage", "Neutrophill_count", "Sitting_height", "Standing_height", "Forced_expiratory_volume_in_1-second_(FEV1),_predicted", "Whole_body_fat-free_mass", "Forced_vital_capacity_(FVC)", "Frequency_of_inability_to_cease_drinking_in_last_year", "Frequency_of_memory_loss_due_to_drinking_alcohol_in_last_year", "Frequency_of_failure_to_fulfil_normal_expectations_due_to_drinking_alcohol_in_last_year", "Recent_changes_in_speed/amount_of_moving_or_speaking", "Cholesterol", "HDL_cholesterol", "Triglycerides", "Alkaline_phosphatase", "Gamma_glutamyltransferase", "Creatinine", "Cystatin_C", "Urate", "SHBG", "IGF-1")

EMBEDDINGS <- c('dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding', 'combined_STRING_EXP', 'combined_STRING')

EMB_SHORT = c("Omics", "STRING Exp", "STRING")

trait_genes_dir <- file.path(DATADIR,"/input_data/disease_gene_eval/disease_gene_table_genebass500k_burden_or_skato.tsv")
trait_genes = fread(trait_genes_dir)

eval_dir <- file.path(DATADIR,"/Evaluation/hyperopt")

# Get average precision from file --------------------------------
get_avg_precision <- function(trait, emb){
  eval_file <- sprintf("%s/prerec_values_genebass500k_burden_or_skato/%s/%s_stratify_num_pub_False_eval.tsv", eval_dir, trait, emb)
  stats <- fread(file = eval_file, sep = '\t', header = TRUE)
  avg = mean(stats[['average_precision']])
  return(avg)
}

get_min_auprc <- function(trait){
  eval_file <- sprintf("%s/prerec_values_genebass500k_burden_or_skato/%s/combined_STRING_stratify_num_pub_False_eval.tsv", eval_dir, trait)
  stats <- fread(file = eval_file, sep = '\t', header = TRUE)
  min_auprc = mean(stats[['min_auprc']])
  return(min_auprc)
}

avg_prec <- matrix(nrow = length(SELECT_TRAITS), ncol = length(EMB_SHORT))
colnames(avg_prec) <- EMB_SHORT

auc_dt_list <- lapply(1:length(SELECT_TRAITS), FUN = function(trait) { 
  emb_auc_list <- lapply(1:length(EMBEDDINGS), FUN = function(emb) {
    get_avg_precision(SELECT_TRAITS[trait],EMBEDDINGS[emb])
  })
  emb_dt <- do.call(data.table, emb_auc_list)
  emb_dt[, `:=` (Random = get_min_auprc(SELECT_TRAITS[trait]), 
                 Trait = SELECT_TRAITS[trait])]
})
avg_prec <- rbindlist(auc_dt_list)

colnames(avg_prec) <- c("Omics", "STRING Exp", "STRING", "Random", "Trait")
avg_prec <- avg_prec[,c("Omics","STRING Exp","STRING","Trait")]
dt_gb <- melt(avg_prec, variable.name = "Embedding", value.name = "auprc")

hm <- ggplot(dt_gb, aes(Embedding, Trait)) + 
  geom_tile(aes(fill = auprc)) + 
  labs(fill = "auPRC") +
  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.justification = c("left", "top"),
        legend.box.just = "left",
        axis.title.y=element_blank(), 
        axis.title.x=element_blank()) 
# guides(x = guide_axis(angle = 45))
hm
