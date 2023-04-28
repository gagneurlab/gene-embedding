library(ggplot2)
library(ggpubr)
library(cowplot)
library(data.table)
library(latex2exp)
library(rstatix)
library(dplyr)
library(patchwork)
library(png)
library(yaml)

config <- read_yaml('config/config.yaml')
DATADIR <- config$datadir[1]

# TRAITS = c("Alanine_aminotransferase","Cholesterol","IGF_1","Albumin","LDL_direct","Alkaline_phosphatase","C_reactive_protein","Total_protein","Apolipoprotein_A","Creatinine","Triglycerides","Apolipoprotein_B","Mean_corpuscular_haemoglobin","Phosphate","Body_mass_index__BMI_","HDL_cholesterol","SHBG", "Urea","Gamma_glutamyltransferase","Glycated_haemoglobin__HbA1c_","Vitamin_D", "Calcium_30680")
TRAITS <- config$SELECT_TRAITS$GWAS

EMBS = c('dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding', 'combined_STRING_EXP', 'combined_STRING','pops_mat_exp','pops_mat', "ZEROS")

EMB_SHORT = c("Omics", "STRING Exp", "STRING", "PoPS Exp", "PoPS", "Zeros", "TRAIT")

embeddings <- list.files(file.path(DATADIR,"/gwas_prediction/"), recursive = TRUE, pattern = "embedding")

rsq_vals <- matrix(nrow = length(TRAITS), ncol = length(EMBS))


rsq_matrix <- lapply(1:length(TRAITS), FUN = function(trait){
  magma_file <- sprintf('%s/magma/%s/%s.genes.out', DATADIR, TRAITS[trait], TRAITS[trait])
  magma_pred <- fread(magma_file)

  rsq_trait <- lapply(1:length(EMBS), FUN = function(emb){
    pred_file <- sprintf("%s/gwas_prediction/%s_ElasticNet/%s/%s_pred_group_cv.tsv", DATADIR, EMBS[emb], TRAITS[trait], TRAITS[trait])
    pred_tab <- fread(pred_file)

    all_tab <- merge(magma_pred,pred_tab, by='GENE', all.x=TRUE)
    all_tab[, pred := ifelse(is.na(pred), mean(pred, na.rm = TRUE), pred)]
    all_tab[, cor(ZSTAT.x, pred)^2]
  })
  dt <- do.call(data.table, rsq_trait)
  dt[, TRAIT := TRAITS[trait]]
})
rsq_vals <- rbindlist(rsq_matrix)
colnames(rsq_vals) <- EMB_SHORT
# rsq_vals

delta_r2 <- rsq_vals
delta_r2[,1:5] <- delta_r2[,1:5] - delta_r2[,Zeros]

colnames(delta_r2) <- EMB_SHORT

emb_list <- c("Omics", "STRING Exp", "STRING", "PoPS Exp", "PoPS") 
color_map <- c("Omics"="#ff7f00", "STRING Exp"="#a6cee3", "STRING"="#1f78b4", "PoPS Exp"="#b2df8a","PoPS"="#33a02c")

delta_r2 <- delta_r2[, ..emb_list]

stats_plot <- melt(delta_r2, variable.name = "Embedding", value.name = "delta_R2")

stat_pvalue_gwas <- stats_plot %>%
  wilcox_test(delta_R2 ~ Embedding, paired = TRUE) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  filter(p.adj.signif == 'ns') %>%
  add_y_position() %>%
  mutate(y.position = seq(min(y.position), max(y.position), length.out = n()))
# stat_pvalue <- filter(stat_pvalue, p.adj.signif == 'ns')


gwas_box_plot <- ggplot(stats_plot[Embedding %in% emb_list], aes(x = reorder(Embedding, delta_R2, FUN = median), y = delta_R2, color=Embedding)) + 
  geom_boxplot(linewidth=1) + 
  theme_cowplot(font_size = 14) +
  ylab(TeX("$\\Delta R^2$")) + 
  scale_y_continuous(minor_breaks = seq(0, 0.2, 0.01)) + 
  theme(axis.title.x=element_blank(),
        legend.position = "none") +
  background_grid() +
  guides(x = guide_axis(angle = 45)) +
  scale_x_discrete(labels = c("STRING\nExp.", "STRING", "Omics", "PoPS\nExp.", "PoPS")) +
  scale_color_manual(values=color_map) +
  stat_pvalue_manual(stat_pvalue_gwas, label = "p.adj.signif")
  # stat_compare_means(comparisons = all_combinations, method = "wilcox.test", paired = TRUE) #, hide.ns = TRUE)
# gwas_box_plot

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# _____________Add runtime path____________
df_time <- fread(file.path(DATADIR,'/gwas_prediction/Embedding_runtime_GWAS.tsv'), sep = '\t', header = TRUE)
dt <- melt(df_time, variable.name = "Embedding", value.name = "Time")

dt$Embedding <- factor(dt$Embedding, levels=c("PoPS","PoPS Exp","Omics","STRING","STRING Exp", "Omics...STRING","Omics...STRING.Exp")) #,"PoPS_XGBoost"))
dt[,Time_hours:= Time/60]

runtime_plot <- ggplot(dt[Embedding %in% c("PoPS","PoPS Exp","Omics","STRING","STRING Exp")], aes(Embedding, Time_hours, color=Embedding)) + 
  geom_boxplot(linewidth=1) + 
  theme_cowplot(font_size = 14) + 
  # scale_y_log10(breaks=c(1,60,600,6000), 
                # labels=c("1 min","1 hour", "10 hours", "100 hours")) + 
  scale_y_log10() +
  theme(axis.title.x=element_blank(),
        legend.position="none") + 
  ylab("Compute time (hours)") +  
  scale_color_manual(values=color_map) +
  background_grid() +
  annotation_logticks(sides = "l") +
  guides(x = guide_axis(angle = 45)) +
  scale_x_discrete(labels = c("PoPS", "PoPS\nExp.", "Omics", "STRING\nExp.", "STRING"))


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# SELECT_TRAITS = c("Total_bilirubin", "Mean_corpuscular_haemoglobin", "Mean_sphered_cell_volume", "Red_blood_cell_(erythrocyte)_distribution_width", "Glycated_haemoglobin_(HbA1c)", "Red_blood_cell_(erythrocyte)_count", "Reticulocyte_count", "Platelet_distribution_width", "Mean_platelet_(thrombocyte)_volume", "Platelet_count", "Lymphocyte_count", "Monocyte_count", "Eosinophill_count", "Lymphocyte_percentage", "Neutrophill_count", "Sitting_height", "Standing_height", "Forced_expiratory_volume_in_1-second_(FEV1),_predicted", "Whole_body_fat-free_mass", "Forced_vital_capacity_(FVC)", "Frequency_of_inability_to_cease_drinking_in_last_year", "Frequency_of_memory_loss_due_to_drinking_alcohol_in_last_year", "Frequency_of_failure_to_fulfil_normal_expectations_due_to_drinking_alcohol_in_last_year", "Recent_changes_in_speed/amount_of_moving_or_speaking", "Cholesterol", "HDL_cholesterol", "Triglycerides", "Alkaline_phosphatase", "Gamma_glutamyltransferase", "Creatinine", "Cystatin_C", "Urate", "SHBG", "IGF-1")

SELECT_TRAITS <- config$SELECT_TRAITS$GENEBASS

EMBEDDINGS <- c('dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding', 'combined_STRING_EXP', 'combined_STRING')

trait_genes = fread(file.path(DATADIR,'/input_data/disease_gene_eval/disease_gene_table_genebass500k_burden_or_skato.tsv'))

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
# avg_prec <- avg_prec[,c("Omics","STRING Exp","STRING","Trait")]
dt_gb <- melt(avg_prec, variable.name = "Embedding", value.name = "auprc")


stat_pvalue <- dt_gb %>%
  wilcox_test(auprc ~ Embedding, paired = TRUE) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  # filter(p.adj.signif != 'ns') %>%
  add_y_position()
  # mutate(y.position = 0.3)#seq(0.29, 0.35, length.out = n()))


genebass_box_plot <- ggplot(dt_gb[Embedding %in% c("Omics", "STRING Exp", "STRING")], aes(Embedding, auprc, color=Embedding)) + 
  geom_boxplot(linewidth=1) +
  theme_cowplot(font_size = 14) +
  ylab("auPRC") +
  theme(axis.title.x=element_blank(),
        legend.position="none") +
  scale_color_manual(values=color_map) +
  background_grid() +
  # stat_pvalue_manual(stat_pvalue) +
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  scale_x_discrete(labels = c("Omics", "STRING Exp.", "STRING")) +
  geom_text(aes(x=1.5, y=0.31, label="*", color="black")) +
  geom_segment(aes(x = 1, y = 0.3, xend = 2, yend = 0.3, color="black")) +
  geom_segment(aes(x = 1, y = 0.27, xend = 1, yend = 0.3, color="black")) +
  geom_segment(aes(x = 2, y = 0.27, xend = 2, yend = 0.3, color="black"))

# genebass_box_plot



ggarrange(genebass_box_plot, pops_zscores, gwas_box_plot, runtime_plot, ncol = 2, nrow = 2, heights = c(2,2), labels = c("A","B","C","D"))

