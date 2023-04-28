library(data.table)
library(ggplot2)
library(glue)
library(patchwork)
library(yaml)

config <- read_yaml('config/config.yaml')
DATADIR <- config$datadir[1]

calculatePrecisionRecall <- function(score, label, decr=TRUE, 
                                     total=sum(label)){
  dt <- data.table(score=score, label=label)
  dt <- dt[order(score, decreasing=decr, na.last=TRUE)]
  dt[,rank      := 1:.N]
  dt[,TP        := cumsum(label)]
  dt[,c('TP', 'rank'):=list(max(TP), max(rank)), by=score]
  dt[,precision := TP/rank]
  dt[,recall    := TP/total]
  dt[,recall_diff := diff(c(0, recall))]
  dt[,au_prc := sum(precision * recall_diff)]
  dt[,recall_diff:=NULL]
  dt
}

prepare_pr_tables <- function(data, group_col) {
  
  pr_tables <- lapply(unique(data[[group_col]]), function(group) {
    dt <- data[get(group_col) == group]
    pr <- dt[, calculatePrecisionRecall(pred, label)]
    out <- cbind(pr, dt[, get(group_col)])
    setnames(out, "V2", "group")
    out
  })
  pr_dt <- rbindlist(pr_tables)
}

files <- list.files(file.path(DATADIR,'/processed_data/emogi_pred/'), full.names = TRUE )
names(files) <- basename(files)

results <- lapply(files, fread)
results <- rbindlist(results, idcol = 'emb')

# rename embeddings
emb_names <- c("emogi_STRING_EXP_test_set.tsv" = "STRING Exp.",
               "emogi_STRING_test_set.tsv" = "STRING", 
               "emogi_test_set.tsv" = "Omics")

results[, emb := emb_names[emb]]

results[method == "emb_only", Method := emb]
#results[method == "and_emb" & emb == "STRING", Method := paste0(emb, ' + EMOGI')]
results[method == "score_only" & emb == "Omics", Method := "EMOGI"]


#embeddings_dt <- prepare_pr_tables(results[!is.na(group)], "Method")
embeddings_dt<- results[!is.na(Method), calculatePrecisionRecall(pred, label), by = "Method"]

color_dict <- c(
  "STRING" = "#1f78b4",
  "STRING Exp." = "#a6cee3",
  "Omics" = '#ff7f00',
  "EMOGI" = '#b15928'
)

cancer_plot <- ggplot(embeddings_dt, aes(recall, precision, color = Method)) + geom_step(direction = 'vh') + 
  theme_cowplot(font_size = 14) + 
  xlab("Recall") + 
  ylab("Precision") + 
  scale_color_manual(values = color_dict) + 
  theme(
    legend.position = c(0.2, 0.2)
    #legend.justification = c("left", "top")
  )


saveRDS(cancer_plot, file = file.path(DATADIR,"figures/cancer/cancer_plot.rds"))

## print auPRC
unique(embeddings_dt[Method %in% c("STRING", "EMOGI"), .(Method, round(au_prc, 2))])

## compute Wilcoxon test
embeddings_dt[Method %in% c("STRING", "EMOGI")] %>% wilcox.test(score ~ Method, data = .)



# embeddings <- ggplot(embeddings_dt, aes(recall, precision, color = group)) + geom_line() + 
#   theme_light() + 
#   scale_color_manual(
#     name = "", 
#     labels = c(
#       glue("STRING:\t\t\t\t\t\tauPRC:  ", format(embeddings_dt[group=="combined_STRING", unique(au_prc)], digits = 3)), 
#       glue("STRING Experimental:\t\t\tauPRC:  ", format(embeddings_dt[group=="combined_STRING_EXP", unique(au_prc)], digits = 3)),
#       glue("Ours:\t\t\t\t\t\tauPRC:  ", format(embeddings_dt[group=="dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding", unique(au_prc)], digits = 3)),
#       glue("PoPS Experimental:\t\t\tauPRC:  ", format(embeddings_dt[group=="pops_mat_exp", unique(au_prc)], digits = 3))
#     ),
#     values = article_colors
#   ) + 
#   xlab("Recall") + 
#   ylab("Precision") + 
#   theme(
#     legend.position = c(0.55, 0.85)
#   )
# 
# methods_colors <- c("#91ED91", "#FFA500", "#6495ED")
# 
# methods <- ggplot(methods_dt, aes(recall, precision, color = group)) + geom_line() + 
#   theme_light() + 
#   scale_color_manual(
#     name = "",
#     labels = c(
#       glue("Emogi + Best Embedding (and_emb):\t\tauPRC:  ", format(methods_dt[group=="and_emb", unique(au_prc)], digits = 3)), 
#       glue("Embedding only (emb_only):\t\t\tauPRC:  ", format(methods_dt[group=="emb_only", unique(au_prc)], digits = 3)), 
#       glue("Emogi only (score_only):\t\t\t\tauPRC:  ", format(methods_dt[group=="score_only", unique(au_prc)], digits = 3))
#       
#     ),
#     values = methods_colors,
#   ) + 
#   xlab("Recall") + 
#   ylab("Precision") + 
#   theme(
#     legend.position = c(0.55, 0.85)
#   )
# 
# embeddings + methods
