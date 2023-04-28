
## plot single Precision-Recall curve for Mito.
library(data.table)
library(glue)
library(yaml)

config <- read_yaml('config/config.yaml')
DATADIR <- config$datadir[1]

files <- sapply(
  c('combined_STRING',
    'combined_STRING_EXP', 
    "dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding",
    "pops_exp_mat",
    "pops_mat"
    ),
  function(x){
    paste0(DATADIR, '/Evaluation/hyperopt/predictions/mito/', x, '.tsv')
  })

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
    pr <- dt[, calculatePrecisionRecall(pred, target)]
    out <- cbind(pr, dt[, get(group_col)])
    setnames(out, "V2", "group")
    out
  })
  pr_dt <- rbindlist(pr_tables)
}




names(files) <- gsub(".tsv", "", basename(files))
tables <- lapply(files, fread)
embeddings_dt <- rbindlist(tables, fill = TRUE, idcol = "filepath")


embeddings_dt <- prepare_pr_tables(embeddings_dt, "filepath")

article_colors <- c('#a6cee3', '#1f78b4', '#ff7f00', '#b2df8a', '#33a02c')
names(article_colors) <- c('combined_STRING_EXP',
                           'combined_STRING',
                           "dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding",
                           "pops_exp_mat",
                           "pops_mat")

embeddings <- ggplot(embeddings_dt, aes(recall, precision, color = group)) + 
  #geom_step(direction = 'vh') + 
  geom_line() + 
  theme_cowplot(font_size = 10) + 
  scale_color_manual(
    name = "", 
    labels = c(
      glue("STRING\t\t\t\t\t\t\t\t\t\t\tauPRC: ", format(embeddings_dt[group=="combined_STRING", unique(au_prc)], digits = 2)), 
      glue("STRING Exp.\t\t\tauPRC: ", format(embeddings_dt[group=="combined_STRING_EXP", unique(au_prc)], digits = 2)),
      glue("Omics\t\t\t\t\t\t\t\t\t\t\t\t\t\tauPRC: ", format(embeddings_dt[group=="dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding", unique(au_prc)], digits = 2)),
      glue("PoPS Exp.\t\t\t\t\t\tauPRC: ", format(embeddings_dt[group=="pops_exp_mat", unique(au_prc)], digits = 2)),
      glue("PoPS\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tauPRC: ", format(embeddings_dt[group=="pops_mat", unique(au_prc)], digits = 2))
    ),
    values = article_colors
  ) + 
  background_grid() + 
  xlab("Recall") + 
  ylab("Precision") + 
  theme(
    #legend.position = "none",
    legend.position = c(1, 1),
    legend.justification = c("right", "top")
  )
embeddings


saveRDS(embeddings, file.path(DATADIR,'/figures/trait_gene_pred/individual_pr_panel/mito.RDS'))



