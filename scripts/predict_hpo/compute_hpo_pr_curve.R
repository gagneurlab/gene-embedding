library(data.table)
library(magrittr)
library(yaml)

config <- read_yaml('config/config.yaml')
DATADIR <- config$datadir[1]

calculatePrecisionRecall <- function(score, label, decr=TRUE, total=sum(label)){
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

# hpo_pred <- fread(file.path(DATADIR, '/processed_data/hpo/term_pred_DepMap18Q3_crispr_emb.gtex_emb.tsv'))
hpo_pred <- fread(snakemake@input[['hpo_pred']])


# this only contains the top predictions per term/gene hence we have to fill it with NA's.
hpo_filler <- fread(file.path(DATADIR, '/input_data/hpo/benchmark/prediction_HPOFiller_temporalpa_20190212_20200608.txt'), 
                    col.names = c('HPO_Filler', 'UniProtID', 'hpo_term'))

labels <- fread(file.path(DATADIR, '/processed_data/hpo/genes_to_phenotype_ensg.tsv'))
labels <- labels[, .(gene_id = Ensembl_ID, hpo_term = HPO_Term_ID, label = TRUE)]
## TODO first merge lables on protein coding genes.


## TODO fix this in the HPO_prediction script
if(colnames(hpo_pred)[1] == "V1"){
  colnames(hpo_pred)[1] <- "gene_id"  
}


hpo_pred <- melt(hpo_pred, 
                 id.vars = 'gene_id', 
                 value.name = 'emb_prediction',
                 variable.name = 'hpo_term'
)

mapping <- fread('/data/ouga/home/ag_gagneur/brechtma/Downloads/mart_export_UNIProt_mapping.txt', sep = '\t')
mapping <- mapping[, .(gene_id = `Gene stable ID`, UniProtID = `UniProtKB/Swiss-Prot ID`)]
mapping <- unique(mapping[UniProtID != ""])
mapping

hpo_filler <- merge(hpo_filler, mapping, by = 'UniProtID', allow.cartesian = T)

# subset to HPO terms predicted by both tools.
hpo_filler <- hpo_filler[hpo_term %in% intersect(hpo_pred[, hpo_term], hpo_filler[, hpo_term])]
hpo_pred <- hpo_pred[hpo_term %in% intersect(hpo_pred[, hpo_term], hpo_filler[, hpo_term])]


# subset HPOFiller predictions to the ones in test.
hpo_filler <- hpo_filler[gene_id %in% hpo_pred[, unique(gene_id)]]

hpo_pred[, uniqueN(gene_id)]
hpo_filler[, uniqueN(gene_id)]

hpo_pred[, uniqueN(hpo_term)]
hpo_filler[, uniqueN(hpo_term)]


hpo_predictions <- merge(hpo_filler, hpo_pred, by = c('gene_id', 'hpo_term'), all = TRUE)

hpo_predictions <- merge(hpo_predictions, labels, all.x = TRUE)
hpo_predictions[is.na(label), label := FALSE]


res <- melt(hpo_predictions, measure.vars = c('HPO_Filler', 'emb_prediction'), value.name = 'pred', variable.name = "Method")

#res[Method == 'prediction', pred := - pred]
res[Method == 'rank', pred := - pred]

# global_pr_curve <- res[, calculatePrecisionRecall(pred, label), by = Method]
# global_pr_curve

## By HPO term
hpo_term_pr_curves <- res[, calculatePrecisionRecall(pred, label), by = c("Method", 'hpo_term')]
print(snakemake@output['pr_curves'])
fwrite(hpo_term_pr_curves, file = snakemake@output[['pr_curves']], sep = '\t')



