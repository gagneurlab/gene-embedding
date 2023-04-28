

library(data.table)
library(magrittr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(yaml)

config <- read_yaml('config/config.yaml')
DATADIR <- config$datadir[1]

# files <- list.files('data/PoPS/magma_emb_pred', '*.tsv', recursive = T, full.names = TRUE)
# files <- list.files('data/PoPS/magma_emb_pred_cov_LD', '*.tsv', recursive = T, full.names = TRUE)
files <- list.files(file.path(DATADIR,'/gwas_prediction/'), '*.tsv', recursive = T, full.names = TRUE)

compute_r2 <- function(filename){
  dt <- fread(filename)
  dt[, cor(pred, ZSTAT)^2]
}

r2_values <- sapply(files, compute_r2)
dt_r2 <- data.table(filename = names(r2_values), r2 = r2_values)


dt_r2 <- separate(dt_r2, 'filename', into = c('A', "B", "C", "D","E", "F", "emb", "trait"), sep = "/")
dt_r2[, A := NULL]
dt_r2[, B := NULL]
dt_r2[, C := NULL]
         
dt_r2         
         
#dt_subset <- dt_r2[emb %in% c("norm-string_norm-stringverse", "norm-stringnolit_norm-stringnolitverse")]
dt_subset <- dcast(dt_r2, trait ~ emb)
dt_subset

ggplot(dt_subset, aes(`STRING_NO_LIT`, `STRING_NO_LIT.crispr_emb.gtex_emb.prot_t5_128_embedding`)) + 
  geom_point() + 
  geom_abline() + 
  geom_text_repel(aes(label = trait)) + 
  theme_cowplot()



dt_subset <- dt_r2[emb %in% c("norm-string_norm-stringverse", "norm-gtex_norm-crisprembnew_norm-stringnolit_norm-stringnolitverse_norm-proteinsmall")]
dt_subset <- dcast(dt_subset, trait ~ emb)

ggplot(dt_subset, aes(`norm-gtex_norm-crisprembnew_norm-stringnolit_norm-stringnolitverse_norm-proteinsmall`, `norm-string_norm-stringverse`)) + 
  geom_point() + 
  geom_abline() + 
  geom_text_repel(aes(label = trait))



ggplot(dt_subset, aes(trait, `norm-gtex_norm-crisprembnew_norm-stringnolit_norm-stringnolitverse_norm-proteinsmall`)) + 
  geom_bar(stat = 'identity') + labs(x = "", y = "R^2 using gene embedding") + coord_flip() + theme_cowplot()



dt_subset <- dt_r2[emb %in% c("norm-string_norm-stringverse", "pca-gtex_pca-crispr_norm-stringnolit_norm-stringnolitverse_norm-proteinsmall")]
dt_subset <- dcast(dt_subset, trait ~ emb)

ggplot(dt_subset, aes(`pca-gtex_pca-crispr_norm-stringnolit_norm-stringnolitverse_norm-proteinsmall`, `norm-string_norm-stringverse`)) + 
  geom_point() + 
  geom_abline() + 
  geom_text_repel(aes(label = trait))

# Compare to PoPS
# 
# pops <- fread('data/PoPS/paper_results/PoPS_FullResults.txt')
# pops
# dt_r2
# 
# 
# compute_pops_r2 <- function(filename){
#   dt <- fread(filename)
#   
#   dt[, cor(pred, ZSTAT)^2]
# }
# dt
# sel_trait <- strsplit(files[[5]], "/")[[1]][5]
# 
# pops[trait == sel_trait]

pops_files <- list.files("data/PoPS/pops_results", '*.results', recursive = TRUE, full.names = TRUE)

pops <- fread(pops_files[[1]])
pops

