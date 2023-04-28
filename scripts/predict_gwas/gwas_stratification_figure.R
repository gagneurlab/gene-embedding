# Figure 5 
options(scipen = 10000)
library(timereg) # for qcut
library(data.table)
library(ggplot2)
library(ggpubr)
library(yaml)

config <- read_yaml('config/config.yaml')
DATADIR <- config$datadir[1]

#disease gene lists
d_gene_path <- sprintf('%s/input_data/disease_gene_eval/disease_gene_table.tsv', DATADIR)
d_gene <- fread(d_gene_path)
# d_gene[, table(disease)]

d_gene <- d_gene[disease %in% c("mito", "epilepsy", "IEI", "Neurology", "neuromuscular", "Ophthalmology")]

# genebass gene lists
gb_gene_path <- sprintf('%s/input_data/disease_gene_eval/disease_gene_table_genebass500k_selecttraits.tsv', DATADIR)
gb_gene <- fread(gb_gene_path)

# HPO genes
hpo_path <- sprintf('%s/processed_data/hpo/genes_to_phenotype_ensg.tsv', DATADIR)
hpo <- fread(hpo_path)

# read emb just to get set of genes
emb_gene_path <- sprintf('%s/embedding/combination/vtf_no_assay_var_pert3_dim_256_vtf_model_depMap.gtex.protT5_run_1_layers_3_embedding.tsv', DATADIR)
emb_genes <- fread(emb_gene_path)[, 'gene_id']

## combined lists
benchmark_genes <- list(
  unique(d_gene[, .(gene_id, benchmark = "Disease\nGene Lists")]),
  unique(hpo[, .(gene_id = Ensembl_ID, benchmark = "HPO")]),
  unique(gb_gene[, .(gene_id, benchmark = "Genebass")]),
  unique(emb_genes[, .(gene_id, benchmark = "Protein\nCoding")])
)

benchmark_genes <- rbindlist(benchmark_genes)
# & !is.na(cuts)

a <- ggplot(benchmark_genes[!is.na(benchmark)], aes(benchmark)) + geom_bar() + theme_cowplot() + 
  labs(x="", y = "Number of genes") + 
  guides(x = guide_axis(angle = 30)) + 
  theme(axis.title.x=element_blank())
a

## read publication counts 
pub_counts_path <- sprintf('%s/processed_data/pub_counts/pub_counts_ensg.tsv',DATADIR)
pub_count <- fread(pub_counts_path)

pub_count <- pub_count[Ensembl_gene_identifier %in% emb_genes[, gene_id]]
pub_count[,cuts := qcut(num_pub, cuts = 5, dig.lab = 5)]


## merge with publication counts
benchmark_genes <- merge(benchmark_genes, pub_count, by.x = 'gene_id', by.y = 'Ensembl_gene_identifier', all.x= T)

benchmark_genes$benchmark  <- factor(benchmark_genes$benchmark, level = c("Disease\nGene Lists","HPO","Genebass","Protein\nCoding"))

b <- ggplot(benchmark_genes[!is.na(benchmark) & !is.na(cuts)], aes(benchmark, fill = cuts)) + 
    geom_bar(position = 'fill')  + 
    labs(fill = "Publication\ncounts binned", x="", y = "Fraction of genes") + 
    theme_cowplot(font_size = 14) + 
    scale_fill_brewer() + 
    # guides(x = guide_axis(angle = 30)) + 
    theme(axis.title.x=element_blank(),
          legend.position = "right") 
    # guides(fill=guide_legend(nrow=2,byrow=TRUE))

b

# b <- ggpubr::ggboxplot(benchmark_genes, x = "benchmark", y = "num_pub + 1") + 
#   scale_y_log10() + background_grid() +
#   annotation_logticks(sides = 'l')
#   #stat_compare_means(comparisons = list(c(3, 4), c(2, 3), c(1,2), c(1,4)))

# b <- ggpubr::ggviolin(benchmark_genes, x = "benchmark", y = "num_pub + 1",
#                       add = "boxplot") + 
#   scale_y_log10() + background_grid() #+
# #stat_compare_means(comparisons = list(c(3, 4), c(2, 3), c(1,2), c(1,4)))


### GWAS Section ###########

# TRAITS = c("Alanine_aminotransferase","Cholesterol","IGF_1","Testosterone","Albumin","LDL_direct","Total_bilirubin","Alkaline_phosphatase","C_reactive_protein","LDL_direct_adjusted_by_medication","Total_protein","Apolipoprotein_A","Creatinine","Lipoprotein_A","Triglycerides","Apolipoprotein_B","Direct_bilirubin","Mean_corpuscular_haemoglobin","Urate","Glucose","Phosphate","Body_mass_index__BMI_","HDL_cholesterol","SHBG")
TRAITS <- config$SELECT_TRAITS$GWAS
names(TRAITS) <- TRAITS 

rsq_all_traits <- function(trait, embedding, string_conf){
  magma_file <- sprintf('%s/magma/%s/%s.genes.out',DATADIR, trait, trait)
  magma_pred <- fread(magma_file)
  
  magma_string <- merge(magma_pred, string_conf, by="GENE")
  
  pred_file <- sprintf("%s/gwas_prediction/%s_ElasticNet/%s/%s_pred_group_cv.tsv",DATADIR, embedding, trait, trait)
  pred_tab <- fread(pred_file)
  
  all_tab <- merge(magma_string, pred_tab, by='GENE', all.x=TRUE)
  all_tab[, pred := ifelse(is.na(pred), mean(pred, na.rm = TRUE), pred)]
  
  return(all_tab[, cor(ZSTAT.x, pred)^2])
}

gene_sets <- list()
for(set in pub_count[, unique(cuts)]){
  gene_sets[[set]] <- pub_count[cuts == set, .(GENE = Ensembl_gene_identifier)]    
}
gene_sets  
lapply(gene_sets, dim)


r2_table <- list()
for(embedding in list("pops_mat", "pops_mat_exp", "combined_STRING", "combined_STRING_EXP", "dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding", "ZEROS")){
  for(set in names(gene_sets)){
    r2_table_tmp <- melt(as.data.table(lapply(TRAITS, FUN = rsq_all_traits, embedding = embedding, string_conf = gene_sets[[set]])),
                         variable.name = 'Trait', value.name = "R2")
    r2_table_tmp[, embedding := embedding]
    r2_table_tmp[, set := set]
    r2_table[[paste0(embedding, "_", set)]] <- r2_table_tmp
  }
}


r2_table <- rbindlist(r2_table)

# this merges the zero R2 to all embeddings
r2_table <- merge(r2_table, r2_table[embedding == "ZEROS", .(Trait, set, R2_zero = R2)])
r2_table[, deltaR2 := R2 - R2_zero]

# Sort publication bins
r2_table[, set := factor(set, levels =  unique(set)[order(as.numeric(unlist(tstrsplit(sub(".", "", unique(set)), ",")[1])))])]
setnames(r2_table, 'set', 'Publication counts binned')
r2_table[, unique(set)]


lables <- c(
  "POPSMat" = "PoPS",
  "pops_mat_exp" = "PoPS Exp.",
  "combined_STRING" = "STRING",
  "combined_STRING_EXP" = "STRING Exp.",
  "dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding"  = "Omics",
  "ZEROS" = "ZEROS"
)

r2_table[, embedding := lables[embedding]]

color_vec <- c(
"STRING" = "#1f78b4",
"STRING Exp." = "#a6cee3",
"Omics" = '#ff7f00',
"PoPS" = "#33a02c",
"PoPS Exp." = '#b2df8a'
)


c <- ggboxplot(r2_table[embedding %in% c("Omics", "STRING", "STRING Exp.")], x = 'Publication counts binned', y = "deltaR2",
               color = "embedding", palette = "npg",
               #add = "jitter",
               #facet.by = "embedding", short.panel.labs = FALSE
               ) +
  # Use only p.format as label. Remove method name.
  #stat_compare_means(
    #aes(label = paste0("p = ", after_stat(p.format))), 
    #method = "wilcoxon", 
  #) + 
  # guides(x = guide_axis(angle = 30)) +
  theme_cowplot(font_size = 14) +
  scale_color_manual(values = color_vec) +
  labs(color='Embedding: ')  +
  background_grid() + 
  ylab(TeX("$\\Delta R^2$")) + 
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

c


### Combine all panels in one figure

ggarrange(b, c, nrow = 2, labels = c("A","B"), align = "v")
# ggarrange(b, c, ncol = 2, labels = c("A","B"), align = "h")
# plot_grid(plot_grid(a, b, labels = c("A", "B")), e, labels = c("", "C"), ncol = 1, rel_heights = c(1, 2))

