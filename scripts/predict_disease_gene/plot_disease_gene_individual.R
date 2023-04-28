
#snakemake <- readRDS('../../../gene-embedding/snakemake.rds')
saveRDS(snakemake, 'snakemake.rds')

library(ggplot2)
library(data.table)
library(dplyr)
library(cowplot)
library(ggpubr)
library(rstatix)


result_files <- snakemake@input[['predictions']]
print(length(result_files))
plot_path <- snakemake@output[['plot_path']]

emb_labels <- snakemake@params[['emb_labels']]
print(emb_labels)
color_dict <- snakemake@params[['color_dict']]

result_files <- lapply(result_files, fread)

print("Files read")

results <- rbindlist(result_files)
print("Rbind done")

results[, mean_average_precision := mean(average_precision), by = c('emb', 'disease')]
print(results)

# disease = c("Mito", "Ophthalmology", "Inborn Errors of Immunity", "Neurology", "Neuromuscular", "Epilepsy"),
# n_genes = c(425, 366, 364, 284, 132, 83),
disease_labels <- c(
  "EP" = "Epilepsy\n(N = 83)",
  "IEI" = "Inborn errors\nof immunity\n(N = 364)",
  "MITO" = "Mito\n(N = 400)",
  "NEU" = "Neurology\n(N = 284)",
  "NM" = "Neuromuscular\n(N = 131)",
  "OPH" = "Ophthalmology\n(N = 365)",
  "Exp" = "Exp",
  "Lit" = "Lit"
)

# order cols

results[, short_name := factor(short_name, levels = c(
                "MITO",#= "Mito\n(N = 425)",
                "OPH",# = "Ophthalmology\n(N = 366)"
                "IEI",# = "Inborn errors\nof immunity\n(N = 364)",
                "NEU",# = "Neurology\n(N = 284)",
                "NM",# = "Neuromuscular\n(N = 132)",
                "EP"# = "Epilepsy\n(N = 83)",
                )
)]


results[, emb_type := emb_type_dict[emb]]

results[, emb := factor(emb, levels = c(
    "combined_STRING",
    "combined_STRING_EXP",
    "pops_mat",                                                 
    "pops_exp_mat",
    "dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding"
)
)]

#emb_type_dict <- c("Lit", "Exp", "Exp", "Exp", "Lit", "Lit", "Exp")
emb_type_dict <- c("combined_STRING" = "Lit", 
                   "combined_STRING_EXP" = "Exp", 
                   "dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding" = "Exp", 
                    "pops_exp_mat" = "Exp", 
                    "pops_mat" = "Lit"
                  )

results[, emb_type := factor(emb_type, levels = c("Lit", "Exp"))]


print("Table ready")

stat.test <- results %>%
  group_by(emb_type, short_name) %>%
  wilcox_test(average_precision ~ emb) %>%
  #adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 
stat.test <- stat.test %>% filter(p.adj.signif != "ns")
stat.test <- stat.test %>% add_xy_position(x = "emb")

stat.test <- as.data.table(stat.test)
stat.test[, xmin:= results[, levels(emb)[xmin]]]
stat.test[, xmax:= results[, levels(emb)[xmax]]]



# p2 <- ggbarplot(results, x = "emb", y="average_precision", fill = "emb", 
#                 add = c("mean", "jitter"), facet = c("emb_type", "short_name")) + 
#   stat_pvalue_manual(stat.test, hide.ns = FALSE, p.adjust.method = "none") + 
#   theme(legend.pos = 'bottom', axis.ticks.x = element_blank(), axis.text.x = element_blank() ) +
#   xlab("") + 
#   facet_wrap(short_name ~ emb_type, scales = "free_x", nrow = 1, labeller = as_labeller(disease_labels)) + 
#   scale_fill_manual(labels = emb_labels, values = color_dict) + 
#   theme_cowplot(font_size = 10) + 
#   background_grid(major = 'y') +
#   theme(legend.pos = 'bottom', axis.ticks.x = element_blank(), axis.text.x = element_blank() ) + 
#   guides(fill = guide_legend(nrow = 1)) +
#   labs(fill='Embedding') +
#   ylab("auPRC") +
#   xlab(NULL)
# p2



p2 <- ggbarplot(results, x = "emb", y="average_precision", fill = "emb", 
                add = c("mean", "jitter"), facet = c("emb_type", "short_name")) + 
  stat_pvalue_manual(stat.test, hide.ns = FALSE, p.adjust.method = "none") + 
  theme(legend.pos = 'bottom', axis.ticks.x = element_blank(), axis.text.x = element_blank() ) +
  xlab("") + 
  facet_nested(~ short_name + emb_type, scales = "free_x", labeller = as_labeller(disease_labels)) + 
  scale_fill_manual(labels = emb_labels, values = color_dict) + 
  theme_cowplot(font_size = 10) + 
  background_grid(major = 'y') +
  theme(legend.pos = 'bottom', axis.ticks.x = element_blank(), axis.text.x = element_blank() ) + 
  guides(fill = guide_legend(nrow = 1)) +
  labs(fill='Embedding') +
  ylab("auPRC") +
  xlab(NULL)
p2





saveRDS(p2, paste0(plot_path, '_ggplot.RDS'))
