
library(data.table)
library(ggplot2)
library(cowplot)
library(rstatix)
library(yaml)

config <- read_yaml('config/config.yaml')
DATADIR <- config$datadir[1]

files <- c( file.path(DATADIR,"/processed_data/hpo/pr_curves/term_pred_combined_STRING.tsv"), file.path(DATADIR,"/processed_data/hpo/pr_curves/term_pred_combined_STRING_EXP.tsv"), file.path(DATADIR,"/processed_data/hpo/pr_curves/term_pred_dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding.tsv"))
names(files) <- basename(files)


files <- lapply(files, fread)
files

plottable <- rbindlist(files, idcol = TRUE)


bxplottable <- unique(plottable[, .(.id, Method, hpo_term, au_prc)])

lables <- c("term_pred_combined_STRING.tsv" = "STRING",
            "term_pred_combined_STRING_EXP.tsv" = "STRING\nExperimental",
            "term_pred_dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding.tsv"  = "Omics"
)

# subset to one 
bxplottable[Method == "emb_prediction"]
bxplottable[Method == "emb_prediction", lable := lables[.id]]

bxplottable[Method == "HPO_Filler" & .id == bxplottable[, .id][[1]], lable := "HPO_Filler"]

bxplottable <- rbind(bxplottable[Method == "emb_prediction"], bxplottable[Method == "HPO_Filler" & .id == bxplottable[, .id][[1]]])

# set order
bxplottable[, lable := factor(lable, levels = c("HPO_Filler", "Omics", "STRING\nExperimental", "STRING"))]

color_vec <- c('#a6cee3', '#1f78b4', '#ff7f00', '#b15928')
names(color_vec) <- c("STRING\nExperimental", "STRING", "Omics", "HPO_Filler")

ggplot(bxplottable, aes(lable, au_prc, color = lable)) + geom_boxplot()  +
  theme_cowplot() + 
  theme(legend.position = "none") + scale_y_log10()


hpo_plot <- ggplot(bxplottable, aes(lable, au_prc, color = lable)) + geom_boxplot()  +
  theme_cowplot(font_size = 14) + 
  scale_color_manual(values = color_vec) + 
  scale_y_log10() + 
  annotation_logticks(sides = "l") + 
  ylab("auPRC") + 
  background_grid(major = "y") + 
  theme(legend.position = "none") + 
  # guides(x = guide_axis(angle = 45)) + 
  scale_x_discrete(labels = c("HPOFiller" ,"Omics", "STRING Exp.", "STRING")) +
  xlab("") #+ 
#  stat_compare_means(comparisons = list(c(1,2), c(2,3), c(3,4), c(1,2), c(2, 3), c(1, 4)))

# hpo_plot
saveRDS(hpo_plot, file = file.path(DATADIR,"/figures/hpo/hpo_plot.rds"))


stat_pvalue <- bxplottable %>%
  wilcox_test(au_prc ~ lable, paired = TRUE) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat_pvalue
  # filter(p.adj.signif == 'ns') %>%
  # add_y_position() %>%
  # mutate(y.position = seq(min(y.position), max(y.position), length.out = n()))

# ggplot(bxplottable, aes(lable, au_prc, fill = lable)) + geom_boxplot()  +
#   theme_cowplot() + 
#   scale_fill_manual(values = color_vec) + 
#   scale_y_log10() + 
#   ylab("auPRC") + 
#   theme(legend.position = "none") + 
#   xlab("") 
# 
# 
# ggplot(bxplottable, aes(lable, au_prc, fill = lable, color = lable)) + geom_boxplot()  +
#   theme_cowplot() + 
#   scale_fill_manual(values = color_vec) + 
#   scale_y_log10() + 
#   ylab("auPRC") + 
#   theme(legend.position = "none") + 
#   xlab("") 
