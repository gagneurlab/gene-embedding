library(ggplot2)
library(ggpubr)
library(magick)
library(yaml)

config <- read_yaml('config/config.yaml')
DATADIR <- config$datadir[1]

# "Code/gene-embedding/general-purpose-gene-embeddings/scripts/predict_cancer_gene/HPO_model.png"
hpo_model <- ggdraw() + 
  draw_image(file.path(DATADIR,"/figures/sketches/HPO_model.png"), scale = .75)

hpo_plot <- readRDS(file.path(DATADIR,"/figures/hpo/hpo_plot.rds"))

cancer_model <- ggdraw() +
  draw_image(file.path(DATADIR,"/figures/sketches/Cancer_model.png"), scale = 0.9)

cancer_plot <- readRDS(file.path(DATADIR, "/figures/cancer/cancer_plot.rds"))

ggarrange(hpo_model, hpo_plot, cancer_model, cancer_plot, ncol = 2, nrow = 2, widths = c(2,3), labels = c("A","B","C","D"))


