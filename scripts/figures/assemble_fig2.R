## Figure 2

library(ggpubr)
library(png)
library(yaml)

config <- read_yaml('config/config.yaml')
DATADIR <- config$datadir[1]
# read panels as ggplot objects from RDS.

# panel A sketch.
panelA <- readPNG(file.path(DATADIR,'/figures/sketches/Figure1A.png'))
panelA <- as.raster(panelA)
panelA <- rasterGrob(panelA, interpolate = FALSE)



# panel B mito pr curve.
mito_pr_curve <- readRDS(file.path(DATADIR,'/figures/trait_gene_pred/individual_pr_panel/mito.RDS'))
mito_pr_curve <- mito_pr_curve + labs(title = "Mito")

# panel C auPRC comparison.

auPRC_panel <- readRDS(file.path(DATADIR,'/figures/trait_gene_pred/trait_gene_pred_individaul.png_ggplot.RDS'))

#auPRC_panel <- auPRC_panel + theme(legend.position = "top") + labs(fill="") + guides(fill = guide_legend(nrow = 3))


panel_ab <- ggarrange(panelA, mito_pr_curve, labels = c("A", "B"))
panel_c <- ggarrange(auPRC_panel, labels = c("C"))
full <- ggarrange(panel_ab, panel_c, nrow = 2, heights = c(2, 3))
full
