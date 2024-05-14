#############################
##### Setup Environment #####
#############################
rm(list = ls())
library("ggplot2"); library("gplots"); library("gridExtra"); library("grid")
library("pheatmap"); library("ggrepel"); library("RColorBrewer"); library("viridis")
library("reshape2"); library("dplyr"); library("rlist")
library("Seurat")
set.seed(0)
options(stringsAsFactors = F)

###############################################################################
##### Visuzelize gene expression and cell-cycle within non-endcrine cells #####
###############################################################################
obj.Sub <- readRDS("Rds/Seurat_object_ex02.subclustering_CHGA (-).Rds")
##### Visualize marker gene expression #####
ls.genes <- c("CHGA", "PDX1", "MKI67", "NKX6-1", "NEUROD1", 
              "FGB", "AGR2", "FGFR1", "FGFR2", "PLK1", 
              "PLK4", "LDHB", "HNF1B", "SOX2", "GATA4", 
              "SOX9", "YAP1")
png("Result/Figure/t-SNE_plot_non-endocrine_cells_gene_expression_1of2.png", width = 9000, height = 3600, res = 300)
FeaturePlot(obj.Sub, ls.genes[1:10], cols.use = c("grey", "blue"), nCol = 5); dev.off()
png("Result/Figure/t-SNE_plot_non-endocrine_cells_gene_expression_2of2.png", width = 9000, height = 3600, res = 300)
FeaturePlot(obj.Sub, ls.genes[11:17], cols.use = c("grey", "blue"), nCol = 5); dev.off()
##### Visuzelize estimated cell-cycle score #####
png("Result/Figure/t-SNE_plot_non-endocrine_cells_cellcycle.png", width = 3600, height = 1800, res = 300)
FeaturePlot(obj.Sub, c("S.Score", "G2M.Score"), cols.use = c("grey", "red"), nCol = 2); dev.off()
png("Result/Figure/violin_plot_non-endocrine_cells_cellcycle.png", width = 3600, height = 1800, res = 300)
VlnPlot(obj.Sub, c("S.Score", "G2M.Score"), same.y.lims = T, size.x.use = 0, point.size.use = NA, nCol = 2); dev.off()

rm(list = ls())
