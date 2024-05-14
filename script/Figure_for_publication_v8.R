#############################
##### Setup Environment #####
#############################
rm(list = ls())
library("ggplot2"); library("gplots"); library("gridExtra"); library("grid")
library("pheatmap"); library("ggrepel"); library("RColorBrewer"); library("viridis")
library("reshape2"); library("dplyr"); library("rlist")
library("Seurat"); library("monocle")
set.seed(0)
options(stringsAsFactors = F)

obj.Seurat <- readRDS("Rds/Seurat_object_merge_1_2_3_4_a_b.Rds")
png("Result/Figure/t-SNE_plot_vitro_iPIC_cellcycle.png", width = 4000, height = 1800, res = 300)
FeaturePlot(obj.Seurat, c("S.Score", "G2M.Score"), cols.use = c("grey", "red"), nCol = 2, no.legend = F); dev.off()

sessionInfo()
rm(list = ls())
