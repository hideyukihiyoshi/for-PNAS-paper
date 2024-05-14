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

####################################
##### t-SNE plot of vitro iPIC #####
####################################
obj.Merge <- readRDS("/Rds/Seurat_object_merge_1_2_3_4_a_b.Rds")
png("Result/Figure/t-SNE_plot_vitro_iPIC_resize.png", width = 1800, height = 1800, res = 300)
TSNEPlot(obj.Merge, do.label = T, label.size = 5.6) + theme(legend.position = "none"); dev.off()

sessionInfo()
rm(list = ls())
