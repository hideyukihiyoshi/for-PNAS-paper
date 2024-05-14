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

#####################################################
##### vitro 2. Trajectory plot among vitro iPIC #####
#####################################################
cds <- readRDS("Rds/CellDataSet_iPIC.Rds")
##### vitro 2-1. trajectory coloured by sample #####
png("Result/Figure/trajectory_vitro_iPIC_sampleID.png", width = 1800, height = 1200, res = 300)
plot_cell_trajectory(cds, color_by = "orig.ident", cell_size = 0.2) + 
  scale_color_manual(values = c("1. 2D compound X(+)" = "#F8766D", "2. 3D compound X(-)" = "#B79F00", "3. 3D compound X(+)" = "#00BA38", "a. 3D iPIC (X+)" = "#619CFF", "b. 3D iPIC (X+) kai" = "#F564E3")) + 
  guides(colour = guide_legend(override.aes = list(size = 3))) + 
  theme(legend.title = element_blank(), legend.position = "top", legend.justification = "left", legend.text = element_text(size = 6), legend.spacing.x = unit(0.2, "line")); dev.off()
##### vitro 2-2. trajectory coloured by state #####
png("Result/Figure/trajectory_vitro_iPIC_state.png", width = 1800, height = 1200, res = 300)
plot_cell_trajectory(cds, color_by = "State", cell_size = 0.2) + 
  guides(colour = guide_legend(override.aes = list(size = 3))) + 
  theme(legend.title = element_blank(), legend.position = "top", legend.justification = "left", legend.text = element_text(size = 6)); dev.off()
##### vitro 2-3. t-SNE plot coloured by state #####
obj.Seurat <- readRDS("Rds/Seurat_object_merge_1_2_3_4_a_b.Rds")
obj.Seurat <- SubsetData(SetAllIdent(obj.Seurat, "orig.ident"), ident.remove = "4. Human islet")
obj.Seurat@meta.data$Pseudotime <- pData(cds)$Pseudotime
obj.Seurat@meta.data$State <- pData(cds)$State
png("Result/Figure/t-SNE_vitro_iPIC_state.png", width = 1900, height = 1800, res = 300)
TSNEPlot(obj.Seurat, group.by = "State"); dev.off()
##### vitro 2-4. trajectory splited by cluster #####
pData(cds)$res.0.6 <- factor(pData(cds)$res.0.6, levels = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19"))
png("Result/Figure/trajectory_vitro_iPIC_cluster.png", width = 2000, height = 1500, res = 300)
plot_cell_trajectory(cds, color_by = "res.0.6", cell_size = 0.1) + 
  scale_color_manual(values = c("0" = "#F8766D", "1" = "#EA8331", "2" = "#D89000", "3" = "#C09B00", "4" = "#A3A500", 
                                "5" = "#7CAE00", "6" = "#39B600", "7" = "#00BB4E", "8" = "#00BF7D", "9" = "#00C1A3", 
                                "10" = "#00BFC4", "11" = "#00BAE0", "12" = "#00B0F6", "13" = "#35A2FF", "14" = "#9590FF", 
                                "15" = "#C77CFF", "16" = "#E76BF3", "17" = "#FA62DB", "18" = "#FF62BC", "19" = "#FF6A98")) + 
  facet_wrap("res.0.6", ncol = 5) + theme(legend.title = element_blank(), legend.position = "none"); dev.off()

rm(list = ls())
