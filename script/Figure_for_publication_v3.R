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

#########################################################
##### vivo Re-3. Dot plot of marker gene expression #####
#########################################################
obj.Merge <- readRDS("Rds/Seurat_object_merge_3to6.Rds")
png("Result/Figure/DotPlot_vivo_marker_genes_v1_modify.png", width = 3600, height = 2400, res = 300)
DotPlot(obj.Merge, genes.plot = rev(c("INS", "GCG", "SST", "PPY", "PRSS1", 
                                      "KRT19", "GHRL", "ESAM", "CHGA", "PDX1", 
                                      "MKI67", "NKX6-1", "MAFA", "MAFB", "G6PC2", 
                                      "UCN3", "IAPP", "FOXA1", "PCSK1", "PCSK2", 
                                      "GCK", "SLC2A2", "FFAR1", "ARX", "PTF1A", 
                                      "AMY1A", "CPA1", "HNF1B", "POU5F1", "SOX2", 
                                      "NANOG", "GATA4", "SOX9", "SPP1", "NEUROD1", 
                                      "NEUROG3", "YAP1", "FGFR1", "FGFR2", "FGFR3", 
                                      "FGFR4", "MUC1", "HES1", "LEPR", "PAX4", 
                                      "NKX2-2", "AFP", "ALB", "CTRB1", "CTRB2", 
                                      "VIM", "SPARC", "COL3A1", "COL1A1", "PDGFRB")), plot.legend = T, x.lab.rot = T, do.return = T) + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5), legend.position = "bottom", legend.title = element_blank()); dev.off()

#########################################################
##### vivo Re-4. Dot plot of marker gene expression #####
#########################################################
png("Result/Figure/DotPlot_vivo_marker_genes_v3_modify.png", width = 1800, height = 1800, res = 300)
DotPlot(SubsetData(obj.Merge, ident.use = c(1, 7, 8)), genes.plot = rev(c("MAFA", "G6PC2", "UCN3", "IAPP", "GCK", "SLC2A2", "PCSK1", "PCSK2")), plot.legend = T, x.lab.rot = T, do.return = T) + theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 8)); dev.off()

rm(list = ls())
