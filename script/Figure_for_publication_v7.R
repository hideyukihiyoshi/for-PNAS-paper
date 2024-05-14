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

###########################################
##### Import processed scRNA-seq data #####
###########################################
obj.Merge <- readRDS("Rds/Seurat_merged_7_samples_R3.6.1.Rds")

#####################################
##### Visualize gene expression #####
#####################################
ls.genes <- c("INS", "GCG", "SST", "PPY", "PRSS1", 
              "KRT19", "GHRL", "ESAM", "CHGA", "PDX1", 
              "MKI67", "NKX6-1", "MAFA", "MAFB", "G6PC", 
              "UCN3", "IAPP", "FOXA1", "PCSK1", "PCSK2", 
              "GCK", "SLC2A3", "FFAR1", "ARX", "PTF1A", 
              "AMY1A", "CPA1", "HNF1B", "POU5F1", "SOX2", 
              "NANOG", "GATA4", "SOX9", "SPP1", "NEUROD1", 
              "NEUROG3", "YAP1", "FGFR1", "FGFR2", "FGFR3", 
              "FGFR4", "MUC1", "HES1", "LEPR", "PAX4", 
              "NKX2-2", "AFP", "ALB", "CTRB1", "CTRB2", 
              "VIM", "SPARC", "COL3A1", "COL1A1", "PDGFRB", 
              "ACTA2", "OGN", "IGF1", "SFRP4")
for(i in c(1:ceiling(length(ls.genes) / 10))){
  tmp.genes <- ls.genes[c((10 * i - 9):min(length(ls.genes), 10 * i))]
  tmp.width <- min(5, length(tmp.genes))
  tmp.height <- ceiling(length(tmp.genes) / 5)
  tmp.plot <- FeaturePlot(obj.Merge, features = tmp.genes, reduction = "umap", combine = F)
  tmp.plot2 <- tmp.plot
  tmp.plot <- lapply(tmp.plot, function(x){ x + theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1)) })
  tmp.plot <- CombinePlots(tmp.plot, ncol = tmp.width)
  png(paste("Result/MergeDataSet/FeaturePlot/T-CiRA00347_UMAP_plot_marker_genes_merged_7_samples_", i, "of", ceiling(length(ls.genes) / 10), ".png", sep = ""), width = 1800 * tmp.width, height = 1800 * tmp.height, res = 300)
  plot(tmp.plot); dev.off(); rm(tmp.plot)
  tmp.plot2 <- lapply(tmp.plot2, function(x){ x + scale_color_gradientn(colours = c("lightgrey", "blue"), limits = c(0, max(obj.Merge@assays$RNA@data))) + theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1), legend.position = "none") })
  tmp.plot2 <- CombinePlots(tmp.plot2, ncol = tmp.width)
  png(paste("Result/MergeDataSet/FeaturePlot/T-CiRA00347_UMAP_plot_marker_genes_merged_7_samples_", i, "of", ceiling(length(ls.genes) / 10), "_scale.png", sep = ""), width = 1800 * tmp.width, height = 1800 * tmp.height, res = 300)
  plot(tmp.plot2); dev.off(); rm(tmp.plot2)
  
  rm(tmp.genes, tmp.width, tmp.height)
}

rm(list = ls())
