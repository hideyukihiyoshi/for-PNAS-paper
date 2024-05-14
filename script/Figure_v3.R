#############################
##### Setup Environment #####
#############################
rm(list = ls())
library("ggplot2"); library("gplots"); library("gridExtra"); library("grid")
library("pheatmap"); library("ggrepel"); library("RColorBrewer"); library("viridis")
library("reshape2"); library("dplyr"); library("rlist"); library("tidyr")
library("Seurat"); library("future")
library("scales")
set.seed(0)
options(stringsAsFactors = F, future.globals.maxSize = 10 * 1024^3)
#plan("multicore", workers = 6)
dir.out <- "Result/Figure/"
theme_set(theme_bw() + theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1), 
                             plot.title = element_text(hjust = 0.5), panel.grid = element_blank(), 
                             strip.background = element_blank()))
##### function #####
integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}
###################
##### Fig. 3D #####
###################
obj <- readRDS("Rds/Seurat_object_merge_3to6.Rds")
obj <- UpdateSeuratObject(obj)
ls.genes <- c("VIM", "COL3A1", "OGN", "IGF1", "SFRP4", "ACTA2")
tmp.plot <- FeaturePlot(obj, ls.genes, reduction = "tsne", combine = F)
tmp.plot <- lapply(tmp.plot, function(x){
  x + scale_colour_gradientn(colours = c("lightgrey", "blue"), breaks = integer_breaks()) + 
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1), 
          legend.position = "bottom", legend.justification = "center")
})
tmp.plot <- CombinePlots(tmp.plot, ncol = length(ls.genes))
png(paste0(dir.out, "Fig3d_t-SNE_plot_scale_bottom.png"), width = 1800 * length(ls.genes), height = 2100, res = 300)
plot(tmp.plot); dev.off(); rm(tmp.plot)

###################
##### Fig. 4G #####
###################
obj <- readRDS("Rds/Seurat_merged_7_samples_R3.6.1.Rds")
obj <- UpdateSeuratObject(obj)
ls.genes <- c("OGN", "COL3A1", "ACTA2")
M <- max(obj@assays$RNA@data[ls.genes,])
tmp.plot <- FeaturePlot(obj, ls.genes, reduction = "umap", combine = F)
tmp.plot <- lapply(tmp.plot, function(x){
  x + scale_colour_gradientn(colours = c("lightgrey", "blue"), breaks = integer_breaks()) + 
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1), 
          legend.position = "bottom", legend.justification = "center")
})
tmp.plot <- CombinePlots(tmp.plot, ncol = length(ls.genes))
png(paste0(dir.out, "Fig4g_UMAP_plot_scale_bottom.png"), width = 1800 * length(ls.genes), height = 2100, res = 300)
plot(tmp.plot); dev.off(); rm(tmp.plot)

###################
##### Fig.S6C #####
###################
ls.genes <- c("PDX1", "CHGA")
M <- max(obj@assays$RNA@data[ls.genes,])
tmp.plot <- FeaturePlot(obj, ls.genes, reduction = "umap", combine = F)
tmp.plot <- lapply(tmp.plot, function(x){
  x + scale_colour_gradientn(colours = c("lightgrey", "blue"), breaks = integer_breaks()) + 
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1), 
          legend.position = "bottom", legend.justification = "center")
})
tmp.plot <- CombinePlots(tmp.plot, ncol = length(ls.genes))
png(paste0(dir.out, "FigS6c_UMAP_plot_scale_bottom.png"), width = 1800 * length(ls.genes), height = 2100, res = 300)
plot(tmp.plot); dev.off(); rm(tmp.plot)

sessionInfo()
rm(list= ls())
