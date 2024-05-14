#############################
##### Setup Environment #####
#############################
rm(list = ls())
library("ggplot2"); library("gplots"); library("gridExtra"); library("grid")
library("pheatmap"); library("ggrepel"); library("RColorBrewer"); library("viridis")
library("reshape2"); library("dplyr"); library("rlist"); library("tidyr")
library("Seurat"); library("future")
set.seed(0)
options(stringsAsFactors = F, future.globals.maxSize = 10 * 1024^3)
#plan("multicore", workers = 6)
dir.out <- "Result/Figure/"
theme_set(theme_bw() + theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1), 
                             plot.title = element_text(hjust = 0.5), panel.grid = element_blank(), 
                             strip.background = element_blank()))

##########################
##### Supple Fig. 1D #####
##########################
obj <- readRDS("Rds/Seurat_object_merge_3to6.Rds")
obj <- UpdateSeuratObject(obj)
png("Result/Figure/Fig.S1D_DotPlot_vivo_iPIC_marker_genes.png", width = 3600, height = 2400, res = 300)
DotPlot(obj, features = c("INS", "GCG", "SST", "PPY", "PRSS1", 
                                      "KRT19", "GHRL", "ESAM", "CHGA", "PDX1", 
                                      "MKI67", "NKX6-1", "MAFA", "MAFB", "G6PC", 
                                      "UCN3", "IAPP", "FOXA1", "PCSK1", "PCSK2", 
                                      "GCK", "SLC2A3", "FFAR1", "ARX", "PTF1A", 
                                      "AMY1A", "CPA1", "HNF1B", "POU5F1", "SOX2", "NCAM1", "GFAP", 
                                      "NANOG", "GATA4", "SOX9", "SPP1", "NEUROD1", 
                                      "NEUROG3", "YAP1", "FGFR1", "FGFR2", "FGFR3", 
                                      "FGFR4", "MUC1", "HES1", "LEPR", "PAX4", 
                                      "NKX2-2", "AFP", "SLC18A1", "TPH1", "LMX1A", "ALB", "CTRB1", "CTRB2", 
                                      "VIM", "SPARC", "COL3A1", "COL1A1", "PDGFRB")) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 16, angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(size = 20), legend.position = "bottom", legend.title = element_blank()); dev.off()

##########################
##### Supple Fig. 1E #####
##########################
tmp.plot <- VlnPlot(obj, c("S.Score", "G2M.Score"), pt.size = 0, same.y.lims = T, combine = F)
tmp.plot <- lapply(tmp.plot, function(x){ x + theme(axis.title = element_blank(), 
                                                    axis.text = element_text(size = 18), 
                                                    axis.text.x = element_text(angle = 0, hjust = 0.5), 
                                                    plot.title = element_text(size = 30), 
                                                    legend.position = "none") })
png("Result/Figure/Fig.S1E_VlnPlot_vivo_iPIC_cell_cycle.png", width = 4000, height = 1800, res = 300)
plot(CombinePlots(tmp.plot, ncol = 2)); dev.off()

###################
##### Fig. 3B #####
###################
out.RCA <- readRDS("Rds/RCA_result_merge_3to6.Rds")
ls.module <- unique(out.RCA$group_labels_color$dynamicColors); ls.module <- ls.module[order(ls.module)]
summ.RCA <- data.frame(matrix(ncol = length(ls.module), nrow = nrow(out.RCA$fpkm_for_clust)))
colnames(summ.RCA) <- ls.module; rownames(summ.RCA) <- rownames(out.RCA$fpkm_for_clust)
ls.ref <- c()
for(i in ls.module){
  cellID <- rownames(subset(out.RCA$group_labels_color, dynamicColors == i))
  score <- apply(out.RCA$fpkm_for_clust[, cellID], 1, mean)
  summ.RCA[,i] <- score
  score <- score[order(score, decreasing = T)]
  ls.ref <- union(ls.ref, names(score)[1:10])
  rm(cellID, score)
}
tmp.melt <- summ.RCA[ls.ref,]; tmp.melt$id <- rownames(tmp.melt); tmp.melt <- melt(tmp.melt, id.vars = "id"); tmp.melt <- subset(tmp.melt, variable == "green" & value > 0)
png("Result/Figure/Fig.3B_vivo_RCA_green_module.png", width = 3600, height = 3600, res = 300)
ggplot(tmp.melt[order(tmp.melt$value),], aes(x = variable, y = factor(id, levels = tmp.melt$id[order(tmp.melt$value)]), fill = value)) + 
  scale_fill_gradient(low = "white", high = "red") + geom_tile() + geom_text(mapping = aes(label = round(value, 2)), size = 10) + 
  ggtitle(paste("Heatmap of RCA\nmerged 4 samples", sep = "")) + labs(fill = "RCA score") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 1), 
        axis.ticks = element_blank(), axis.line = element_blank(), plot.title = element_text(size = 30), legend.position = "right", legend.title = element_text(size = 20), legend.text = element_text(size = 20), panel.background = element_blank(), panel.border = element_blank()); dev.off()

###################
##### Fig. 3C #####
###################
png("Result/Figure/Fig.3C_DotPlot_vivo_iPIC_marker_genes.png", width = 2400, height = 2400, res = 300)
DotPlot(obj, features = c("FBLN1", "DCN", "SPARCL1", "HTRA3", "ELN", 
                          "ASPN", "OGN", "FOXS1", "TGFB3", "PLAT", 
                          "PLXDC1", "LRRC17", "SSPN", "SFRP4", "PLAC9", 
                          "RGS5", "ACTA2", "IGF1", "MXRA8", "DPT")) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 16, angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(size = 20), legend.position = "bottom", legend.title = element_blank()); dev.off()

###################
##### Fig. 3G #####
###################
obj <- BuildClusterTree(obj, features = rownames(obj))
png("Result/Figure/Fig.3G_dendrogram_vivo_iPIC.png", width = 1200, height = 1200, res = 300)
ape::plot.phylo(Tool(object = obj, slot = "BuildClusterTree"), direction = "downwards", 
                show.node.label = F, no.margin = T, cex = 1.5); dev.off()

##########################
##### Supple Fig. 2A #####
##########################
ls.module <- unique(out.RCA$group_labels_color$dynamicColors); ls.module <- ls.module[order(ls.module)]
summ.RCA <- data.frame(matrix(ncol = length(ls.module), nrow = nrow(out.RCA$fpkm_for_clust)))
colnames(summ.RCA) <- ls.module; rownames(summ.RCA) <- rownames(out.RCA$fpkm_for_clust)
ls.ref <- c()
for(i in ls.module){
  cellID <- rownames(subset(out.RCA$group_labels_color, dynamicColors == i))
  score <- apply(out.RCA$fpkm_for_clust[, cellID], 1, mean)
  summ.RCA[,i] <- score
  score <- score[order(score, decreasing = T)]
  ls.ref <- union(ls.ref, names(score)[1:10])
  rm(cellID, score)
}
tmp.melt <- summ.RCA[ls.ref,]; tmp.melt$id <- rownames(tmp.melt); tmp.melt <- melt(tmp.melt, id.vars = "id")
png("Result/Figure/Fig.S2A_heatmap_vivo_iPIC.png", width = 3600, height = 3600, res = 300)
ggplot(tmp.melt[order(tmp.melt$value),], aes(x = variable, y = id, fill = value)) + scale_fill_gradient(low = "white", high = "red") + geom_tile() + geom_text(mapping = aes(label = round(value, 2)), size = 2) + 
  ggtitle(paste("Heatmap of RCA\nmerged 4 samples", sep = "")) + labs(fill = "RCA score") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 12), 
        axis.ticks = element_blank(), axis.line = element_blank(), plot.title = element_text(size = 20), legend.title = element_text(size = 16), legend.text = element_text(size = 12), panel.background = element_blank(), panel.border = element_blank()); dev.off()
rm(out.RCA, ls.module, summ.RCA, tmp.melt)

###################
##### Fig. 4F #####
###################
out.RCA <- readRDS("Rds/RCA_merged_7_samples.Rds")
ls.module <- unique(obj.Merge@meta.data$RCA.merge); ls.module <- ls.module[order(ls.module)]
summ.RCA <- data.frame(matrix(ncol = length(ls.module), nrow = nrow(out.RCA$fpkm_for_clust)))
colnames(summ.RCA) <- ls.module; rownames(summ.RCA) <- rownames(out.RCA$fpkm_for_clust)
ls.ref <- c()
for(i in ls.module){
  cellID <- rownames(subset(out.RCA$group_labels_color, dynamicColors == i))
  ls.score <- apply(out.RCA$fpkm_for_clust[, cellID], 1, mean)
  summ.RCA[,i] <- ls.score
  ls.score <- ls.score[order(ls.score, decreasing = T)]
  ls.ref <- union(ls.ref, names(ls.score)[1:10])
  rm(cellID, ls.score)
}
tmp.melt <- summ.RCA[ls.ref,]
tmp.melt <- tmp.melt[((tmp.melt$pink > 0) & (tmp.melt$brown > 0)),]
tmp.melt <- tmp.melt[c("pink", "brown")]
tmp.melt <- tmp.melt[order(apply(tmp.melt, 1, mean), decreasing = T),]
ls.order <- rownames(tmp.melt)
tmp.melt$id <- rownames(tmp.melt); tmp.melt <- melt(tmp.melt, id.vars = "id")
png("Result/Figure/Fig.4F_heatmap_merged_7_samples_pink_brown_module.png", width = 3600, height = 3600, res = 300)
ggplot(tmp.melt, aes(x = variable, y = factor(id, levels = rev(ls.order)), fill = value)) + 
  ggtitle("Heatmap of RCA\nmerged 7 samples") + geom_tile() + scale_fill_gradient(low = "white", high = "red") + 
  geom_text(aes(label = round(value, 2)), size = 8) + labs(fill = "RCA score") + theme_bw() + 
  theme(axis.title = element_blank(), axis.text = element_text(size = 20), axis.ticks = element_blank(), axis.line = element_blank(), plot.title = element_text(hjust = 0.5, size = 25), legend.title = element_text(size = 20), legend.text = element_text(size = 18), panel.background = element_blank(), panel.grid = element_blank(), panel.border = element_blank()); dev.off()

##########################
##### Supple Fig. 6B #####
##########################
obj <- readRDS("Rds/Seurat_merged_7_samples_R3.6.1.Rds")
png("Result/Figure/Fig.S6B_DotPlot_merged_7_samples_marker_genes.png", width = 3600, height = 2400, res = 300)
DotPlot(obj, features = c("INS", "GCG", "SST", "PPY", "PRSS1", 
                          "KRT19", "GHRL", "ESAM", "CHGA", "PDX1", 
                          "MKI67", "NKX6-1", "MAFA", "MAFB", "G6PC", 
                          "UCN3", "IAPP", "FOXA1", "PCSK1", "PCSK2", 
                          "GCK", "SLC2A3", "FFAR1", "ARX", "PTF1A", 
                          "AMY1A", "CPA1", "HNF1B", "POU5F1", "SOX2", "NCAM1", "GFAP", 
                          "NANOG", "GATA4", "SOX9", "SPP1", "NEUROD1", 
                          "NEUROG3", "YAP1", "FGFR1", "FGFR2", "FGFR3", 
                          "FGFR4", "MUC1", "HES1", "LEPR", "PAX4", 
                          "NKX2-2", "AFP", "SLC18A1", "TPH1", "LMX1A", "ALB", "CTRB1", "CTRB2", 
                          "VIM", "SPARC", "COL3A1", "COL1A1", "PDGFRB")) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 16, angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(size = 20), legend.position = "bottom", legend.title = element_blank()); dev.off()

###################
##### Fig. 4D #####
###################
tmp.plot <- FeaturePlot(obj, c("PDX1", "CHGA"), combine = F)
tmp.plot <- lapply(tmp.plot, function(x){ x + theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1), 
                                                    legend.position = "none") })
png("Result/Figure/Fig.4D_FeaturePlot_merged_7_samples.png", width = 3600, height = 1800, res = 300)
plot(CombinePlots(tmp.plot, ncol = 2)); dev.off()
obj@meta.data$Sample <- ifelse(obj@meta.data$orig.ident %in% c("L1. 3D iPIC (X+) vitro long term culture -1", 
                                                               "L2. 3D iPIC (X+) vitro long term culture -2"), 
                               "long-term culture", "others")
Idents(obj) <- "Sample"
png("Result/Figure/Fig.4D_FeaturePlot_merged_7_samples_PDX1.png", width = 3600, height = 2000, res = 300)
FeaturePlot(obj, c("PDX1")) + facet_wrap("ident", ncol = 2) + theme_bw() + 
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1), 
        legend.position = "none", strip.text = element_text(size = 20), plot.title =element_text(size = 30, hjust = 0.5), 
        strip.background = element_blank(), panel.grid = element_blank()); dev.off()
png("Result/Figure/Fig.4D_FeaturePlot_merged_7_samples_CHGA.png", width = 3600, height = 2000, res = 300)
FeaturePlot(obj, c("CHGA")) + facet_wrap("ident", ncol = 2) + theme_bw() + 
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1), 
        legend.position = "none", strip.text = element_text(size = 20), plot.title =element_text(size = 30, hjust = 0.5), 
        strip.background = element_blank(), panel.grid = element_blank()); dev.off()

##########################
##### Supple Fig. 6C #####
##########################
out.RCA <- readRDS("Rds/RCA_merged_7_samples.Rds")
ls.module <- unique(obj@meta.data$RCA.merge); ls.module <- ls.module[order(ls.module)]
summ.RCA <- data.frame(matrix(ncol = length(ls.module), nrow = nrow(out.RCA$fpkm_for_clust)))
colnames(summ.RCA) <- ls.module; rownames(summ.RCA) <- rownames(out.RCA$fpkm_for_clust)
ls.ref <- c()
for(i in ls.module){
  cellID <- rownames(subset(out.RCA$group_labels_color, dynamicColors == i))
  ls.score <- apply(out.RCA$fpkm_for_clust[, cellID], 1, mean)
  summ.RCA[,i] <- ls.score
  ls.score <- ls.score[order(ls.score, decreasing = T)]
  ls.ref <- union(ls.ref, names(ls.score)[1:10])
  rm(cellID, ls.score)
}
tmp.melt <- summ.RCA[ls.ref,]
tmp.melt <- tmp.melt[order(apply(tmp.melt, 1, mean), decreasing = T),]
ls.order <- rownames(tmp.melt)
tmp.melt$id <- rownames(tmp.melt); tmp.melt <- melt(tmp.melt, id.vars = "id")
png(paste(dir.out, "tmp_heatmap.png", sep = ""))
tmp.heatmap <- heatmap(as.matrix(summ.RCA[ls.ref,]), col = bluered(256)); dev.off()
tmp.melt <- summ.RCA[ls.ref,]; tmp.melt <- tmp.melt[tmp.heatmap$rowInd, tmp.heatmap$colInd]
tmp.melt$id <- rownames(tmp.melt); tmp.melt <- melt(tmp.melt, id.vars = "id")
png("Result/Figure/Fig.6C_heatmap_merged_7_samples.png", width = 3600, height = 3600, res = 300)
ggplot(tmp.melt, aes(x = variable, y = id, fill = value)) + ggtitle("Heatmap of Reference Component Analysis\nmerged 7 samples") + geom_tile() + geom_text(aes(label = round(value, 2)), size = 2) + scale_fill_gradient(low = "white", high = "red") + labs(fill = "RCA score") + theme_minimal() + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 12), 
        axis.ticks = element_blank(), axis.line = element_blank(), plot.title = element_text(size = 20, hjust = 0.5), legend.title = element_text(size = 16), legend.text = element_text(size = 12), panel.background = element_blank(), panel.border = element_blank()); dev.off()

##########################
##### Supple Fig. 1A #####
##########################


sessionInfo()
rm(list = ls())
