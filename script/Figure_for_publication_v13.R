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

##########################################################
##### RCA HeatMap of merged 6 saples (1to4, a and b) #####
##########################################################
out.RCA <- readRDS("Rds/RCA_result_merge_1_2_3_4_a_b.Rds")
ls.module <- unique(out.RCA$group_labels_color$dynamicColors); ls.module <- ls.module[order(ls.module)]
summ.RCA <- data.frame(matrix(ncol = length(ls.module), nrow = nrow(out.RCA$fpkm_for_clust)), row.names = rownames(out.RCA$fpkm_for_clust))
colnames(summ.RCA) <- ls.module
ls.ref <- c()
for(i in ls.module){
  ls.cell <- rownames(subset(out.RCA$group_labels_color, dynamicColors == i))
  tmp <- out.RCA$fpkm_for_clust[ls.cell]
  ls.score <- apply(tmp, 1, mean)
  summ.RCA[,i] <- ls.score
  ls.score <- ls.score[order(ls.score, decreasing = T)]
  ls.ref <- union(ls.ref, names(ls.score)[1:10])
  rm(ls.cell, tmp, ls.score)
}
tmp1 <- read.table("Data/temporary_RCAheatmap1_order_210519.tsv", header = F, sep = "\t")
tmp.melt <- summ.RCA[ls.ref,]; tmp.melt <- tmp.melt[ls.ref[order(ls.ref)],]; tmp.melt$id <- rownames(tmp.melt); tmp.melt <- melt(tmp.melt, id.vars = "id")
png("Result/Figure/heatmap_vitro_RCA_order.png", width = 2700, height = 2700, res = 300)
ggplot(tmp.melt, aes(x = variable, y = factor(id, rev(tmp1$V2)), fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(value,2)), size = 1.8) + theme_bw() + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "none", panel.border = element_blank(), panel.grid = element_blank()); dev.off()
tmp2 <- read.table("Data/temporary_RCAheatmap2_order_210519.tsv", header = T, sep = "\t", row.names = 1)
tmp.melt <- summ.RCA[rownames(tmp2),colnames(tmp2)]
tmp.melt$id <- rownames(tmp.melt); tmp.melt <- melt(tmp.melt, id.vars = "id")
png("Result/Figure/heatmap_vitro_RCA_filter_reference_coord_flip_order.png", width = 3000, height = 1200, res = 300)
ggplot(tmp.melt, aes(x = factor(variable, rev(colnames(tmp2))), y = factor(id, rownames(tmp2)), fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(value,2)), size = 3.6) + coord_flip() + theme_bw() + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "none", panel.border = element_blank(), panel.grid = element_blank()); dev.off()
png("Result/Figure/heatmap_vitro_RCA_filter_reference_order.png", width = 2700, height = 1800, res = 300)
ggplot(tmp.melt, aes(x = factor(variable, colnames(tmp2)), y = factor(id, rev(rownames(tmp2))), fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(value,2)), size = 3.6) + theme_bw() + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "none", panel.border = element_blank(), panel.grid = element_blank()); dev.off()
rm(summ.RCA, ls.module, ls.ref, tmp1, tmp2)

###################################
##### t-SNE plot of RCA score #####
###################################
obj.Seurat <- readRDS("Rds/Seurat_object_merge_1_2_3_4_a_b.Rds")
obj.tmp <- obj.Seurat
obj.tmp@meta.data$`Fetal Liver` <- t(out.RCA$fpkm_for_clust["Fetal Liver",])
obj.tmp@meta.data$`Colon` <- t(out.RCA$fpkm_for_clust["Colon",])

obj.Seurat@meta.data$RCA.module <- out.RCA$group_labels_color$dynamicColors
ls.module <- unique(out.RCA$group_labels_color$dynamicColors); ls.module <- ls.module[order(ls.module)]
png("Result/Figure/t-SNE_plot_vitro_iPIC_RCA_facet.png", width = 9000, height = 9000, res = 300)
TSNEPlot(obj.Seurat, group.by = "RCA.module", colors.use = ls.module) + facet_wrap("ident", ncol = 5) + 
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(size = 25)); dev.off()
png("Result/Figure/t-SNE_plot_vitro_iPIC_RCA_facet_no_color.png", width = 9000, height = 9000, res = 300)
TSNEPlot(obj.Seurat, group.by = "RCA.module", colors.use = rep("black", 21)) + facet_wrap("ident", ncol = 5) + 
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(size = 25)); dev.off()


sessionInfo()
rm(list = ls())
