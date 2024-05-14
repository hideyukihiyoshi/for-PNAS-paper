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
##### t-SNE plot of merged 6 samples (1to4, a and b) #####
##########################################################
obj.Seurat <- readRDS("Rds/Seurat_object_merge_1_2_3_4_a_b.Rds")
out.RCA <- readRDS("Rds/RCA_result_merge_1_2_3_4_a_b.Rds")
obj.Seurat@meta.data$RCA.module <- out.RCA$group_labels_color$dynamicColors
ls.module <- unique(out.RCA$group_labels_color$dynamicColors); ls.module <- ls.module[order(ls.module)]
png("Result/Figure/t-SNE_plot_vitro_iPIC_RCA_facet.png", width = 9000, height = 9000, res = 300)
TSNEPlot(obj.Seurat, group.by = "RCA.module", colors.use = ls.module) + facet_wrap("ident", ncol = 5) + 
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(size = 25)); dev.off()
png("Result/Figure/t-SNE_plot_vitro_iPIC_RCA_facet_no_color.png", width = 9000, height = 9000, res = 300)
TSNEPlot(obj.Seurat, group.by = "RCA.module", colors.use = rep("black", 21)) + facet_wrap("ident", ncol = 5) + 
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(size = 25)); dev.off()

##########################################################
##### RCA HeatMap of merged 6 saples (1to4, a and b) #####
##########################################################
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
tmp.melt <- summ.RCA[ls.ref,]; tmp.melt[ls.ref[order(ls.ref)],]
tmp.melt <- tmp.melt[apply(tmp.melt[c("brown", "green", "magenta", "red")], 1, max) >= 2.5,]
tmp.melt$id <- rownames(tmp.melt); tmp.melt <- melt(tmp.melt, id.vars = "id")
png("Result/Figure/heatmap_vitro_RCA_filter_reference.png", width = 2700, height = 1800, res = 300)
ggplot(tmp.melt, aes(x = variable, y = id, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(value,2)), size = 1.8) + theme_bw() + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "none", panel.border = element_blank(), panel.grid = element_blank()); dev.off()

sessionInfo()
rm(list = ls())
