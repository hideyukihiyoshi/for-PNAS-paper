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
obj.Merge <- readRDS("Result/Rds/Seurat_merged_7_samples.Rds")
out.RCA <- readRDS("Result/Rds/RCA_merged_7_samples.Rds")

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
png("Result/MergeDataSet/DotPlot/DotPlot_merged_7_samples_marker_genes_v2.png", width = 3600, height = 2400, res = 300)
DotPlot(obj.Merge, genes.plot = rev(ls.genes), plot.legend = T, x.lab.rot = T, do.return = T) + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5), legend.position = "bottom", legend.title = element_blank()); dev.off()

################################
##### Visualize RCA result #####
################################
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
png("Result/MergeDataSet/RCA/heatmap_merged_7_samples_pink_brown_module.png", width = 1200, height = 1800, res = 300)
ggplot(tmp.melt, aes(x = variable, y = factor(id, levels = rev(ls.order)), fill = value)) + ggtitle("Heatmap of RCA\nmerged 7 samples") + geom_tile() + scale_fill_gradient(low = "white", high = "red") + 
  geom_text(aes(label = round(value, 2)), size = 1.5) + theme_bw() + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none", panel.background = element_blank(), panel.grid = element_blank(), panel.border = element_blank()); dev.off()

rm(list = ls())
