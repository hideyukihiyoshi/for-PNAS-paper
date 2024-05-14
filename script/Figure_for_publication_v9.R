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
obj.Seurat <- readRDS("Rds/Seurat_object_merge_1_2_3_4_a_b.Rds")
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
tmp.melt <- summ.RCA[ls.ref,]; tmp.melt <- tmp.melt[ls.ref[order(ls.ref)],]; tmp.melt$id <- rownames(tmp.melt); tmp.melt <- melt(tmp.melt, id.vars = "id")
png("Result/Figure/heatmap_vitro_RCA.png", width = 2700, height = 2700, res = 300)
ggplot(tmp.melt, aes(x = variable, y = id, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "red") + geom_text(aes(label = round(value,2)), size = 1.8) + theme_bw() + 
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "none", panel.border = element_blank(), panel.grid = element_blank()); dev.off()
rm(out.RCA, ls.module, summ.RCA, ls.ref, tmp.melt)

#########################################################
##### Perform differential gene expression analysis #####
#########################################################
out.FindAllMarkers <- readRDS("Rds/out_FindAllMarkers_merge_1_2_3_4_a_b_cluster.Rds")
out.FindMarkers.5_10_13 <- readRDS("Rds/out_FindMarkers_merge_1_2_3_4_a_b_cluster_5_10_13.Rds")
out.FindMarkers.5_10_13_16 <- readRDS("Rds/out_FindMarkers_merge_1_2_3_4_a_b_cluster_5_10_13_16.Rds")
ls.gene <- c("FGFR1", "FGFR2", "FGFR3", "FGFR4", "PLK1", "PLK2", "PLK3", "PLK4", "AURKA", "AURKB", "LDHA", "LDHB")
df.DEG <- data.frame(matrix(ncol = 24, nrow = 0), check.names = F)
for(i in c(5, 10, 13, 16)){
  tmp.FindAllMarkers <- subset(out.FindAllMarkers, cluster == i)
  ls.record <- c()
  for(j in ls.gene){
    ls.record <- c(ls.record, subset(tmp.FindAllMarkers, gene == j)$p_val, subset(tmp.FindAllMarkers, gene == j)$avg_diff)
  }
  df.DEG <- rbind(df.DEG, ls.record)
  rm(tmp.FindAllMarkers, ls.record)
}
ls.record1 <- c()
ls.record2 <- c()
for(i in ls.gene){
  ls.record1 <- c(ls.record1, out.FindMarkers.5_10_13[i,"p_val"], out.FindMarkers.5_10_13[i,"avg_diff"])
  ls.record2 <- c(ls.record2, out.FindMarkers.5_10_13_16[i,"p_val"], out.FindMarkers.5_10_13_16[i,"avg_diff"])
}
df.DEG <- rbind(df.DEG, ls.record1, ls.record2)
rownames(df.DEG) <- c("5", "10", "13", "16", "5 + 10 + 13", "5 + 10 + 13 + 16")
ls.col <- c()
for(i in 1:24){
  ls.col <- c(ls.col, ls.gene[ceiling(i/2)])
}
ls.col <- paste(ls.col, rep(c("p-value", "logFC"), 12), sep = "_")
colnames(df.DEG) <- ls.col
write.csv(df.DEG, "Result/Figure/table_differential_gene_expression_analysis_vitro_iPIC_cluster.csv", row.names = T)

sessionInfo()
rm(list = ls())
