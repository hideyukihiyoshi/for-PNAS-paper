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
summ.RCA <- summ.RCA[ls.ref,]
summ.RCA <- summ.RCA[order(rownames(summ.RCA)),]
write.csv(summ.RCA, "Result/Figure/table_vitro_RCA_reference.csv")

sessionInfo()
rm(list = ls())
