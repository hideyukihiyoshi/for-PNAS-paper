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

###############################################
##### Visualize RCA result in t-SNE space #####
###############################################
obj.Merge <- readRDS("Rds/Seurat_object_merge_1_2_3_4_a_b.Rds")
out.RCA <- readRDS("Rds/RCA_result_merge_1_2_3_4_a_b.Rds")
obj.Merge@meta.data$L7_Neuroepithelial <- t(out.RCA$fpkm_for_clust["L7_Neuroepithelial",])

tmp.plot <- FeaturePlot(obj.Merge, c("L7_Neuroepithelial"), cols.use = c("lightgrey", "blue"), no.legend = F, do.return = T); dev.off()
tmp.plot <- lapply(tmp.plot, function(x){ x + labs(colour = element_blank()) + theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1), legend.title = element_blank()) })
png("Result/Figure/t-SNE_vitro_iPIC_RCA_Neuroepithelial.png", width = 1800, height = 1800, res = 300)
tmp.plot$L7_Neuroepithelial; dev.off()
#grid.arrange(tmp.plot$Fetalliver, tmp.plot$Leukemialymphoblastic.MOLT.4, tmp.plot$colon, tmp.plot$Pancreas, tmp.plot$PancreaticIslet, ncol = 3); dev.off()

##################################################
##### Summarize DGEanalysis of LDHA and LDHB #####
##################################################
out.FindAllMarkers <- readRDS("Rds/out_FindAllMarkers_merge_1_2_3_4_a_b_cluster.Rds")
out.FindMarkers.c5_13 <- readRDS("Rds/out_FindMarkers_merge_1_2_3_4_a_b_cluster_5_13.Rds")
out.FindMarkers.c5_10_13 <- readRDS("Rds/out_FindMarkers_merge_1_2_3_4_a_b_cluster_5_10_13.Rds")
ls.gene <- c("FGFR1", "FGFR2", "MKI67", "TUBA1B", "PLK1", "PLK4", "LDHA", "LDHB")
summ.DGEanalysis <- data.frame(matrix(ncol = 16, nrow = 5))
line <- 1
for(i in c(5, 10, 13)){
  tmp.FindAllMarkers <- subset(out.FindAllMarkers, cluster == i)
  for(j in c(1:8)){
    summ.DGEanalysis[line, (2 * j - 1)] <- subset(tmp.FindAllMarkers, gene == ls.gene[j])$p_val
    summ.DGEanalysis[line, 2 * j] <- subset(tmp.FindAllMarkers, gene == ls.gene[j])$avg_diff
  }
  line <- line + 1
  rm(tmp.FindAllMarkers)
}
for(i in c(1:8)){
  summ.DGEanalysis[4, (2 * i - 1)] <- out.FindMarkers.c5_13[ls.gene[i], "p_val"]
  summ.DGEanalysis[4, 2 * i] <- out.FindMarkers.c5_13[ls.gene[i], "avg_diff"]
  summ.DGEanalysis[5, (2 * i - 1)] <- out.FindMarkers.c5_10_13[ls.gene[i], "p_val"]
  summ.DGEanalysis[5, 2 * i] <- out.FindMarkers.c5_10_13[ls.gene[i], "avg_diff"]
}
ls.colname <- c()
for(i in c(1:8)){
  ls.colname <- c(ls.colname, ls.gene[i], ls.gene[i])
}
ls.colname <- paste(ls.colname, rep(c("p-value", "logFC"), 8), sep = "_")
colnames(summ.DGEanalysis) <- ls.colname
rownames(summ.DGEanalysis) <- c("5", "10", "13", "5 + 10", "5 + 10 + 13")
write.csv(summ.DGEanalysis, "Result/Figure/table_differentaial_gene_expression_analysis_vitro_iPIC.csv", row.names = T)

rm(list = ls())
