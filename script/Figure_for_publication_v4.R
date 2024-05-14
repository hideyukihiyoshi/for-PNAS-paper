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
#obj.Merge <- UpdateSeuratObject(obj.Merge)
out.RCA <- readRDS("Rds/RCA_result_merge_1_2_3_4_a_b.Rds")
obj.Merge@meta.data$Fetalliver <- t(out.RCA$fpkm_for_clust["Fetalliver",])
obj.Merge@meta.data$Leukemialymphoblastic.MOLT.4 <- t(out.RCA$fpkm_for_clust["Leukemialymphoblastic.MOLT.4",])
obj.Merge@meta.data$colon <- t(out.RCA$fpkm_for_clust["colon",])
obj.Merge@meta.data$Pancreas <- t(out.RCA$fpkm_for_clust["Pancreas",])
obj.Merge@meta.data$PancreaticIslet <- t(out.RCA$fpkm_for_clust["PancreaticIslet",])

#tmp.plot <- FeaturePlot(obj.Merge, features = c("Fetalliver", "Leukemialymphoblastic.MOLT.4", "colon", "Pancreas", "PancreaticIslet"), reduction = "tsne", combine = F)
#tmp.plot <- lapply(tmp.plot, function(x){ x + theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1)) })
#tmp.plot <- CombinePlots(tmp.plot, ncol = 3)
#png("Result/Figure/t-SNE_vitro_iPIC_RCA.png", width = 5400, height = 3600, res = 300)
#plot(tmp.plot); dev.off()

tmp.plot <- FeaturePlot(obj.Merge, c("Fetalliver", "Leukemialymphoblastic.MOLT.4", "colon", "Pancreas", "PancreaticIslet"), cols.use = c("lightgrey", "blue"), no.legend = F, do.return = T); dev.off()
tmp.plot <- lapply(tmp.plot, function(x){ x + labs(colour = element_blank()) + theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1), legend.title = element_blank()) })
png("Result/Figure/t-SNE_vitro_iPIC_RCA.png", width = 5400, height = 3600, res = 300)
grid.arrange(tmp.plot$Fetalliver, tmp.plot$Leukemialymphoblastic.MOLT.4, tmp.plot$colon, tmp.plot$Pancreas, tmp.plot$PancreaticIslet, ncol = 3); dev.off()

rm(list = ls())
