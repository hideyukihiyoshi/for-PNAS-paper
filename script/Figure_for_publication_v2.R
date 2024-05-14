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

##########################################################################
##### vivo 1. Scatter plot of human cells and mouse cells (slide 11) #####
##########################################################################
obj.5.Hs <- CreateSeuratObject(Read10X("Data/cellranger/GRCh38-and-mm10_SI-GA-E5/outs/filtered_gene_bc_matrices/GRCh38/"), project = "5. 3D iPIC (X+) implantation 2 month")
obj.5.Mm <- CreateSeuratObject(Read10X("Data/cellranger/GRCh38-and-mm10_SI-GA-E5/outs/filtered_gene_bc_matrices/mm10/"))
obj.6.Hs <- CreateSeuratObject(Read10X("Data/cellranger/GRCH38-and-mm10_SI-GA-F5/outs/filtered_gene_bc_matrices/GRCh38/"), project = "6. 3D iPIC (X+) implantation 6 month")
obj.6.Mm <- CreateSeuratObject(Read10X("Data/cellranger/GRCH38-and-mm10_SI-GA-F5/outs/filtered_gene_bc_matrices/mm10/"))
for(i in c(5, 6)){
  eval(parse(text = paste("tmp.Hs <- obj.", i, ".Hs", sep = "")))
  eval(parse(text = paste("tmp.Mm <- obj.", i, ".Mm", sep = "")))
  df.UMI <- data.frame(Barcode = union(rownames(tmp.Hs@meta.data), rownames(tmp.Mm@meta.data)))
  df.UMI$UMI_Hs <- ifelse(df.UMI$Barcode %in% rownames(tmp.Hs@meta.data), df.UMI$Barcode, 0)
  df.UMI$UMI_Mm <- ifelse(df.UMI$Barcode %in% rownames(tmp.Mm@meta.data), df.UMI$Barcode, 0)
  for(j in c(1:nrow(df.UMI))){
    if(df.UMI$UMI_Hs[j] != 0){
      df.UMI$UMI_Hs[j] <- as.numeric(tmp.Hs@meta.data[df.UMI$UMI_Hs[j], "nUMI"])
    }
    if(df.UMI$UMI_Mm[j] != 0){
      df.UMI$UMI_Mm[j] <- as.numeric(tmp.Mm@meta.data[df.UMI$UMI_Mm[j], "nUMI"])
    }
  }
  df.UMI$UMI_Hs <- as.numeric(df.UMI$UMI_Hs); df.UMI$UMI_Mm <- as.numeric(df.UMI$UMI_Mm)
  df.UMI$Colour <- ifelse(df.UMI$UMI_Hs > 0 & df.UMI$UMI_Mm > 0, paste("Multiplet: ", nrow(subset(df.UMI, UMI_Hs > 0 & UMI_Mm > 0)), " cells", sep = ""), 
                          ifelse(df.UMI$UMI_Hs > 0, paste("Human: ", nrow(subset(df.UMI, UMI_Hs > 0 & UMI_Mm == 0)), " cells", sep = ""), paste("Mouse: ", nrow(subset(df.UMI, UMI_Hs == 0 & UMI_Mm > 0)), " cells", sep = "")))
  max.UMI <- max(df.UMI$UMI_Hs, df.UMI$UMI_Mm, na.rm = T)
  tmp.plot <- ggplot(df.UMI, aes(x = UMI_Hs, y = UMI_Mm, colour = Colour)) + ggtitle(paste("Distribution of UMI counts mapped to reference\n", tmp.Hs@project.name, sep = "")) + xlab("UMI counts mapped to human reference (GRCh38)") + xlim(0, max.UMI) + ylab("UMI counts mapped to mouse reference (mm10)") + ylim(0, max.UMI) + 
    geom_abline(intercept = 0, slope = 1, lty = 3) + geom_point() + guides(colour = guide_legend(override.aes = list(size = 3))) + theme_bw() + theme(axis.text.y = element_text(angle = 90, hjust = 0.5), plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.position = "top", panel.grid = element_blank())
  png(paste("Result/Figure/plot_vivo_UMI_mapped_GRCh38_and_mm10_", tmp.Hs@project.name, ".png", sep = ""), width = 1800, height = 1800, res = 300)
  plot(tmp.plot); dev.off(); rm(tmp.plot)
  
  eval(parse(text = paste("ls.exclude.", i, " <- df.UMI$Barcode[df.UMI$UMI_Mm > 0]", sep = "")))
  
  rm(tmp.Hs, tmp.Mm, df.UMI, max.UMI)
}
rm(obj.5.Hs, obj.5.Mm, obj.6.Hs, obj.6.Mm)

###################################################################
##### vivo 2. Scatter plot of nUMI and percent.mito (slide 1) #####
###################################################################
obj.3 <- CreateSeuratObject(Read10X("Data/cellranger/SI-GA-C1/outs/filtered_gene_bc_matrices/GRCh38/"), project = "3. 3D compound X(+)")
obj.4 <- CreateSeuratObject(Read10X("Data/cellranger/SI-GA-C2/outs/filtered_gene_bc_matrices/GRCh38/"), project = "4. Human islet")
obj.5 <- CreateSeuratObject(Read10X("Data/cellranger/SI-GA-E5/outs/filtered_gene_bc_matrices/GRCh38/"), project = "5. 3D iPIC (X+) implantation 2 month")
obj.5 <- SubsetData(obj.5, cells.use = rownames(obj.5@meta.data)[!(rownames(obj.5@meta.data) %in% ls.exclude.5)])
obj.5@raw.data <- obj.5@raw.data[, (colnames(obj.5@raw.data) %in% rownames(obj.5@meta.data))]
obj.6 <- CreateSeuratObject(Read10X("Data/cellranger/SI-GA-F5/outs/filtered_gene_bc_matrices/GRCh38/"), project = "6. 3D iPIC (X+) implantation 6 month")
obj.6 <- SubsetData(obj.6, cells.use = rownames(obj.6@meta.data)[!(rownames(obj.6@meta.data) %in% ls.exclude.6)])
obj.6@raw.data <- obj.6@raw.data[, (colnames(obj.6@raw.data) %in% rownames(obj.6@meta.data))]
max.UMI <- max(c(obj.3@meta.data$nUMI, obj.4@meta.data$nUMI, obj.5@meta.data$nUMI, obj.6@meta.data$nUMI), na.rm = T)
for(i in c(3:6)){
  eval(parse(text = paste("tmp <- obj.", i, sep = "")))
  mito.genes <- grep("^MT-", rownames(tmp@raw.data), value = T)
  percent.mito <- colSums(tmp@raw.data[mito.genes,]) / colSums(tmp@raw.data)
  tmp <- AddMetaData(tmp, metadata = percent.mito, col.name = "percent.mito")
  tmp.plot <- ggplot(tmp@meta.data, aes(x = nUMI, y = percent.mito)) + xlim(0, max.UMI) + ylim(0, 1) + geom_point() + theme_bw() + theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1), panel.grid = element_blank())
  png(paste("Result/Figure/plot_vivo_nUMI_percent.mito_", tmp@project.name, ".png", sep = ""), width = 1800, height = 1800, res = 300)
  plot(tmp.plot); dev.off(); rm(tmp.plot)
  
  rm(tmp, mito.genes, percent.mito)
}
rm(obj.3, obj.4, obj.5, obj.6)

#########################################################
##### vivo 3. Table of summary statistics (slide 1) #####
#########################################################
df.SummStat.vitro <- read.csv("Result/summary_meta-data.csv", header = T)
df.SummStat.vitro <- df.SummStat.vitro[c("sampleID", "QC", "num_cells", "median.nUMI", "median.nGene")]
colnames(df.SummStat.vitro) <- c("Sample ID", "processing status", "# of cells", "median UMI counts per cell", "median genes per cell")
df.SummStat.vivo <- read.csv("Result/summary_mata-data.csv", header = T)
df.SummStat.vivo <- df.SummStat.vivo[c("sampleID", "QC", "number_of_cells", "nUMI.median", "nGene.median")]
colnames(df.SummStat.vivo) <- c("Sample ID", "processing status", "# of cells", "median UMI counts per cell", "median genes per cell")
df.SummStat <- rbind(df.SummStat.vitro, df.SummStat.vivo)
df.SummStat <- df.SummStat[c(3, 9, 13, 15, 17, 14, 16, 18, 4, 10),]
write.csv(df.SummStat, "Result/Figure/table_vivo_summary_statistics.csv", row.names = F)
rm(df.SummStat.vitro, df.SummStat.vivo, df.SummStat)

################################################################
##### vivo 3. Dot plot of marker gene expression (slide 2) #####
################################################################
obj.Merge <- readRDS("Rds/Seurat_object_merge_3to6.Rds")
png("Result/Figure/DotPlot_vivo_marker_genes_v1.png", width = 3600, height = 2400, res = 300)
DotPlot(obj.Merge, genes.plot = rev(c("INS", "GCG", "SST", "PPY", "PRSS1", 
                                      "KRT19", "GHRL", "ESAM", "CHGA", "PDX1", 
                                      "MKI67", "NKX6-1", "MAFA", "MAFB", "G6PC", 
                                      "UCN3", "IAPP", "FOXA1", "PCSK1", "PCSK2", 
                                      "GCK", "SLC2A3", "FFAR1", "ARX", "PTF1A", 
                                      "AMY1A", "CPA1", "HNF1B", "POU5F1", "SOX2", 
                                      "NANOG", "GATA4", "SOX9", "SPP1", "NEUROD1", 
                                      "NEUROG3", "YAP1", "FGFR1", "FGFR2", "FGFR3", 
                                      "FGFR4", "MUC1", "HES1", "LEPR", "PAX4", 
                                      "NKX2-2", "AFP", "ALB", "CTRB1", "CTRB2", 
                                      "VIM", "SPARC", "COL3A1", "COL1A1", "PDGFRB")), plot.legend = T, x.lab.rot = T, do.return = T) + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5), legend.position = "bottom", legend.title = element_blank()); dev.off()
png("Result/Figure/DotPlot_vivo_marker_genes_v2.png", width = 2400, height = 2400, res = 300)
DotPlot(obj.Merge, genes.plot = rev(c("FBLN1", "DCN", "SPARCL1", "HTRA3", "ELN", 
                                      "ASPN", "OGN", "FOXS1", "TGFB3", "PLAT", 
                                      "PLXDC1", "LRRC17", "SSPN", "SFRP4", "PLAC9", 
                                      "RGS5", "ACTA2", "IGF1", "MXRA8", "DPT")), plot.legend = T, x.lab.rot = T, do.return = T) + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 8)); dev.off()

################################################################
##### vivo 4. Dot plot of marker gene expression (slide 8) #####
################################################################
png("Result/Figure/DotPlot_vivo_marker_genes_v3.png", width = 1800, height = 1800, res = 300)
DotPlot(SubsetData(obj.Merge, ident.use = c(1, 7, 8)), genes.plot = rev(c("MAFA", "G6PC", "UCN3", "IAPP", "GCK", "SLC2A2", "PCSK1", "PCSK2")), plot.legend = T, x.lab.rot = T, do.return = T) + theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 8)); dev.off()

############################################################################
##### vivo 5. Dendrogram of average expression among cluster (slide 7) #####
############################################################################
png("Result/Figure/dendrogram_vivo_cluster.png", width = 1200, height = 1200, res = 300)
obj.Merge <- BuildClusterTree(obj.Merge, genes.use = rownames(obj.Merge@raw.data)); dev.off()
png("Result/Figure/dendrogram_vivo_cluster.png", width = 1200, height = 1200, res = 300)
ape::plot.phylo(obj.Merge@cluster.tree[[1]], direction = "downwards", show.node.label = F, no.margin = T); dev.off()
rm(obj.Merge)

############################################
##### vivo 6. Heatmap of RCA (slide 4) #####
############################################
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
png("Result/Figure/heatmap_vivo_RCA_green_module.png", width = 1200, height = 1800, res = 300)
ggplot(tmp.melt[order(tmp.melt$value),], aes(x = variable, y = factor(id, levels = tmp.melt$id[order(tmp.melt$value)]), fill = value)) + scale_fill_gradient(low = "white", high = "red") + geom_tile() + geom_text(mapping = aes(label = round(value, 2)), size = 2) + 
  ggtitle(paste("Heatmap of RCA\nmerged 4 samples", sep = "")) + labs(fill = "RCA score") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 6), axis.ticks = element_blank(), axis.line = element_blank(), legend.position = "none", panel.background = element_blank()); dev.off()
rm(out.RCA, ls.module, summ.RCA, tmp.melt)

#######################################################################################
##### vivo 7. Heatmap of gene expression among pseudotime (slide 5, 6, 9, and 10) #####
#######################################################################################
##### function #####
ggColorHue <- function(n, l=65) {
  hues <- seq(15, 375, length=n+1)
  hcl(h=hues, l=l, c=100)[1:n]
}
table.ramp <- function(n, mid = 0.5, sill = 0.5, base = 1, height = 1)
{
  x <- seq(0, 1, length.out = n)
  y <- rep(0, length(x))
  sill.min <- max(c(1, round((n - 1) * (mid - sill / 2)) + 1))
  sill.max <- min(c(n, round((n - 1) * (mid + sill / 2)) + 1))
  y[sill.min:sill.max] <- 1
  base.min <- round((n - 1) * (mid - base / 2)) + 1
  base.max <- round((n - 1) * (mid + base / 2)) + 1
  xi <- base.min:sill.min
  yi <- seq(0, 1, length.out = length(xi))
  i <- which(xi > 0 & xi <= n)
  y[xi[i]] <- yi[i]
  xi <- sill.max:base.max
  yi <- seq(1, 0, length.out = length(xi))
  i <- which(xi > 0 & xi <= n)
  y[xi[i]] <- yi[i]
  height * y
}
rgb.tables <- function(n,
                       red = c(0.75, 0.25, 1),
                       green = c(0.5, 0.25, 1),
                       blue = c(0.25, 0.25, 1))
{
  rr <- do.call("table.ramp", as.list(c(n, red)))
  gr <- do.call("table.ramp", as.list(c(n, green)))
  br <- do.call("table.ramp", as.list(c(n, blue)))
  rgb(rr, gr, br)
}
matlab.like <- function(n) rgb.tables(n)
matlab.like2 <- function(n)
  rgb.tables(n,
             red = c(0.8, 0.2, 1),
             green = c(0.5, 0.4, 0.8),
             blue = c(0.2, 0.2, 1))
blue2green2red <- matlab.like2
ls.cell <- c("aciner", "alpha", "beta", "duct", "enterochromaffin", "mesenchyme")
for(i in ls.cell){
  cds <- readRDS(paste("Rds/CellDataSet_merged_mouse_implantation_", i, ".Rds", sep = ""))
  out.differentialGeneTest <- readRDS(paste("Rds/out_differentialGeneTest_pseudotime_", i, ".Rds", sep = ""))
  out.differentialGeneTest <- subset(out.differentialGeneTest, qval < 0.05)
  cds <- cds[rownames(out.differentialGeneTest),]
  
  newdata <- data.frame(Pseudotime = seq(min(pData(cds)$Pseudotime), max(pData(cds)$Pseudotime),length.out = 100)) 
  m <- genSmoothCurves(cds, trend_formula = "~sm.ns(Pseudotime, df = 3)", relative_expr = T, new_data = newdata)
  m <- m[!apply(m, 1, sum) == 0,]
  m <- log10(m + 1)
  m <- m[!apply(m, 1, sd) == 0,]
  m <- t(scale(t(m), center = T))
  m <- m[is.na(rownames(m)) == F,]
  m[is.nan(m)] <- 0
  m[m > 3] <- 3
  m[m < -3] <- -3
  
  hm <- m
  rownames(hm) <- as.character(fData(cds)[rownames(hm), "gene_short_name"])
  dist.row <- as.dist((1 - cor(t(hm))) / 2)
  dist.row[is.na(dist.row)] <- 1
  bks <- seq(-3.1,3.1, by = 0.1)
  hmcols <- blue2green2red(length(bks) - 1)
  
  ph <- pheatmap(hm, useRaster = T, cluster_cols = F, cluster_rows = T, show_rownames = T, show_colnames = F, 
                 clustering_distance_rows = dist.row, clustering_method = "ward.D2", cutree_rows = 6, 
                 silent = T, filename = NA, breaks = bks, border_color = NA, color = hmcols)
  #ann.row <- data.frame(Cluster = factor(cutree(ph$tree_row, 6)))
  tmp.ann <- data.frame(gene = ph$gtable$grobs[[3]]$label, Order = c(1:length(ph$gtable$grobs[[3]]$label)), row.names = ph$gtable$grobs[[3]]$label)
  tmp.ann <- tmp.ann[rownames(hm),]
  tmp.ann$Cluster <- cutree(ph$tree_row, 6)
  tmp.ann$NewCluster <- rep(NA, nrow(tmp.ann))
  for(j in c(1:length(unique(tmp.ann[order(tmp.ann$Order),]$Cluster)))){
    tmp.ann$NewCluster <- ifelse(tmp.ann$Cluster == unique(tmp.ann[order(tmp.ann$Order),]$Cluster)[j], paste("gene cluster ", j, sep = ""), tmp.ann$NewCluster)
  }
  ann.row <- data.frame(GeneCluster = tmp.ann$NewCluster, row.names = tmp.ann$gene)
  #ann.col <- data.frame(Pseudotime = seq(1:100))
  ann.col <- data.frame(Pseudotime = seq(1:100))
  #ann.color <- list()
  tmp.color1 <- ggColorHue(6)
  names(tmp.color1) <- paste("gene cluster ", c(1:6), sep = "")
  tmp.color2 <- colorRampPalette(c("cyan", "magenta"))(100)
  names(tmp.color2) <- 1:100
  ann.color <- list("GeneCluster" = tmp.color1, "Pseudotime" = tmp.color2)
  
  ph2 <- pheatmap(hm[,], useRaster = T, cluster_cols = F, cluster_rows = T, show_rownames = F, show_colnames = F, 
                  clustering_distance_rows = dist.row, clustering_method = "ward.D2", cutree_rows = 6, 
                  annotation_row = ann.row, annotation_col = ann.col, annotation_colors = ann.color, annotation_names_col = F, annotation_names_row = F, 
                  treeheight_row = 20, breaks = bks, fontsize = 6, color = hmcols, border_color = NA, silent = T, filename = NA)
  ls.grobs <- ph2$gtable$grobs[[5]]$childrenOrder
  ls.grobs <- ls.grobs[grepl("Pseudotime", names(ls.grobs))]
  
  grid.rect(gp = gpar("fill", col = NA))
  png(paste("Result/Figure/pseudotime_heatmap_vivo_", i, ".png", sep = ""), width = 1800, height = 1200, res = 300)
  grid.draw(ph2$gtable)
  grid.ls(grid.force())
  for(j in ls.grobs){
    grid.gedit(j, gp = gpar(col = "white", fill = "white"))
  }
  dev.off()
  
  ann.row <- ann.row[order(rownames(ann.row)),]
  out.differentialGeneTest <- read.csv(paste("Result/out_differentialGeneTest_pseudotime_q0.05_", i, ".csv", sep = ""), header = T)
  rownames(out.differentialGeneTest) <- out.differentialGeneTest$gene_short_name
  out.differentialGeneTest <- out.differentialGeneTest[order(rownames(out.differentialGeneTest)),]
  ann.row <- cbind(ann.row, out.differentialGeneTest["cluster"])
  colnames(ann.row) <- c("GeneCluster", "old_cluster")
  write.csv(ann.row, paste("Result/Figure/table_pseudotime_cluster_vivo_", i," .csv", sep = ""), row.names = T)
  
  rm(cds, out.differentialGeneTest, newdata, m, hm, dist.row, bks, hmcols, ph, tmp.ann, ann.row, ann.col, tmp.color1, tmp.color2, ann.color, ph2, ls.grobs)
}

rm(list = ls())
