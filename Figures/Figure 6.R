library(Seurat)
library(Signac)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(qlcMatrix)
library(reshape2)
library(scales)
library(Matrix)
library(WGCNA)
library(flashClust)
library(Hmisc)
library(qlcMatrix)
library(reshape2)
library(scales)
library(Matrix)
library(stringr)


#require "Preprocessing.R"
#require processing by "recall ATAC peaks and build motifs.R"
#require processing by "build WGCNA.R"
#require "Antibody_stats_HIVpos.csv"

#Fig. 6A
DimPlot(HIVposCells, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, raster = FALSE, pt.size = 2) 

#Fig. 6B
HIVposCells$HIV_RNA_DNA2 <- HIVposCells$HIV_RNA_DNA
Idents(HIVposCells) <- "HIV_RNA_DNA2"
HIVposCells <-RenameIdents(HIVposCells, 'positive_positive' = 'HIV RNA and double positive', 
                       'positive_negative' = 'HIV RNA and double positive', 'negative_positive' = 'HIV DNA positive',
                       'negative_negative' = 'HIV negative') 
HIVposCells$HIV_RNA_DNA2 <- Idents(HIVposCells)
HIVposCells$Condition_HIV_RNA_DNA2 <- str_c(HIVposCells$Condition, '_', HIVposCells$HIV_RNA_DNA2)

DimPlot(HIVposCells, reduction = "umap.rna", group.by = "Condition_HIV_RNA_DNA2", raster = FALSE, pt.size = 2) 

#Fig. 6C
HIVposCells$Participants <- HIVposCells$MULTI_ID
Idents(HIVposCells) <- "Participants"
HIVposCells <-RenameIdents(HIVposCells, 'hashtag1-TotalA' = '236', 'hashtag2-TotalA' = '640', 'hashtag3-TotalA' = '839', 'hashtag4-TotalA' = '739', 'hashtag5-TotalA' = '799', 'hashtag6-TotalA' = '910') 
HIVposCells$Participants <- Idents(HIVposCells)

DimPlot(HIVposCells, reduction = "umap.rna", group.by = "Participants", raster = FALSE, pt.size = 2) 

#Fig. 6D
Idents(HIVposCells) <- "seurat_clusters"; DefaultAssay(HIVposCells) <- "RNA"
heatmap.HIVposCells.markers <- FindAllMarkers(HIVposCells, test.use = "wilcox", only.pos = TRUE, min.pct = 0.15, logfc.threshold = 0.25, pseudocount.use = 1)
heatmap.HIVposCells.markers$BH <- p.adjust(heatmap.HIVposCells.markers$p_val, 
                                          method = "BH", n = nrow(heatmap.HIVposCells.markers))

heatmap.HIVposCells.markers %>%
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC) -> heatmap.HIVposCells.markers

DoHeatmap(HIVposCells, features = heatmap.HIVposCells.markers$gene, label = TRUE, raster = FALSE) + NoLegend()

#Fig. 6E
DefaultAssay(HIVposCells.ATAC) <- 'chromvar'; Idents(HIVposCells.ATAC) <- "seurat_clusters"
motif.names <- HIVposCells.ATAC@assays$ATAC@motifs@motif.names
HIVposCells.ATAC@assays$chromvar@counts <- HIVposCells.ATAC@assays$chromvar@data

heatmap.HIVpos.Chromvar.markers <- FindAllMarkers(HIVposCells.ATAC,  
                                                  only.pos = TRUE, min.pct = 0,
                                                  logfc.threshold = 0.3, mean.fxn = rowMeans, fc.name = "avg_diff", pseudocount.use = 1)
heatmap.HIVpos.Chromvar.markers$BH <- p.adjust(heatmap.HIVpos.Chromvar.markers$p_val, 
                                               method = "BH", n = nrow(heatmap.HIVpos.Chromvar.markers))
heatmap.HIVpos.Chromvar.markers$gene <- ConvertMotifID(motif.names, id = rownames(heatmap.HIVpos.Chromvar.markers))
heatmap.HIVpos.Chromvar.markers$motif <- rownames(heatmap.HIVpos.Chromvar.markers)
heatmap.HIVpos.Chromvar.markers <- heatmap.HIVpos.Chromvar.markers[!is.na(heatmap.HIVpos.Chromvar.markers$gene),]

heatmap.HIVpos.Chromvar.markers %>%
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_diff) -> top10.heatmap.HIVpos.Chromvar.markers

avg_exp.HIVpos.chromvar <- AverageExpression(HIVposCells.ATAC, return.seurat = TRUE,  group.by = 'seurat_clusters', assays = "chromvar", slot = "counts")
avg_exp.HIVpos.chromvar$orig.ident <- colnames(avg_exp.HIVpos.chromvar)
DefaultAssay(avg_exp.HIVpos.chromvar) <- "chromvar"; Idents(avg_exp.HIVpos.chromvar) <- "orig.ident"
Idents(avg_exp.HIVpos.chromvar) <- factor(Idents(avg_exp.HIVpos.chromvar), levels = c("IRF", "cytotoxic", "AP-1", "MT"))

DoHeatmap(object = avg_exp.HIVpos.chromvar, features = rev(top10.heatmap.HIVpos.Chromvar.markers$motif), slot = "counts",
          draw.lines = FALSE) + scale_y_discrete(labels = top10.heatmap.HIVpos.Chromvar.markers$gene) + scale_fill_gradientn(colors = c("blue", "white", "red"))

#Fig. 6F
Antibody_stats_HIVpos <- data.frame(antibody = HIVposCells[["Antibody"]]@data)
Antibody_stats_HIVpos$mean <- rowMeans(Antibody_stats_HIVpos, na.rm=TRUE)
Antibody_stats_HIVpos$sd <- apply(Antibody_stats_HIVpos, 1, sd, na.rm=TRUE)
Antibody_stats_HIVpos$names <- rownames(Antibody_stats_HIVpos)
Antibody_stats_HIVpos <- Antibody_stats_HIVpos[,c(523, 524, 525)]
#match  each row of protein rows with its specific isotype control (using e.g., excel; required columns are protein and standard deviations and mean expressions) then compute 2 sample z score

Antibody_stats_HIVpos$X2.sample.z.score <- (Antibody_stats_HIVpos$mean - 
                                              Antibody_stats_HIVpos$isotype_mean)/sqrt(Antibody_stats_HIVpos$sd^2/522 +
                                                                                         Antibody_stats_HIVpos$isotype_sd^2/522)
Antibody_stats_HIVpos <- Antibody_stats_HIVpos[Antibody_stats_HIVpos$X2.sample.z.score > 2,]
rownames(Antibody_stats_HIVpos) <- Antibody_stats_HIVpos$names
#see Antibody_stats_HIVpos.csv

DefaultAssay(HIVposCells) <- "CITE_dsb"; Idents(HIVposCells) <- "seurat_clusters"
heatmap.HIVposCells.Antibody.markers  <- FindAllMarkers(object = HIVposCells, min.pct = 0.15, 
                                                   only.pos = TRUE, test.use = 'wilcox', features =rownames(Antibody_stats_HIVpos))
heatmap.HIVposCells.Antibody.markers$BH <- p.adjust(heatmap.HIVposCells.Antibody.markers$p_val, method = "BH", n = nrow(heatmap.HIVposCells.Antibody.markers), pseudocount.use = 1)
heatmap.HIVposCells.Antibody.markers <- heatmap.HIVposCells.Antibody.markers[!duplicated(heatmap.HIVposCells.Antibody.markers$gene),]

avg_exp.HIVposCells <- AverageExpression(HIVposCells, return.seurat = T, group.by = 'seurat_clusters')
avg_exp.HIVposCells$orig.ident <- colnames(avg_exp.HIVposCells)
Idents(avg_exp.HIVpos.cite) <- "orig.ident"
avg_exp.HIVpos.cite <- ScaleData(avg_exp.HIVpos.cite)
DefaultAssay(avg_exp.HIVpos.cite) <- "CITE_dsb"; Idents(avg_exp.HIVpos.cite) <- "orig.ident"
Idents(avg_exp.HIVpos.cite) <- factor(Idents(avg_exp.HIVpos.cite), levels = c("IRF", "Cytotoxic" ,"AP-1", "MT"))

DoHeatmap(object = avg_exp.HIVpos.cite, features = rev(heatmap.HIVposCells.Antibody.markers$gene), group.bar.height = 0.02, size = 4, angle = 0, 
          draw.lines = FALSE) + scale_fill_gradientn(colors = c("blue", "white", "red"))

#Fig.6H
splitcorheatmap(corHIVpos_IRFMod_corHIVpos_MTMod, HIVpos_IRF, "yIRF_xMT_WGCNAmod_IRF") 

#FIg. 6
splitcorheatmap(corHIVpos_cytotoxicMod_corHIVpos_MTMod, HIVpos_cytotoxic, "ycytotoxic_xMT_WGCNAmod_cytotoxic") 

#Fig. 6J
splitcorheatmap(corHIVpos_AP1Mod_corHIVpos_MTMod, HIVpos_AP1, "yAP1_xMT_WGCNAmod_AP1")
