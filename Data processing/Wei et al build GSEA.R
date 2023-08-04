library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(tidyverse)
library(fgsea)
library(msigdbr)
library(reshape2)
library(tibble)

####make GSEA sets####
m_df_H<- msigdbr(species = "Homo sapiens", category = "H")
m_df_H<- rbind(msigdbr(species = "Homo sapiens", category = "C2"), m_df_H)
m_df_H<- rbind(msigdbr(species = "Homo sapiens", category = "C7"), m_df_H)
fgsea_sets<- m_df_H %>% split(x = .$gene_symbol, f = .$gs_name)

####GSEA sets enriched in Viremic condition####
Idents(Combined) <- "Condition"; DefaultAssay(Combined) <- "RNA"
DEG.Combined.GSEA <- FindAllMarkers(object = Combined, only.pos = TRUE, min.pct = 0.1,
                                      logfc.threshold = 0.25, test.use = "wilcox")
DEG.Combined.GSEA$BH <- p.adjust(DEG.Combined.GSEA$p_val, method = "BH", n = nrow(DEG.Combined.GSEA))
DEG.Combined.GSEA <- DEG.Combined.GSEA[DEG.Combined.GSEA$BH < 0.05,]

DEG.Combined.GSEA.Viremic <- DEG.Combined.GSEA[DEG.Combined.GSEA$cluster == "Viremic",]
DEG.Combined.GSEA.Viremic <- DEG.Combined.GSEA.Viremic[!duplicated(DEG.Combined.GSEA.Viremic$gene),]
DEG.Combined.GSEA.Viremic<- DEG.Combined.GSEA.Viremic[order(DEG.Combined.GSEA.Viremic$gene.name),]

Rank.DEG.Combined.GSEA.Viremic <- DEG.Combined.GSEA.Viremic$avg_log2FC
names(Rank.DEG.Combined.GSEA.Viremic) <- rownames(DEG.Combined.GSEA.Viremic)
fgsea_Combined.GSEA.Viremic <- fgsea(pathways = fgsea_sets, 
                                                stats = Rank.DEG.Combined.GSEA.Viremic,
                                                eps   = 0.0, minSize=15, maxSize=500, scoreType = "pos")
#find top GSEA
head(fgsea_Combined.GSEA.Viremic[order(-abs(NES)), ], n=70)


####GSEA enriched in Viremic condition in Proliferating cluster###
Idents(Proliferating) <- "Condition"; DefaultAssay(Combined) <- "RNA"
DEG.Proliferating.GSEA <- FindAllMarkers(object = Proliferating, only.pos = TRUE, min.pct = 0.1,
                                    logfc.threshold = 0.25, test.use = "wilcox")
DEG.Proliferating.GSEA$BH <- p.adjust(DEG.Proliferating.GSEA$p_val, method = "BH", n = nrow(DEG.Proliferating.GSEA))
DEG.Proliferating.GSEA <- DEG.Combined.GSEA[DEG.Proliferating.GSEA$BH < 0.05,]

DEG.Proliferating.Sample.Viremic <- DEG.Proliferating.GSEA[DEG.Proliferating.GSEA$cluster == 'Viremic',]

GSEA.DEG.Proliferating.Sample.Viremic <- DEG.Proliferating.Sample.Viremic$avg_log2FC
names(GSEA.DEG.Proliferating.Sample.Viremic) <- rownames(DEG.Proliferating.Sample.Viremic)
fgsea_GSEA.DEG.Proliferating.Sample.Viremic <- fgsea(pathways = fgsea_sets, 
                                                     stats = GSEA.DEG.Proliferating.Sample.Viremic,
                                                     eps   = 0.0, minSize=15, maxSize=500, scoreType = "pos")
#find top GSEA
head(fgsea_GSEA.DEG.Proliferating.Sample.Viremic[order(-abs(NES)), ], n=70)


####GSEA enriched in Suppressed condition in Proliferating cluster####
DEG.Proliferating.Sample.Suppressed <- DEG.Proliferating.GSEA[DEG.Proliferating.GSEA$cluster == 'Suppressed',]

GSEA.DEG.Proliferating.Sample.Suppressed <- DEG.Proliferating.Sample.Suppressed$avg_log2FC
names(GSEA.DEG.Proliferating.Sample.Suppressed) <- rownames(DEG.Proliferating.Sample.Suppressed)
fgsea_GSEA.DEG.Proliferating.Sample.Suppressed <- fgsea(pathways = fgsea_sets, 
                                                        stats = GSEA.DEG.Proliferating.Sample.Suppressed,
                                                        eps   = 0.0, minSize=15, maxSize=500, scoreType = "pos")
#find top GSEA
head(fgsea_GSEA.DEG.Proliferating.Sample.Suppressed[order(-abs(NES)), ], n=70)
