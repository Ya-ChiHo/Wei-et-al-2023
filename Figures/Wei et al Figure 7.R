library(Seurat)
library(Signac)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(reshape2)
library(stringr)

#require "Preprocessing.R"
#require ""recall ATAC peaks and build motifs.R"

#Fig. 7A
Combined$Condition_WNNcelltype_ATAC <- str_c(Combined$Condition,'_', Combined$WNNcelltype_ATAC)

Idents(Combined) <- "Condition_WNNcelltype_ATAC"
Idents(Combined) <- factor(Idents(Combined), levels = c("Healthy_Th1", "Viremic_Th1", "Suppressed_Th1",
                                                        "Healthy_Th2", "Viremic_Th2", "Suppressed_Th2",
                                                        "Healthy_Th17", "Viremic_Th17", "Suppressed_Th17",
                                                        "Healthy_Treg", "Viremic_Treg", "Suppressed_Treg",
                                                        "Healthy_Proliferating", "Viremic_Proliferating", "Suppressed_Proliferating",
                                                        "Healthy_memory", "Viremic_memory", "Suppressed_memory"))
DefaultAssay(Combined) <- "RNA"
p1 <- DotPlot(Combined, features = c("IKZF3")) 

Combined.ATAC$Condition_WNNcelltype_ATAC <- str_c(Combined.ATAC$Condition,'_', Combined.ATAC$WNNcelltype_ATAC)

Idents(Combined.ATAC) <- "Condition_WNNcelltype_ATAC"
Idents(Combined.ATAC) <- factor(Idents(Combined.ATAC), levels = c("Healthy_Th1", "Viremic_Th1", "Suppressed_Th1",
                                                        "Healthy_Th2", "Viremic_Th2", "Suppressed_Th2",
                                                        "Healthy_Th17", "Viremic_Th17", "Suppressed_Th17",
                                                        "Healthy_Treg", "Viremic_Treg", "Suppressed_Treg",
                                                        "Healthy_Proliferating", "Viremic_Proliferating", "Suppressed_Proliferating",
                                                        "Healthy_memory", "Viremic_memory", "Suppressed_memory"))
DefaultAssay(Combined.ATAC) <- "activities"
p2 <- DotPlot(Combined.ATAC, features = c("IKZF3"))
p1+p2

#Fig. 7B
Viremic$WNNcelltype_ATAC_HIV_RNA_DNA2 <- str_c(Viremic$WNNcelltype_ATAC,'_', Viremic$HIV_RNA_DNA2)
Idents(Viremic) <- "WNNcelltype_ATAC_HIV_RNA_DNA2"
Idents(Viremic) <- factor(Idents(Viremic), levels = c("Th1_negative", "Th1_DNA positive", "Th1_RNA and double positive", 
                                                      "Th2_negative", "Th2_DNA positive", "Th2_RNA and double positive",
                                                      "Th17_negative", "Th17_DNA positive", "Th17_RNA and double positive",
                                                      "Treg_negative", "Treg_DNA positive", "Treg_RNA and double positive",
                                                      "Proliferating_negative", "Proliferating_DNA positive", "Proliferating_RNA and double positive",
                                                      "memory_negative", "memory_DNA positive", "memory_RNA and double positive"))
DefaultAssay(Viremic) <- "RNA"; 
p2 <- DotPlot(Viremic, features = c("IKZF3"))

Viremic.ATAC$WNNcelltype_ATAC_HIV_RNA_DNA2 <- str_c(Viremic.ATAC$WNNcelltype_ATAC,'_', Viremic.ATAC$HIV_RNA_DNA2)
Idents(Viremic.ATAC) <- "WNNcelltype_ATAC_HIV_RNA_DNA2"
Idents(Viremic.ATAC) <- factor(Idents(Viremic.ATAC), levels = c("Th1_negative", "Th1_DNA positive", "Th1_RNA and double positive", 
                                                      "Th2_negative", "Th2_DNA positive", "Th2_RNA and double positive",
                                                      "Th17_negative", "Th17_DNA positive", "Th17_RNA and double positive",
                                                      "Treg_negative", "Treg_DNA positive", "Treg_RNA and double positive",
                                                      "Proliferating_negative", "Proliferating_DNA positive", "Proliferating_RNA and double positive",
                                                      "memory_negative", "memory_DNA positive", "memory_RNA and double positive"))
DefaultAssay(Viremic.ATAC) <- "activities"; 
p2 <- DotPlot(Viremic.ATAC, features = c("IKZF3"))
p1+p2

#Fig. 7C
DefaultAssay(Viremic.ATAC) <- "ATAC"; Idents(Viremic.ATAC) <- "HIV_RNA_DNA2"
CoveragePlot(Viremic.ATAC, region = "IKZF3") 

#Fig. 7D
Idents(Viremic.ATAC) <- "HIV_RNA_DNA2"; DefaultAssay(Viremic.ATAC) <- "activities"
Idents(Viremic) <- "HIV_RNA_DNA2"; DefaultAssay(Viremic) <- "RNA"

p1 <- ggplot(data.frame(IKZF3 = Viremic.ATAC[["activities"]]@data["IKZF3",], 
             cluster = factor(Idents(Viremic.ATAC), levels = c("DNA positive", "RNA and double positive", "negative"))), 
             aes(x = cluster, y = IKZF3)) + geom_violin(aes(fill = cluster), trim=FALSE, scale = "width") 
p2 <- ggplot(data.frame(IKZF3 = Viremic[["RNA"]]@data["IKZF3",], 
             cluster = factor(Idents(Viremic), levels = c("DNA positive", "RNA and double positive", "negative"))), 
             aes(x = cluster, y = IKZF3)) + geom_violin(aes(fill = cluster), trim=FALSE, scale = "width") 
p1+p2