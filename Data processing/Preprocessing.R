library(Seurat)
library(Signac)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(dsb)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(GenomeInfoDb)
library(TFBSTools)
library(JASPAR2022)
library(harmony)
library(clustree)
library(stringr)
library(future)

#require "HIV_RNApos_cells.csv"
#require "HIV_DNApos_cells.csv"
#require "bc_ATAC_HIV.txt"
#require "bc_RNA_HIV.txt"

####build Seurat objects(RNA + protein)####
HD1 <- Read10X("~/HD1_MAH_cite/filtered_feature_bc_matrix")
YW1 <- Read10X("~/YW1_MAH_cite/filtered_feature_bc_matrix")
YW2 <- Read10X("~/YW2_MAH_cite/filtered_feature_bc_matrix")
YW3 <- Read10X("~/YW3_MAH_cite/filtered_feature_bc_matrix")
YW4 <- Read10X("~/YW4_MAH_cite/filtered_feature_bc_matrix")
YW5 <- Read10X("~/YW5_MAH_cite/filtered_feature_bc_matrix")
YW6 <- Read10X("~/YW6_MAH_cite/filtered_feature_bc_matrix")
YW8 <- Read10X("~/YW8_MAH_cite/filtered_feature_bc_matrix")
YW9 <- Read10X("~/YW9_MAH_cite/filtered_feature_bc_matrix")
YW10 <- Read10X("~/YW10_MAH_cite/filtered_feature_bc_matrix")

RNA_HD1 <- HD1$'Gene Expression'; antibody_HD1 <- HD1$'Antibody Capture'
RNA_YW1 <- YW1$'Gene Expression'; antibody_YW1 <- YW1$'Antibody Capture'
RNA_YW2 <- YW2$'Gene Expression'; antibody_YW2 <- YW2$'Antibody Capture'
RNA_YW3 <- YW3$'Gene Expression'; antibody_YW3 <- YW3$'Antibody Capture'
RNA_YW4 <- YW4$'Gene Expression'; antibody_YW4 <- YW4$'Antibody Capture'
RNA_YW5 <- YW5$'Gene Expression'; antibody_YW5 <- YW5$'Antibody Capture'
RNA_YW6 <- YW6$'Gene Expression'; antibody_YW6 <- YW6$'Antibody Capture'
RNA_YW8 <- YW8$'Gene Expression'; antibody_YW8 <- YW8$'Antibody Capture'
RNA_YW9 <- YW9$'Gene Expression'; antibody_YW9 <- YW9$'Antibody Capture'
RNA_YW10 <- YW10$'Gene Expression'; antibody_YW10 <- YW10$'Antibody Capture'

HD1 <- CreateSeuratObject(counts = RNA_HD1); antibody_assay_HD1 <- CreateAssayObject(counts = antibody_HD1); HD1[["Antibody"]] <- antibody_assay_HD1
YW1 <- CreateSeuratObject(counts = RNA_YW1); antibody_assay_YW1 <- CreateAssayObject(counts = antibody_YW1); YW1[["Antibody"]] <- antibody_assay_YW1
YW2 <- CreateSeuratObject(counts = RNA_YW2); antibody_assay_YW2 <- CreateAssayObject(counts = antibody_YW2); YW2[["Antibody"]] <- antibody_assay_YW2
YW3 <- CreateSeuratObject(counts = RNA_YW3); antibody_assay_YW3 <- CreateAssayObject(counts = antibody_YW3); YW3[["Antibody"]] <- antibody_assay_YW3
YW4 <- CreateSeuratObject(counts = RNA_YW4); antibody_assay_YW4 <- CreateAssayObject(counts = antibody_YW4); YW4[["Antibody"]] <- antibody_assay_YW4
YW5 <- CreateSeuratObject(counts = RNA_YW5); antibody_assay_YW5 <- CreateAssayObject(counts = antibody_YW5); YW5[["Antibody"]] <- antibody_assay_YW5
YW6 <- CreateSeuratObject(counts = RNA_YW6); antibody_assay_YW6 <- CreateAssayObject(counts = antibody_YW6); YW6[["Antibody"]] <- antibody_assay_YW6
YW8 <- CreateSeuratObject(counts = RNA_YW8); antibody_assay_YW8 <- CreateAssayObject(counts = antibody_YW8); YW8[["Antibody"]] <- antibody_assay_YW8
YW9 <- CreateSeuratObject(counts = RNA_YW9); antibody_assay_YW9 <- CreateAssayObject(counts = antibody_YW9); YW9[["Antibody"]] <- antibody_assay_YW9
YW10 <- CreateSeuratObject(counts = RNA_YW10); antibody_assay_YW10 <- CreateAssayObject(counts = antibody_YW10); YW10[["Antibody"]] <- antibody_assay_YW10

####add per cell percentage of mitochondrial gene to objects####
HD1[["percent.mt"]] <- PercentageFeatureSet(HD1, pattern = "^MT-")
YW1[["percent.mt"]] <- PercentageFeatureSet(YW1, pattern = "^MT-")
YW2[["percent.mt"]] <- PercentageFeatureSet(YW2, pattern = "^MT-")
YW3[["percent.mt"]] <- PercentageFeatureSet(YW3, pattern = "^MT-")
YW4[["percent.mt"]] <- PercentageFeatureSet(YW4, pattern = "^MT-")
YW5[["percent.mt"]] <- PercentageFeatureSet(YW5, pattern = "^MT-")
YW6[["percent.mt"]] <- PercentageFeatureSet(YW6, pattern = "^MT-")
YW8[["percent.mt"]] <- PercentageFeatureSet(YW8, pattern = "^MT-")
YW9[["percent.mt"]] <- PercentageFeatureSet(YW9, pattern = "^MT-")
YW10[["percent.mt"]] <- PercentageFeatureSet(YW10, pattern = "^MT-")

####add to each objects a hashtag assay####
hashtag_HD1 <- antibody_HD1[rownames(antibody_HD1) %in% c("hashtag1_TotalA", "hashtag2_TotalA", "hashtag3_TotalA", "hashtag4_TotalA", "hashtag5_TotalA", "hashtag6_TotalA"),]
hashtag_YW1 <- antibody_YW1[rownames(antibody_YW1) %in% c("hashtag1_TotalA", "hashtag2_TotalA", "hashtag3_TotalA", "hashtag4_TotalA", "hashtag5_TotalA", "hashtag6_TotalA"),]
hashtag_YW2 <- antibody_YW2[rownames(antibody_YW2) %in% c("hashtag1_TotalA", "hashtag2_TotalA", "hashtag3_TotalA", "hashtag4_TotalA", "hashtag5_TotalA", "hashtag6_TotalA"),]
hashtag_YW3 <- antibody_YW3[rownames(antibody_YW3) %in% c("hashtag1_TotalA", "hashtag2_TotalA", "hashtag3_TotalA", "hashtag4_TotalA", "hashtag5_TotalA", "hashtag6_TotalA"),]
hashtag_YW4 <- antibody_YW4[rownames(antibody_YW4) %in% c("hashtag1_TotalA", "hashtag2_TotalA", "hashtag3_TotalA", "hashtag4_TotalA", "hashtag5_TotalA", "hashtag6_TotalA"),]
hashtag_YW5 <- antibody_YW5[rownames(antibody_YW5) %in% c("hashtag1_TotalA", "hashtag2_TotalA", "hashtag3_TotalA", "hashtag4_TotalA", "hashtag5_TotalA", "hashtag6_TotalA"),]
hashtag_YW6 <- antibody_YW6[rownames(antibody_YW6) %in% c("hashtag1_TotalA", "hashtag2_TotalA", "hashtag3_TotalA", "hashtag4_TotalA", "hashtag5_TotalA", "hashtag6_TotalA"),]
hashtag_YW8 <- antibody_YW8[rownames(antibody_YW8) %in% c("hashtag1_TotalA", "hashtag2_TotalA", "hashtag3_TotalA", "hashtag4_TotalA", "hashtag5_TotalA", "hashtag6_TotalA"),]
hashtag_YW9 <- antibody_YW9[rownames(antibody_YW9) %in% c("hashtag1_TotalA", "hashtag2_TotalA", "hashtag3_TotalA", "hashtag4_TotalA", "hashtag5_TotalA", "hashtag6_TotalA"),]
hashtag_YW10 <- antibody_YW10[rownames(antibody_YW10) %in% c("hashtag1_TotalA", "hashtag2_TotalA", "hashtag3_TotalA", "hashtag4_TotalA", "hashtag5_TotalA", "hashtag6_TotalA"),]

HD1[["hashtag"]] <- CreateAssayObject(counts = hashtag_ID)
YW1[["hashtag"]] <- CreateAssayObject(counts = hashtag_YW1)
YW2[["hashtag"]] <- CreateAssayObject(counts = hashtag_YW2)
YW3[["hashtag"]] <- CreateAssayObject(counts = hashtag_YW3)
YW4[["hashtag"]] <- CreateAssayObject(counts = hashtag_YW4)
YW5[["hashtag"]] <- CreateAssayObject(counts = hashtag_YW5)
YW6[["hashtag"]] <- CreateAssayObject(counts = hashtag_YW6)
YW8[["hashtag"]] <- CreateAssayObject(counts = hashtag_YW8)
YW9[["hashtag"]] <- CreateAssayObject(counts = hashtag_YW9)
YW10[["hashtag"]] <- CreateAssayObject(counts = hashtag_YW10)

####detecting doublets####
HD1 <- NormalizeData(HD1, assay = "hashtag", normalization.method = "CLR")
YW1 <- NormalizeData(YW1, assay = "hashtag", normalization.method = "CLR")
YW2 <- NormalizeData(YW2, assay = "hashtag", normalization.method = "CLR")
YW3 <- NormalizeData(YW3, assay = "hashtag", normalization.method = "CLR")
YW4 <- NormalizeData(YW4, assay = "hashtag", normalization.method = "CLR")
YW5 <- NormalizeData(YW5, assay = "hashtag", normalization.method = "CLR")
YW6 <- NormalizeData(YW6, assay = "hashtag", normalization.method = "CLR")
YW8 <- NormalizeData(YW8, assay = "hashtag", normalization.method = "CLR")
YW9 <- NormalizeData(YW9, assay = "hashtag", normalization.method = "CLR")
YW10 <- NormalizeData(YW10, assay = "hashtag", normalization.method = "CLR")

HD1 <- MULTIseqDemux(HD1, assay = "hashtag",  autoThresh  = TRUE, maxiter = 10, qrange = seq(from = 0.1, to = 0.999, by = 0.001), verbose = FALSE)
YW1 <- MULTIseqDemux(YW1, assay = "hashtag",  autoThresh  = TRUE, maxiter = 10, qrange = seq(from = 0.1, to = 0.999, by = 0.001), verbose = FALSE)
YW2 <- MULTIseqDemux(YW2, assay = "hashtag",  autoThresh  = TRUE, maxiter = 10, qrange = seq(from = 0.1, to = 0.999, by = 0.001), verbose = FALSE)
YW3 <- MULTIseqDemux(YW3, assay = "hashtag",  autoThresh  = TRUE, maxiter = 10, qrange = seq(from = 0.1, to = 0.999, by = 0.001), verbose = FALSE)
YW4 <- MULTIseqDemux(YW4, assay = "hashtag",  autoThresh  = TRUE, maxiter = 10, qrange = seq(from = 0.1, to = 0.999, by = 0.001), verbose = FALSE)
YW5 <- MULTIseqDemux(YW5, assay = "hashtag",  autoThresh  = TRUE, maxiter = 10, qrange = seq(from = 0.1, to = 0.999, by = 0.001), verbose = FALSE)
YW6 <- MULTIseqDemux(YW6, assay = "hashtag",  autoThresh  = TRUE, maxiter = 10, qrange = seq(from = 0.1, to = 0.999, by = 0.001), verbose = FALSE)
YW8 <- MULTIseqDemux(YW8, assay = "hashtag",  autoThresh  = TRUE, maxiter = 10, qrange = seq(from = 0.1, to = 0.999, by = 0.001), verbose = FALSE)
YW9 <- MULTIseqDemux(YW9, assay = "hashtag",  autoThresh  = TRUE, maxiter = 10, qrange = seq(from = 0.1, to = 0.999, by = 0.001), verbose = FALSE)
YW10 <- MULTIseqDemux(YW10, assay = "hashtag",  autoThresh  = TRUE, maxiter = 10, qrange = seq(from = 0.1, to = 0.999, by = 0.001), verbose = FALSE)


####merge objects and remove doublets####
Combined.cite <- merge(x = HD1, y = list(YW1, YW2, YW3, YW4, YW5, YW6, YW8, YW9, YW10),
                            add.cell.ids = c("HD1", "YW1", "YW2", "YW3", "YW4", "YW5", "YW6", "YW8", "YW9", "YW10" ))
Combined.cite$orig.ident <- unlist(sapply(X = strsplit(colnames(Combined.cite), split = "_"), FUN = "[", 1))

Combined.cite$MULTI_ID.global <- Combined.cite$MULTI_ID
Idents(Combined.cite) <- "MULTI_ID.global"
Combined.cite <-RenameIdents(Combined.cite, 'hashtag1-TotalA' = 'Singlet', 'hashtag2-TotalA' = 'Singlet', 
                                  'hashtag3-TotalA' = 'Singlet', 'hashtag4-TotalA' = 'Singlet','hashtag5-TotalA' = 'Singlet',
                                  'hashtag6-TotalA' = 'Singlet')
Combined.cite$MULTI_ID.global <- Idents(Combined.cite)
Idents(Combined.cite) <- "MULTI_ID.global"
Combined.cite <-subset(Combined.cite, idents = c("Singlet"))

####adding HIV RNA and HIV DNA positive cells to Object####
#See Supplemental Table S3 (Wei et al. Immunity, 2023) for HIV-1-infected cell barcodes (with HIV read copies per barcode >= 2)
HIV_DNA <- read.table("~/bc_ATAC_HIV.txt", head = T)
HIV_RNA <- read.table("~/bc_RNA_HIV.txt", head = T)

Combined.cite$HIV_DNA <- ifelse(colnames(Combined.cite) %in% HIV_DNA$cell, "positive", "negative")
Combined.cite$HIV_DNA_copies <- ifelse(colnames(Combined.cite) %in% HIV_DNA$cell, HIV_DNA$copies[match(colnames(Combined.cite), HIV_DNA$cell)],0)

Combined.cite$HIV_RNA <- ifelse(colnames(Combined.cite) %in% HIV_RNA$cell, "positive", "negative")
Combined.cite$HIV_RNA_copies <- ifelse(colnames(Combined.cite) %in% HIV_RNA$cell, HIV_RNA$Copies[match(colnames(Combined.cite), HIV_RNA$cell)],0)

Combined.cite$HIV_RNA_DNA <- str_c(Combined.cite$HIV_RNA, '_', Combined.cite$HIV_DNA)

####remove low quality cells by RNA metrics####
Combined.cite <- subset(Combined.cite, nFeature_RNA > 200 & percent.mt < 27 & nCount_RNA < 10000 & nCount_RNA > 500)

####compute DSB-normalized protein expressions####
DefaultAssay(Combined.cite) <- "Antibody"
Combined_cite_mtx <- Combined.cite@assays$Antibody@counts
dsb.norm <- ModelNegativeADTnorm(cell_protein_matrix =  Combined_cite_mtx,
                                 denoise.counts = TRUE, use.isotype.control = TRUE, 
                                 isotype.control.name.vec = c("Mouse-IgG1-isotype-TotalA", 
                                                              "Mouse-IgG2a-isotype-TotalA", "Mouse-IgG2b-isotype-TotalA",
                                                              "Mouse-IgG2b-isotype-TotalA-1", "RatIgG1-TotalA", "RatIgG1-TotalA-i",
                                                              "RatIgG2a-TotalA", "RatIgG2c-TotalA", "ArmenianHamsterIgG-TotalA",
                                                              "RatIgG1k-TotalA", "RatIgG1l-TotalA", "Rat-IgG2b-isotype-TotalA"))
Combined.cite[["CITE_dsb"]] = Seurat::CreateAssayObject(data = dsb.norm)

####make common peak set, create ATAC Seurat Objects, and merge Objects####
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2)

HD1.peaks <- read.table(file = "~/HD1_MAH_cellranger/atac_peaks.bed", col.names = c("chr", "start", "end"))
YW1.peaks <- read.table(file = "~/YW1_MAH_cellranger/atac_peaks.bed", col.names = c("chr", "start", "end"))
YW2.peaks <- read.table(file = "~/YW2_MAH_cellranger/atac_peaks.bed", col.names = c("chr", "start", "end"))
YW3.peaks <- read.table(file = "~/YW3_MAH_cellranger/atac_peaks.bed", col.names = c("chr", "start", "end"))
YW4.peaks <- read.table(file = "~/YW4_MAH_cellranger/atac_peaks.bed", col.names = c("chr", "start", "end"))
YW5.peaks <- read.table(file = "~/YW5_MAH_cellranger/atac_peaks.bed", col.names = c("chr", "start", "end"))
YW6.peaks <- read.table(file = "~/YW6_MAH_cellranger/atac_peaks.bed", col.names = c("chr", "start", "end"))
YW8.peaks <- read.table(file = "~/YW8_MAH_cellranger/atac_peaks.bed", col.names = c("chr", "start", "end"))
YW9.peaks <- read.table(file = "~/YW9_MAH_cellranger/atac_peaks.bed", col.names = c("chr", "start", "end"))
YW10.peaks <- read.table(file = "~/YW10_MAH_cellranger/atac_peaks.bed", col.names = c("chr", "start", "end"))

HD1.gr <- makeGRangesFromDataFrame(HD1.peaks)
YW1.gr <- makeGRangesFromDataFrame(YW1.peaks)
YW2.gr <- makeGRangesFromDataFrame(YW2.peaks)
YW3.gr <- makeGRangesFromDataFrame(YW3.peaks)
YW4.gr <- makeGRangesFromDataFrame(YW4.peaks)
YW5.gr <- makeGRangesFromDataFrame(YW5.peaks)
YW6.gr <- makeGRangesFromDataFrame(YW6.peaks)
YW8.gr <- makeGRangesFromDataFrame(YW8.peaks)
YW9.gr <- makeGRangesFromDataFrame(YW9.peaks)
YW10.gr <- makeGRangesFromDataFrame(YW10.peaks)

combined.peaks <- reduce(x = c(HD1.gr, YW1.gr, YW2.gr, YW3.gr, YW4.gr, YW5.gr, YW6.gr, YW8.gr, YW9.gr, YW10.gr))
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]

HD1.md <- read.table(file = "~/HD1_MAH_cellranger/per_barcode_metrics.csv",
                     stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)[-1,]
YW1.md <- read.table(file = "~/YW1_MAH_cellranger/per_barcode_metrics.csv",
                     stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)[-1,]
YW2.md <- read.table(file = "~/YW2_MAH_cellranger/per_barcode_metrics.csv",
                     stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)[-1,]
YW3.md <- read.table(file = "~/YW3_MAH_cellranger/per_barcode_metrics.csv",
                     stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)[-1,]
YW4.md <- read.table(file = "~/YW4_MAH_cellranger/per_barcode_metrics.csv",
                     stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)[-1,]
YW5.md <- read.table(file = "~/YW5_MAH_cellranger/per_barcode_metrics.csv",
                     stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)[-1,]
YW6.md <- read.table(file = "~/YW6_MAH_cellranger/per_barcode_metrics.csv",
                     stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)[-1,]
YW8.md <- read.table(file = "~/YW8_MAH_cellranger/per_barcode_metrics.csv",
                     stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)[-1,]
YW9.md <- read.table(file = "~/YW9_MAH_cellranger/per_barcode_metrics.csv",
                     stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)[-1,]
YW10.md <- read.table(file = "~/YW10_MAH_cellranger/per_barcode_metrics.csv",
                      stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)[-1,]

HD1.frags <- CreateFragmentObject(path = "~/HD1_MAH_cellranger/atac_fragments.tsv.gz", cells = rownames(HD1.md))
YW1.frags <- CreateFragmentObject(path = "~/YW1_MAH_cellranger/atac_fragments.tsv.gz", cells = rownames(YW1.md))
YW2.frags <- CreateFragmentObject(path = "~/YW2_MAH_cellranger/atac_fragments.tsv.gz", cells = rownames(YW2.md))
YW3.frags <- CreateFragmentObject(path = "~/YW3_MAH_cellranger/atac_fragments.tsv.gz", cells = rownames(YW3.md))
YW4.frags <- CreateFragmentObject(path = "~/YW4_MAH_cellranger/atac_fragments.tsv.gz", cells = rownames(YW4.md))
YW5.frags <- CreateFragmentObject(path = "~/YW5_MAH_cellranger/atac_fragments.tsv.gz", cells = rownames(YW5.md))
YW6.frags <- CreateFragmentObject(path = "~/YW6_MAH_cellranger/atac_fragments.tsv.gz", cells = rownames(YW6.md))
YW8.frags <- CreateFragmentObject(path = "~/YW8_MAH_cellranger/atac_fragments.tsv.gz", cells = rownames(YW8.md))
YW9.frags <- CreateFragmentObject(path = "~/YW9_MAH_cellranger/atac_fragments.tsv.gz", cells = rownames(YW9.md))
YW10.frags <- CreateFragmentObject(path = "~/YW10_MAH_cellranger/atac_fragments.tsv.gz", cells = rownames(YW10.md))

HD1.counts <- FeatureMatrix(fragments = HD1.frags, features = combined.peaks, cells = rownames(HD1.md))
YW1.counts <- FeatureMatrix(fragments = YW1.frags, features = combined.peaks, cells = rownames(YW1.md))
YW2.counts <- FeatureMatrix(fragments = YW2.frags, features = combined.peaks, cells = rownames(YW2.md))
YW3.counts <- FeatureMatrix(fragments = YW3.frags, features = combined.peaks, cells = rownames(YW3.md))
YW4.counts <- FeatureMatrix(fragments = YW4.frags, features = combined.peaks, cells = rownames(YW4.md))
YW5.counts <- FeatureMatrix(fragments = YW5.frags, features = combined.peaks, cells = rownames(YW5.md))
YW6.counts <- FeatureMatrix(fragments = YW6.frags, features = combined.peaks, cells = rownames(YW6.md))
YW8.counts <- FeatureMatrix(fragments = YW8.frags, features = combined.peaks, cells = rownames(YW8.md))
YW9.counts <- FeatureMatrix(fragments = YW9.frags, features = combined.peaks, cells = rownames(YW9.md))
YW10.counts <- FeatureMatrix(fragments = YW10.frags, features = combined.peaks, cells = rownames(YW10.md))

HD1_chrom_assay <- CreateChromatinAssay(counts = HD1.counts, fragments = HD1.frags)
YW1_chrom_assay <- CreateChromatinAssay(counts = YW1.counts, fragments = YW1.frags)
YW2_chrom_assay <- CreateChromatinAssay(counts = YW2.counts, fragments = YW2.frags)
YW3_chrom_assay <- CreateChromatinAssay(counts = YW3.counts, fragments = YW3.frags)
YW4_chrom_assay <- CreateChromatinAssay(counts = YW4.counts, fragments = YW4.frags)
YW5_chrom_assay <- CreateChromatinAssay(counts = YW5.counts, fragments = YW5.frags)
YW6_chrom_assay <- CreateChromatinAssay(counts = YW6.counts, fragments = YW6.frags)
YW8_chrom_assay <- CreateChromatinAssay(counts = YW8.counts, fragments = YW8.frags)
YW9_chrom_assay <- CreateChromatinAssay(counts = YW9.counts, fragments = YW9.frags)
YW10_chrom_assay <- CreateChromatinAssay(counts = YW10.counts, fragments = YW10.frags)

HD1 <- CreateSeuratObject(counts = HD1_chrom_assay, assay = "ATAC", meta.data = HD1.md)
YW1 <- CreateSeuratObject(counts = YW1_chrom_assay, assay = "ATAC", meta.data = YW1.md)
YW2 <- CreateSeuratObject(counts = YW2_chrom_assay, assay = "ATAC", meta.data = YW2.md)
YW3 <- CreateSeuratObject(counts = YW3_chrom_assay, assay = "ATAC", meta.data = YW3.md)
YW4 <- CreateSeuratObject(counts = YW4_chrom_assay, assay = "ATAC", meta.data = YW4.md)
YW5 <- CreateSeuratObject(counts = YW5_chrom_assay, assay = "ATAC", meta.data = YW5.md)
YW6 <- CreateSeuratObject(counts = YW6_chrom_assay, assay = "ATAC", meta.data = YW6.md)
YW8 <- CreateSeuratObject(counts = YW8_chrom_assay, assay = "ATAC", meta.data = YW8.md)
YW9 <- CreateSeuratObject(counts = YW9_chrom_assay, assay = "ATAC", meta.data = YW9.md)
YW10 <- CreateSeuratObject(counts = YW10_chrom_assay, assay = "ATAC", meta.data = YW10.md)

HD1$orig.ident <- 'HD1'; YW1$orig.ident <- 'YW1'; YW2$orig.ident <- 'YW2'; YW3$orig.ident <- 'YW3'; YW4$orig.ident <- 'YW4'; YW5$orig.ident <- 'YW5'; YW6$orig.ident <- 'YW6'; YW8$orig.ident <- 'YW8'; YW9$orig.ident <- 'YW9'; YW10$orig.ident <- 'YW10'

Combined.ATAC <- merge(x = HD1, y = list(YW1, YW2, YW3, YW4, YW5, YW6, YW8, YW9, YW10), 
                       add.cell.ids = c("HD1", "YW1", "YW2", "YW3", "YW4", "YW5", "YW6", "YW8", "YW9", "YW10"))

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(Combined.ATAC) <- annotations

####remove poor quality cells by ATAC metrics####
Combined.ATAC <- NucleosomeSignal(object = Combined.ATAC)
Combined.ATAC <- TSSEnrichment(object = Combined.ATAC, fast = FALSE)

Idents(Combined.ATAC) <- "is_cell"
Combined.ATAC <-subset(Combined.ATAC, idents = c("1"))
Combined.ATAC$high.tss <- ifelse(Combined.ATAC$TSS.enrichment > 2, 'High', 'Low')
Combined.ATAC$nucleosome_group <- ifelse(Combined.ATAC$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

Combined.ATAC <- subset(Combined.ATAC, subset = nCount_ATAC < 100000 & nCount_ATAC > 500 & nucleosome_signal < 1 & TSS.enrichment > 2)

####combine RNA + protein Seurat Object with ATAC Seurat Object, keeping only cell barcodes that pass QC metrics by both RNA and ATAC####
all.equal(colnames(Combined.ATAC), colnames(Combined.cite))
joint.bcs<- intersect(colnames(Combined.ATAC), colnames(Combined.cite))
Idents(Combined.ATAC) <- colnames(Combined.ATAC)
Combined <-subset(Combined.ATAC, idents = joint.bcs)
Idents(Combined.cite) <- colnames(Combined.cite)
Combined.cite <- subset(Combined.cite, idents = joint.bcs)

Combined[["RNA"]] <- Combined.cite[["RNA"]] 
Combined[["Antibody"]] <- Combined.cite[["Antibody"]] 

Combined$HIV_DNA <- Combined.cite$HIV_DNA
Combined$HIV_DNA_copies <- Combined.cite$HIV_DNA_copies 
Combined$HIV_RNA <- Combined.cite$HIV_RNA  
Combined$HIV_RNA_copies <- Combined.cite$HIV_RNA_copies 
Combined$HIV_RNA_DNA <- Combined.cite$HIV_RNA_DNA  
Combined$MULTI_classification <- Combined.cite$MULTI_classification


####ATAC integration (batch effect correction) and normalization and UMAP generation####
DefaultAssay(Combined) <- "ATAC"
Combined <- RunTFIDF(Combined)
Combined <- FindTopFeatures(Combined, min.cutoff = 20)
Combined <- RunSVD(Combined)
Combined <- RunUMAP(Combined, dims = 2:30, reduction = 'lsi',
                         reduction.name = "umap.atac", reduction.key = "atacUMAP_")
DepthCor(Combined)
Idents(Combined) <- "orig.ident"
#correct batch effects across 10X runs
Combined.atac_anchors <- FindIntegrationAnchors(object.list = SplitObject(Combined, split.by = "orig.ident"),
                                                anchor.features = rownames(subset(Combined, idents = c("HD1"))), 
                                                reduction = "rlsi", dims = 2:30)
Combined.ATAC.integrated <- IntegrateEmbeddings(anchorset = Combined.atac_anchors, 
                                                reductions = Combined[["lsi"]],
                                                new.reduction.name = "integrated_lsi", 
                                                dims.to.integrate = 1:30) 

Combined.ATAC.integrated <- RunUMAP(Combined.ATAC.integrated, reduction = "integrated_lsi", dims = 2:30, 
                                    reduction.name = "umap.atac", reduction.key = "atacUMAP_")

Combined.ATAC.integrated <- RunTFIDF(Combined.ATAC.integrated)
Combined.ATAC.integrated <- FindTopFeatures(Combined.ATAC.integrated, min.cutoff = 20)
Combined.ATAC.integrated <- RunSVD(Combined.ATAC.integrated)

####add motifs to integrated ATAC object####
pfm <- getMatrixSet(x = JASPAR2022, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
DefaultAssay(Combined.ATAC.integrated) <- "ATAC"
keep.peaks <- as.logical(seqnames(granges(Combined.ATAC.integrated)) %in% main.chroms)
Combined.ATAC.integrated <- Combined.ATAC.integrated[keep.peaks,]
Combined.ATAC.integrated <- AddMotifs(Combined.ATAC.integrated, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)

####add gene accessibility score to integrated ATAC object####
gene.activities <- GeneActivity(Combined.ATAC.integrated)
Combined.ATAC.integrated[['activities']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(Combined.ATAC.integrated) <- "activities"
Combined.ATAC.integrated <- NormalizeData(object = Combined.ATAC.integrated, assay = 'activities', 
                                          normalization.method = 'LogNormalize', 
                                          scale.factor = median(Combined.ATAC.integrated$nCount_activities))
Combined.ATAC.integrated <- FindVariableFeatures(Combined.ATAC.integrated, selection.method = 'vst', nfeatures = 2000)
Combined.ATAC.integrated <- ScaleData(Combined.ATAC.integrated)


####protein integration (batch effect correction) and normalization and UMAP generation####
DefaultAssay(Combined) <- "Antibody"
rownames(Combined)
#remove hashtag and isotype antibodies from list of antibodies to make list of variable features
antibody_features <- rownames(Combined)[-c(31,32,33,34,85,86,87,88,89,166,167,168,169,170,171,172,173,174)]

DefaultAssay(Combined) <- "Antibody"
Idents(Combined) <- "orig.ident"
#correct batch effects across 10X runs
Combine.list <- SplitObject(Combined, split.by = "orig.ident")
Combine.list <- lapply(X = Combine.list, FUN = function(x){
  x <- NormalizeData(x, normalization.method = 'CLR', margin = 2)
  x <- ScaleData(x, features = antibody_features, verbose = FALSE)
  x <- RunPCA(x, features = antibody_features, verbose = FALSE)
})

immune.anchors <- FindIntegrationAnchors(object.list = Combine.list,
                                         anchor.features = antibody_features,
                                         reduction = "rpca", k.anchor = 20)
Combined <- IntegrateData(anchorset = immune.anchors, new.assay.name = "antibodyintegrated")

DefaultAssay(Combined) <- "antibodyintegrated"
Combined <- ScaleData(Combined, verbose = FALSE)
Combined <- RunPCA(Combined, npcs = 15, verbose = FALSE, reduction.name = 'apca', assay = 'antibodyintegrated')
Combined <- RunUMAP(Combined, reduction = "apca", reduction.name = "umap.antibody", 
                    reduction.key = "antibodyUMAP_", dims = 1:15, assay = 'antibodyintegrated')
Combined <- FindNeighbors(Combined, reduction = "apca", dims = 1:15)
Combined <- FindClusters(Combined, algorithm = 3, resolution = 0.5)

####RNA integration (batch effect correction) and normalization and UMAP generation####
DefaultAssay(Combined) <- "RNA"
#if using Seurat v5: Combined <- JoinLayers(Combined)

#correct batch effects across 10X runs
Idents(Combined) <- "orig.ident"
Combined <- NormalizeData(Combined, normalization.method = 'LogNormalize') %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>% RunPCA(npcs = 30) %>% RunHarmony("orig.ident", plot_convergence = TRUE) %>%
  RunUMAP(reduction = "harmony", reduction.name = "umap.rna", reduction.key = "rnaUMAP_", dims = 1:30, assay = 'RNA') %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.5) %>%
  identity()

####transfer integrated ATAC object to integrated RNA and protein object####
Combined[["ATAC"]] <- Combined.ATAC.integrated[["ATAC"]]
Combined[["activities"]] <- Combined.ATAC.integrated[["activities"]]

Combined@reductions$integrated_lsi <- Combined.ATAC.integrated@reductions$integrated_lsi
Combined@reductions$umap.atac <- Combined.ATAC.integrated@reductions$umap.atac
Combined@reductions$lsi <- Combined.ATAC.integrated@reductions$lsi

####generate 3WNN UMAP####
Combined <- FindMultiModalNeighbors(Combined, reduction.list = list("harmony", "apca", "integrated_lsi"), 
                                    dims.list = list(1:13, 1:10, 2:30), modality.weight.name = "3WNN.weight")
Combined <- RunUMAP(Combined, nn.name = "weighted.nn", reduction.name = "umap.wnn", reduction.key = "wnnUMAP_")
Combined <- FindClusters(Combined, graph.name = "wsnn", algorithm = 3, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), verbose = TRUE)
#use clustree to find best number of clusters
clustree(Combined, prefix = "wsnn_res.")

Combined$seurat_clusters <- Combined$wsnn_res.0.6
Idents(Combined) <- "seurat_clusters"
table(Combined$seurat_clusters)
#group clusters with <= 2 cells into nearest large clusters
#e.g., Combined <-RenameIdents(Combined, '33' = '0') if cluster 33 has 2 cells and nearest large cluster is 0

DimPlot(Combined, reduction = "umap.rna", group.by = "orig.ident", label = FALSE, raster = FALSE) + ggtitle("RNA UMAP res 0.6") 
DimPlot(Combined, reduction = "umap.antibody", group.by = "orig.ident", label = FALSE, raster = FALSE) + ggtitle("Antibody UMAP res 0.6") 
DimPlot(Combined, reduction = "umap.atac", group.by = "orig.ident", label = FALSE, raster = FALSE) + ggtitle("Antibody UMAP res 0.6") 
DimPlot(Combined, reduction = "umap.wnn", group.by = "orig.ident", label = FALSE, raster = FALSE) + ggtitle("WNN UMAP res 0.6")

####manual annotation by marker gene expression (RNA)####
DefaultAssay(Combined) <- "RNA"

#refer to Table S2 "cluster annotation" in Wei et al. 2023, Immunity for celltypes associated with markers below
CD4.RNA.markers <- c("percent.mt", "CCR7", "HLA-DRA", "HLA-DRB1","SELL",
                     "MKI67", 
                     "CTLA4", "PDCD1",
                     "IL2RA", "FOXP3", "IKZF2",
                     "CCR6", "ITGAE", 
                     "GATA3",
                     "TBX21","CCL5","GNLY", "GZMK", "GZMB", "GZMH")

DotPlot(Combined, features = rev(CD4.RNA.markers), dot.scale = 20, cols = c("grey","blue")) + scale_y_discrete(limits = rev)
#rename/merge cluster to celltypes in new metadata column "WNNcelltype_RNA" by high RNA expression of marker genes above, refer to Supplemental Figure S1E, S1H in Wei et al. 2023, Immunity. 

####manual annotation by transcription factor gene accessibility (ATAC)####
DefaultAssay(Combined) <- "activities"

#refer to Table S2 "cluster annotation" in Wei et al. 2023, Immunity for celltypes associated with markers below
CD4.ATAC.markers <- c("TBX21", "EMOMES", "GATA3", "RORC", "FOXP3")

DotPlot(Combined, features = c(CD4.ATAC.markers), dot.scale =30, cols = c("grey","blue")) + scale_y_discrete(limits = rev)
#rename/merge cluster to celltypes in new metadata column "WNNcelltype_ATAC" by high RNA expression of marker genes above, refer to Supplemental Figure S1D, S1G in Wei et al. 2023, Immunity. 

####manual annotation by surface protein markers + RNA + ATAC####
DefaultAssay(Combined) <- "Antibody"

#refer to Table S2 "cluster annotation" in Wei et al. 2023, Immunity for celltypes associated with markers below
CD4.antibody.markers <- c(
  "CD45RA-TotalA", "CD45RO-TotalA", "CD197-TotalA", 
  "CD38-TotalA", "CD278-TotalA", "HLA-DR-TotalA", "CD279-TotalA",
  "CD152-TotalA", "CD223-TotalA", "CD183-TotalA", 
  "CD194-TotalA", "CD196-TotalA", "CD195-TotalA", "CD185-TotalA",
  "KLRG1-TotalA", "GPR56-TotalA", "CD25-TotalA",
  "CD49d-TotalA", "ITGB7-TotalA", "CD103-TotalA")

DotPlot(Combined, features = c(CD4.antibody.markers), dot.scale = 18, cols = c("grey","blue")) + scale_y_discrete(limits = rev)
#expand cluster naming in "WNNcelltype_RNA" to include additional celltypes identifed by high protein expression of marker genes, in new metadata column "WNNcelltype", refer to Figure 1A, S1I in Wei et al. 2023, Immunity. 

####subset Seurat Object by conditions (viremic, suppressed, uninfected) and re-scale subset####
Combined$Condition <- Combined$orig.ident

Idents(Combined) <- "Condition"
Combined <-RenameIdents(Combined, 'HD1' = 'Uninfected', 'YW1' = 'Suppressed', 'YW2' = 'Suppressed', 'YW3' = 'Suppressed',
                        'YW4' = 'Suppressed', 'YW5' = 'Suppressed', 'YW6' = 'Suppressed', 'YW8' = 'Viremic', 'YW9' = 'Viremic', 'YW10' = 'Viremic')
Combined$Condition <- Idents(Combined)

Uninfected <-subset(Combined, idents = c("Uninfected"))
Viremic <-subset(Combined, idents = c("Viremic"))
Suppressed <-subset(Combined, idents = c("Suppressed"))

DefaultAssay(Uninfected) <- "ATAC"; DefaultAssay(Viremic) <- "ATAC"; DefaultAssay(Suppressed) <- "ATAC"
Uninfected <- RunTFIDF(Uninfected)
Uninfected <- FindTopFeatures(Uninfected, min.cutoff = 20)
Uninfected <- RunSVD(Uninfected)
Viremic <- RunTFIDF(Viremic)
Viremic <- FindTopFeatures(Viremic, min.cutoff = 20)
Viremic <- RunSVD(Viremic)
Suppressed <- RunTFIDF(Suppressed)
Suppressed <- FindTopFeatures(Suppressed, min.cutoff = 20)
Suppressed <- RunSVD(Suppressed)

DefaultAssay(Uninfected) <- "Antibody"; DefaultAssay(Viremic) <- "Antibody"; DefaultAssay(Suppressed) <- "Antibody"
Uninfected <- ScaleData(Uninfected); Viremic <- ScaleData(Viremic); Suppressed <- ScaleData(Suppressed)

DefaultAssay(Uninfected) <- "RNA"; DefaultAssay(Viremic) <- "RNA"; DefaultAssay(Suppressed) <- "RNA"
Uninfected <- ScaleData(Uninfected); Viremic <- ScaleData(Viremic); Suppressed <- ScaleData(Suppressed)

####subset Proliferating cluster and re-scale data####
Idents(Combined) <- "WNNcelltype"
Proliferating <- subset(Combined, idents = c("Proliferating"))

DefaultAssay(Proliferating) <- "ATAC"
Proliferating <- RunTFIDF(Proliferating)
Proliferating <- FindTopFeatures(Proliferating, min.cutoff = 20)
Proliferating <- RunSVD(Proliferating)

DefaultAssay(Proliferating) <- "Antibody"
Proliferating <- ScaleData(Proliferating)

DefaultAssay(Proliferating) <- "RNA"
Proliferating <- ScaleData(Proliferating)

####subset Seurat Objects by HIV-1 infection and re-scale subset####
Idents(Viremic) <- "HIV_RNA_DNA"
Viremic_HIV_RNA_DNA_pos <- subset(Viremic, idents = c("positive_positive"))
Viremic_HIV_RNA_pos <- subset(Viremic, idents = c("positive_negative", "positive_positive"))
Viremic_HIV_DNA_pos <- subset(Viremic, idents = c("negative_positive"))
Viremic_HIV_neg <- subset(Viremic, idents = c("negative_negative"))

Idents(Suppressed) <- "HIV_RNA_DNA"
Suppressed_HIV_RNA_DNA_pos <- subset(Suppressed, idents = c("positive_positive"))
Suppressed_HIV_RNA_pos <- subset(Suppressed, idents = c("positive_negative", "positive_positive"))
Suppressed_HIV_DNA_pos <- subset(Suppressed, idents = c("negative_positive"))
Suppressed_HIV_neg <- subset(Suppressed, idents = c("negative_negative"))

Idents(Combined) <- "HIV_RNA_DNA"
HIVposCells <- subset(Combined, idents = c("positive_positive", "negative_positive", "positive_negative"))

DefaultAssay(Viremic_HIV_RNA_DNA_pos) <- "ATAC"; DefaultAssay(Viremic_HIV_DNA_pos) <- "ATAC"; DefaultAssay(Viremic_HIV_RNA_pos) <- "ATAC"; DefaultAssay(Viremic_HIV_neg) <- "ATAC"
DefaultAssay(Suppressed_HIV_RNA_DNA_pos) <- "ATAC"; DefaultAssay(Suppressed_HIV_DNA_pos) <- "ATAC"; DefaultAssay(Suppressed_HIV_RNA_pos) <- "ATAC"; DefaultAssay(Suppressed_HIV_neg) <- "ATAC"
DefaultAssay(HIVposCells) <- "ATAC"

Viremic_HIV_RNA_DNA_pos <- RunTFIDF(Viremic_HIV_RNA_DNA_pos)
Viremic_HIV_RNA_DNA_pos <- FindTopFeatures(Viremic_HIV_RNA_DNA_pos, min.cutoff = 20)
Viremic_HIV_RNA_DNA_pos <- RunSVD(Viremic_HIV_RNA_DNA_pos)
Viremic_HIV_RNA_pos <- RunTFIDF(Viremic_HIV_RNA_DNA_pos)
Viremic_HIV_RNA_pos <- FindTopFeatures(Viremic_HIV_RNA_DNA_pos, min.cutoff = 20)
Viremic_HIV_RNA_pos <- RunSVD(Viremic_HIV_RNA_DNA_pos)
Viremic_HIV_DNA_pos <- RunTFIDF(Viremic_HIV_DNA_pos)
Viremic_HIV_DNA_pos <- FindTopFeatures(Viremic_HIV_DNA_pos, min.cutoff = 20)
Viremic_HIV_DNA_pos <- RunSVD(Viremic_HIV_DNA_pos)
Viremic_HIV_neg <- RunTFIDF(Viremic_HIV_neg)
Viremic_HIV_neg <- FindTopFeatures(Viremic_HIV_neg, min.cutoff = 20)
Viremic_HIV_neg <- RunSVD(Viremic_HIV_neg)
Suppressed_HIV_RNA_DNA_pos <- RunTFIDF(Suppressed_HIV_RNA_DNA_pos)
Suppressed_HIV_RNA_DNA_pos <- FindTopFeatures(Suppressed_HIV_RNA_DNA_pos, min.cutoff = 20)
Suppressed_HIV_RNA_DNA_pos <- RunSVD(Suppressed_HIV_RNA_DNA_pos)
Suppressed_HIV_RNA_pos <- RunTFIDF(Suppressed_HIV_RNA_pos)
Suppressed_HIV_RNA_pos <- FindTopFeatures(Suppressed_HIV_RNA_pos, min.cutoff = 20)
Suppressed_HIV_RNA_pos <- RunSVD(Suppressed_HIV_RNA_pos)
Suppressed_HIV_DNA_pos <- RunTFIDF(Suppressed_HIV_DNA_pos)
Suppressed_HIV_DNA_pos <- FindTopFeatures(Suppressed_HIV_DNA_pos, min.cutoff = 20)
Suppressed_HIV_DNA_pos <- RunSVD(Suppressed_HIV_DNA_pos)
Suppressed_HIV_neg <- RunTFIDF(Suppressed_HIV_neg)
Suppressed_HIV_neg <- FindTopFeatures(Suppressed_HIV_neg, min.cutoff = 20)
Suppressed_HIV_neg <- RunSVD(Suppressed_HIV_neg)
HIVposCells <- RunTFIDF(HIVposCells)
HIVposCells <- FindTopFeatures(HIVposCells, min.cutoff = 20)
HIVposCells <- RunSVD(HIVposCells)

DefaultAssay(Viremic_HIV_RNA_DNA_pos) <- "Antibody"; DefaultAssay(Viremic_HIV_RNA_pos) <- "Antibody"; DefaultAssay(Viremic_HIV_DNA_pos) <- "Antibody"; DefaultAssay(Viremic_HIV_neg) <- "Antibody"
Viremic_HIV_RNA_DNA_pos <- ScaleData(Viremic_HIV_RNA_DNA_pos); Viremic_HIV_RNA_pos <- ScaleData(Viremic_HIV_RNA_pos); Viremic_HIV_DNA_pos <- ScaleData(Viremic_HIV_DNA_pos); Viremic_HIV_neg <- ScaleData(Viremic_HIV_neg)
DefaultAssay(Suppressed_HIV_RNA_DNA_pos) <- "Antibody"; DefaultAssay(Suppressed_HIV_DNA_pos) <- "Antibody"; DefaultAssay(Suppressed_HIV_neg) <- "Antibody"
Suppressed_HIV_RNA_DNA_pos <- ScaleData(Suppressed_HIV_RNA_DNA_pos); Suppressed_HIV_RNA_pos <- ScaleData(Suppressed_HIV_RNA_pos); Suppressed_HIV_DNA_pos <- ScaleData(Suppressed_HIV_DNA_pos); Suppressed_HIV_neg <- ScaleData(Suppressed_HIV_neg)

DefaultAssay(HIVposCells) <- "Antibody"
HIVposCells <- ScaleData(HIVposCells)

DefaultAssay(Viremic_HIV_RNA_DNA_pos) <- "RNA"; DefaultAssay(Viremic_HIV_RNA_pos) <- "RNA"; DefaultAssay(Viremic_HIV_DNA_pos) <- "RNA"; DefaultAssay(Viremic_HIV_neg) <- "RNA"
Viremic_HIV_RNA_DNA_pos <- ScaleData(Viremic_HIV_RNA_DNA_pos); Viremic_HIV_RNA_pos <- ScaleData(Viremic_HIV_RNA_pos); Viremic_HIV_DNA_pos <- ScaleData(Viremic_HIV_DNA_pos); Viremic_HIV_neg <- ScaleData(Viremic_HIV_neg)
DefaultAssay(Suppressed_HIV_RNA_DNA_pos) <- "RNA"; DefaultAssay(Suppressed_HIV_RNA_pos) <- "RNA"; DefaultAssay(Suppressed_HIV_DNA_pos) <- "RNA"; DefaultAssay(Suppressed_HIV_neg) <- "RNA"
Suppressed_HIV_RNA_DNA_pos <- ScaleData(Suppressed_HIV_RNA_DNA_pos); Suppressed_HIV_RNA_pos <- ScaleData(Suppressed_HIV_RNA_pos); Suppressed_HIV_DNA_pos <- ScaleData(Suppressed_HIV_DNA_pos); Suppressed_HIV_neg <- ScaleData(Suppressed_HIV_neg)

DefaultAssay(HIVposCells) <- "RNA"
HIVposCells <- ScaleData(HIVposCells)


#####RNA integration (batch effect correction) and normalization and UMAP generation for HIVposCells####
HIVposCells <- NormalizeData(HIVposCells, normalization.method = 'LogNormalize') %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>% RunPCA(npcs = 6) %>% RunHarmony("MULTI_ID", plot_convergence = TRUE) %>%
  RunUMAP(reduction = "harmony", reduction.name = "umap.rna", reduction.key = "rnaUMAP_", dims = 1:6, assay = 'RNA') %>%
  FindNeighbors(reduction = "harmony", dims = 1:6) %>%
  FindClusters(resolution = 0.3) %>%
  identity()

Idents(HIVposCells) <- "seurat_clusters"
HIVposCells <- RenameIdents(HIVposCells, '0' = 'AP-1', '1' = 'MT', '3' = 'Cytotoxic', '2' = 'IRF')
HIVposCells$seurat_clusters <- Idents(HIVposCells)
