library(Seurat)
library(Signac)
library(rlang)
library(tidyverse)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(GenomeInfoDb)
library(ggplot2)
library(cowplot)
library(dplyr)
library(JASPAR2022)
library(TFBSTools)
library(chromVAR)
library(stringr)
library(BiocParallel)
library(chromVARmotifs)
library(MACSr)
library(basilisk)


####rebuild ATAC assay with MACS3 recalled peaks by celltype####
DefaultAssay(Combined) <- "ATAC"
frags <- Fragments(Combined)
macs2peak <- CallPeaks(Combined, group.by = "WNNcelltypeRNA", macs3.path = "~/MyPythonEnv/bin/macs3",outdir = '~/macs3peaks')

plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2)

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

HD1.counts <- FeatureMatrix(fragments = HD1.frags, features = macs3peak, cells = rownames(HD1.md))
YW1.counts <- FeatureMatrix(fragments = YW1.frags, features = macs3peak, cells = rownames(YW1.md))
YW2.counts <- FeatureMatrix(fragments = YW2.frags, features = macs3peak, cells = rownames(YW2.md))
YW3.counts <- FeatureMatrix(fragments = YW3.frags, features = macs3peak, cells = rownames(YW3.md))
YW4.counts <- FeatureMatrix(fragments = YW4.frags, features = macs3peak, cells = rownames(YW4.md))
YW5.counts <- FeatureMatrix(fragments = YW5.frags, features = macs3peak, cells = rownames(YW5.md))
YW6.counts <- FeatureMatrix(fragments = YW6.frags, features = macs3peak, cells = rownames(YW6.md))
YW8.counts <- FeatureMatrix(fragments = YW8.frags, features = macs3peak, cells = rownames(YW8.md))
YW9.counts <- FeatureMatrix(fragments = YW9.frags, features = macs3peak, cells = rownames(YW9.md))
YW10.counts <- FeatureMatrix(fragments = YW10.frags, features = macs3peak, cells = rownames(YW10.md))

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

Idents(Combined.ATAC) <- "is_cell"
Combined.ATAC <-subset(Combined.ATAC, idents = c("1"))

Combined.ATAC <- RunTFIDF(Combined.ATAC)
Combined.ATAC <- FindTopFeatures(Combined.ATAC, min.cutoff = 'q0')
Combined.ATAC <- RunSVD(Combined.ATAC)

####transfer cell bc info from processed Combined Seurat Object to new Combined.ATAC object####
Combined.ATAC$gex_barcode <- colnames(Combined.ATAC)
all.equal(colnames(Combined.ATAC), colnames(Combined))
joint.bcs<- intersect(colnames(Combined.ATAC), colnames(Combined))
Combined.ATAC <-subset(Combined.ATAC, idents = joint.bcs)

Combined.ATAC$Condition <- Combined$Condition
Combined.ATAC$MULTI_classification <- Combined$MULTI_classification
Combined.ATAC$WNNcelltype <- Combined$WNNcelltype
Combined.ATAC$WNNcelltype_ATAC <- Combined$WNNcelltype_ATAC
Combined.ATAC$HIV_DNA <- Combined$HIV_DNA
Combined.ATAC$HIV_RNA <- Combined$HIV_RNA
Combined.ATAC$HIV_RNA_DNA <- Combined$HIV_RNA_DNA

####Add JASPAR2022 motif info####
pfm <- getMatrixSet(x = JASPAR2022, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- as.logical(seqnames(granges(Combined.ATAC)) %in% main.chroms)
Combined.ATAC <- Combined.ATAC[keep.peaks,]
Combined.ATAC <- AddMotifs(Combined.ATAC, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
motif.names <- Combined.ATAC@assays$ATAC@motifs@motif.names

####rebuild gene accessibility####
gene.activities <- GeneActivity(Combined.ATAC)
Combined.ATAC[['activities']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(Combined.ATAC) <- "activities"
Combined.ATAC <- NormalizeData(object = Combined.ATAC, assay = 'activities', 
                               normalization.method = 'LogNormalize', 
                               scale.factor = median(Combined.ATAC$nCount_activities))

####add chromVAR scores####
register(SerialParam(stop.on.error = TRUE, log = TRUE,
                     threshold = "INFO", logdir = NA_character_, progressbar = TRUE))
DefaultAssay(Combined.ATAC) <- "ATAC"
Combined.ATAC <- RunChromVAR(Combined.ATAC, genome = BSgenome.Hsapiens.UCSC.hg38)

####annotate peaks to gene by HOMER####
table.macs3peak <- data.frame(seqnames = macs3peak@seqnames, macs3peak@ranges)
table.macs3peak <- data.frame(rownames(table.macs3peak), seqnames = macs3peak@seqnames, macs3peak@ranges, macs3peak@strand)
table.macs3peak <- subset(table.macs3peak, select = -width)
table.macs3peak$macs3peak.strand <- 0
colnames(table.macs3peak) = c("ID", "chromosome", "start", "end", "strand")
write.table(table.macs3peak,"~/table_macs3peak.txt",sep="\t", row.names=FALSE, quote = FALSE)

#on shell: annotatePeaks.pl ~/table_macs3peak.txt hg38 > ~/homer_macs3peaks.txt

homer_annotation <- read.delim(file = "~/homer_macs3peaks.txt", header = TRUE, sep = "\t")
homer_annotation$coord <- str_c(homer_annotation$Chr,'-', homer_annotation$Start, '-', homer_annotation$End)

####subset ATAC Object by conditions (viremic, suppressed, uninfected) and re-scale subset####
Combined.ATAC$Condition <- Combined.ATAC$orig.ident

Uninfected.ATAC <-subset(Combined.ATAC, idents = c("Uninfected"))
Viremic.ATAC <-subset(Combined.ATAC, idents = c("Viremic"))
Suppressed.ATAC <-subset(Combined.ATAC, idents = c("Suppressed"))

DefaultAssay(Uninfected.ATAC) <- "ATAC"; DefaultAssay(Viremic.ATAC) <- "ATAC"; DefaultAssay(Suppressed.ATAC) <- "ATAC"
Uninfected.ATAC <- RunTFIDF(Uninfected)
Uninfected.ATAC <- FindTopFeatures(Uninfected, min.cutoff = 20)
Uninfected.ATAC <- RunSVD(Uninfected)
Viremic.ATAC <- RunTFIDF(Viremic)
Viremic.ATAC <- FindTopFeatures(Viremic, min.cutoff = 20)
Viremic.ATAC <- RunSVD(Viremic)
Suppressed.ATAC <- RunTFIDF(Suppressed)
Suppressed.ATAC <- FindTopFeatures(Suppressed, min.cutoff = 20)
Suppressed.ATAC <- RunSVD(Suppressed)

####subset ATAC Proliferating cluster and re-scale data####
Idents(Combined.ATAC) <- "WNNcelltype"
Proliferating.ATAC <- subset(Combined.ATAC, idents = c("Proliferating"))

DefaultAssay(Proliferating.ATAC) <- "ATAC"
Proliferating.ATAC <- RunTFIDF(Proliferating.ATAC)
Proliferating.ATAC <- FindTopFeatures(Proliferating.ATAC, min.cutoff = 20)
Proliferating.ATAC <- RunSVD(Proliferating.ATAC)

####subset ATAC Objects by HIV-1 infection and re-scale subset####
Idents(Viremic.ATAC) <- "HIV_RNA_DNA"
Viremic_HIV_RNA_DNA_pos.ATAC <- subset(Viremic.ATAC, idents = c("positive_positive"))
Viremic_HIV_RNA_pos.ATAC <- subset(Viremic.ATAC, idents = c("positive_negative"))
Viremic_HIV_DNA_pos.ATAC <- subset(Viremic.ATAC, idents = c("negative_positive"))
Viremic_HIV_neg.ATAC <- subset(Viremic.ATAC, idents = c("negative_negative"))

Idents(Suppressed.ATAC) <- "HIV_RNA_DNA"
Suppressed_HIV_RNA_DNA_pos.ATAC <- subset(Suppressed.ATAC, idents = c("positive_positive"))
Suppressed_HIV_RNA_pos.ATAC <- subset(Suppressed.ATAC, idents = c("positive_negative"))
Suppressed_HIV_DNA_pos.ATAC <- subset(Suppressed.ATAC, idents = c("negative_positive"))
Suppressed_HIV_neg.ATAC <- subset(Suppressed.ATAC, idents = c("negative_negative"))

Idents(Combined.ATAC) <- "HIV_RNA_DNA"
HIVposCells.ATAC <- subset(Combined.ATAC, idents = c("positive_positive", "negative_positive", "positive_negative"))

Viremic_HIV_RNA_DNA_pos.ATAC <- RunTFIDF(Viremic_HIV_RNA_DNA_pos.ATAC)
Viremic_HIV_RNA_DNA_pos.ATAC <- FindTopFeatures(Viremic_HIV_RNA_DNA_pos.ATAC, min.cutoff = 20)
Viremic_HIV_RNA_DNA_pos.ATAC <- RunSVD(Viremic_HIV_RNA_DNA_pos.ATAC)
Viremic_HIV_RNA_pos.ATAC <- RunTFIDF(Viremic_HIV_RNA_pos.ATAC)
Viremic_HIV_RNA_pos.ATAC <- FindTopFeatures(Viremic_HIV_RNA_pos.ATAC, min.cutoff = 20)
Viremic_HIV_RNA_pos.ATAC <- RunSVD(Viremic_HIV_RNA_pos.ATAC)
Viremic_HIV_DNA_pos.ATAC <- RunTFIDF(Viremic_HIV_DNA_pos.ATAC)
Viremic_HIV_DNA_pos.ATAC <- FindTopFeatures(Viremic_HIV_DNA_pos.ATAC, min.cutoff = 20)
Viremic_HIV_DNA_pos.ATAC <- RunSVD(Viremic_HIV_DNA_pos.ATAC)
Viremic_HIV_neg.ATAC <- RunTFIDF(Viremic_HIV_neg.ATAC)
Viremic_HIV_neg.ATAC <- FindTopFeatures(Viremic_HIV_neg.ATAC, min.cutoff = 20)
Viremic_HIV_neg.ATAC <- RunSVD(Viremic_HIV_neg.ATAC)
HIVposCells.ATAC <- RunTFIDF(HIVposCells.ATAC)
HIVposCells.ATAC <- FindTopFeatures(HIVposCells.ATAC, min.cutoff = 20)
HIVposCells.ATAC <- RunSVD(HIVposCells.ATAC)

####TF footprints####
Idents(Combined.ATAC) <- "WNNcelltype"
Combined.ATAC <- Footprint(Combined.ATAC, motif.name = c("IRF4", "IRF7", "STAT1::STAT2"), genome = BSgenome.Hsapiens.UCSC.hg38, in.peaks = TRUE)


