# Wei-et-al-2023
R Scripts related to the Wei et al 2023 manuscript. 

Please see STAR Methods in Wei et al 2023 manuscript for software and package versions. 

These scripts are minimally commented to be functional and reproducible, they are not tutorials. Please see the Methods section in the manuscript for additional details and non-R related analyses.

Scripts in 'Data processing' describe steps taken to generate R objects, functions and processes that are required for various downstream analyses that appear in the manuscript.

Contents of 'Data processing':
Preprocessing.R: Building Seurat objects, doublet removal, low quality cell filter, determining HIV-1-infected cells in filtered bc matrices, RNA/ATAC/protein data normalization and batch effect removal, celltype annotation, creating intermediate objedcts for downstream analyses.

recall ATAC peaks and building motifs.R - require objects generated by 'Preprocessing.R': recall peaks per celltype using MACS3, add JASPAR2022 TF-binding motifs, computing gene accessibilities, adding chromVAR scores, creating TF footprints, creating intermediate objects for downstream analyses.

build WGCNA.R - require objects generated by 'Preprocessing.R': functions from Kazer et al. 2019, Nature Medicine and from Collora et al. 2022, Immunity to perform WGCNA, gene-gene correlation matrices, WGCNA-identified module gene lists, correlation plot functions and module scores.

build GSEA.R - require objects generated by 'Preprocessing.R': GSEA gene sets, intermediates for GSEA (e.g., gene ranks by differential expression, gene set normalized enrichment scores, leading edge genes, etc.).

All scripts in the 'Figures' folder require processed objects and functions in one or more items in the 'Data processing' folder, see scripts for detailed lists.
Scripts in 'Figures' generate the plots as shown in manuscript main figures and, when relevant, the functions/analyses/data processing steps required to generate these plots.

All plots generated in Supplemental figures also require objects generated by scripts in the 'Data processing' folder. Scripts in 'Figures' folder may be modified to reproduce all supplemental figures. Some analyses (related to supplemental figures) may not be entirely representable by scripts in that repository. If there are particular aspects of the analyses you would like to see that are not here, or have any other questions, please email Yulong.Wei@yale.edu or Ya-Chi.Ho@yale.edu.

Raw reads and initial data (e.g., cellranger filtered feature bc matrices, hashtag keys, fragment.tsv etc.) are hosted on GEO:GSE239916

Comment on Seurat v5 and Signac for Seurat v5 (2023-08-09):
Note that some functions have been removed/replaced, default parameters changed, and data structure different, in Seurat v5 vs Seurat v4. Version difference may lead to some scripts breaking and may lead to slight graphical differences (e.g., volcanoplot because of change in FC in FindMarkers), but do not otherwise affect the results (e.g., list of differentially expressed genes).


