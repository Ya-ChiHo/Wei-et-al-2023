library(WGCNA)
library(flashClust)
library(Hmisc)
library(dplyr)
library(doParallel)
registerDoParallel(cores=10)
library(openxlsx)
library(ggplot2)
library(cowplot)
library(Seurat)
library(qlcMatrix)
library(reshape2)
library(scales)
library(ggplot2)
library(Matrix)
library(patchwork)

####required WGCNA functions from Kazer et al. 2019, Nature Medicine####
# Choosing the appropriate power for generating the adjacency matrix
FindPower <- function(datExpr){
  #choose soft-threshold power
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  sft=pickSoftThreshold(datExpr,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,
                        corOptions = list(use = 'p', method = "pearson"),networkType = "signed")
  
  # Plot the results
  par(mfrow = c(1,2));
  cex1 = 0.9;
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
  
  # Red line corresponds to using an R^2 cut-off
  abline(h=0.80,col="red")
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}
# Generating the adjacency matrix and performing clustering
ClusterTOM <- function(datExpr, softPower){
  #dev.off()
  #Calclute the adjacency matrix
  adj= adjacency(datExpr,type = "signed", power = softPower);
  
  #Turn adjacency matrix into a topological overlap matrix (TOM) to minimize the effects of noise and spurious associations.
  TOM=TOMsimilarityFromExpr(datExpr,networkType = "signed", TOMType = "signed", power = softPower, corType="bicor");
  
  colnames(TOM) = rownames(TOM) = colnames(datExpr)
  dissTOM=1-TOM
  
  #Hierarchical clustering of the genes based on the TOM dissimilarity measure
  geneTree = flashClust(as.dist(dissTOM),method="complete");
  
  #Plot the resulting clustering tree (dendrogram)
  plot(geneTree, xlab="", sub="",cex=0.3);
  
  return(list(dissTOM = dissTOM, geneTree = geneTree)) #returns list with dissimilarity TOM, and the clustered gene tree.
}
# Cut the resulting clustering dendrogram using the "tree" method for cutreeDynamic. Minimum module size can be specified.
CutTOMTree <- function(datExpr, dissTOM, geneTree, minModuleSize = 10){
  #dev.off()
  # Module identification using dynamic tree cut, you can also choose the hybrid method
  dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
  #dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
  
  #Get the module labels and the size of each module. Lable 0 is reserved for unassigned genes
  print(table(dynamicMods))
  
  #Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
  
  #Set the diagonal of the dissimilarity to NA 
  diag(dissTOM) = NA;
  
  #extract modules
  module_colors= setdiff(unique(dynamicColors), "grey")
  modules = lapply(module_colors, function(x){colnames(datExpr)[which(dynamicColors==x)]})
  names(modules) = module_colors
  return(list(dyanmicColors = dynamicColors, modules = modules)) #returns list with module colors, and the modules themselves
}
# Merge modules with low dissimilarity. Cutoff for dissimilarity merge can be specified
MergeSimilarModules <- function(datExpr, dynamicColors, geneTree, MEDissThres = 0.5){
  #cacluate eigengenes
  MEList = moduleEigengenes(datExpr, colors=dynamicColors)
  MEs = MEList$eigengenes
  
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs);
  
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");
  
  # Plot the result
  #sizeGrWindow(7, 6)
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  abline(h = MEDissThres, lwd=2, col="red")
  
  # Call an automatic merging function
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors;
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs;
  #plot showing how merged modules exist
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
  #extract merged modules
  merged_module_colors= setdiff(unique(mergedColors), "grey")
  merged_modules = lapply(merged_module_colors, function(x){colnames(datExpr)[which(mergedColors==x)]})
  names(merged_modules) = merged_module_colors
  
  return(list(mergedColors = mergedColors, merged_modules = merged_modules)) #returns list with merged colors, and the merged modules themselves
}
# Test to determine if the genes within the module are truly the least dissimilar compared to randomly generated modules of the same size.
TestModuleSignificance <- function(mod, dissTOM, expr.data, n_perm = 10000, pval = 0.05, n.bin = 10){
  #vectorize the actual distribution of (dis)similarities, and remove zeros!
  true.diss = as.vector(dissTOM[mod,mod])
  true.diss = true.diss[-which(true.diss == 0)]
  
  #size of module for permutations
  mod.size = length(mod)
  
  #bin all genes by expression
  expr.avg = rowMeans(expr.data)
  expr.avg = expr.avg[order(expr.avg)]
  expr.avg.cut = as.numeric(x = cut2(x = expr.avg, m=round(length(expr.avg)/n.bin)))
  names(expr.avg.cut) = names(expr.avg)
  
  #create a table of binnings of all genes and of our module genes
  all.bin.table = table(expr.avg.cut)
  mod.bin.table = table(expr.avg.cut[mod])
  
  #randomly generate module with same expression binning structure and run t.test and record results
  test.results = data.frame(statistic = rep(NA, n_perm), pval = rep(NA, n_perm)) #container for results
  
  for (i in 1:n_perm){ #by permutation
    random.mod = list() #create an empty list we will fill with gene names (each element will be gene names by bin)
    
    #for each bin of the mod bin table, randomly select that number of genes from the full set with that expression
    for (j in names(mod.bin.table)){ 
      bin.genes = sample(names(expr.avg.cut)[which(expr.avg.cut == as.numeric(j))], mod.bin.table[j], replace = FALSE)
      random.mod[[as.numeric(j)]] = bin.genes #stick those genes into the random.mod list
    }
    #unlist and vectorize the distribution of (dis)similarities (remember to remove zeros)
    random.mod = unlist(random.mod)
    random.diss = as.vector(dissTOM[random.mod,random.mod])
    random.diss = random.diss[-which(random.diss == 0)]
    
    #perform a one-sided wilcox.test and record the statistic and p-value.
    #Note, IMPORTANT: here we perform the test asking if the true diss is LESS THAN the random diss, as we are trying to minimize dissimilarity
    test = wilcox.test(x = true.diss, y = random.diss, alternative = "less")
    test.results[i,] = c(test$statistic, test$p.value)
  }
  
  #correct for multiple hypothesis testing, and then report the proportion of bad tests
  test.results$FDR = p.adjust(test.results$pval, method = "fdr")
  num.failed.tests = sum(test.results$FDR > pval)
  print(paste(paste(num.failed.tests, n_perm, sep="/"), "permutations failed the Mann-Whitney test.", sep=" "))
  
  #is the percentage of failed tests less than or equal to the p-val specified?
  return(num.failed.tests/n_perm <= pval) #returns a vector of booleans indicating if each module was significant based on the specific p.val
}

####wrapper function for above functions with prompts, from Collora et al. 2022, Immunity####
modid<-function(seuratobj, cluster, prefix="modules"){
  #initial ensuring we're scaled and PCA ready
  test_clus<-subset(seuratobj, idents=cluster)
  test_clus<-FindVariableFeatures(test_clus)
  test_clus<-ScaleData(test_clus)
  test_clus<-RunPCA(test_clus, reduction.key = "PCWGCNA_")
  print(ElbowPlot(test_clus))
  nPCS<-readline("how many NPCs?")
  nPCS<-as.integer(nPCS)
  #expanding the dim predicted list
  test_clus<-ProjectDim(test_clus)
  
  genes<-c()
  for (i in 1:nPCS){genelist<-TopFeatures(test_clus, dim = i, nfeatures = 50,balanced = T, projected = T )
  for (i in 1:length(genelist)){genes<-c(genes, genelist[[i]])}}
  genes<-unique(genes)
  test_clus<-as.matrix(test_clus@assays$RNA@data[genes,])
  FindPower(datExpr=t(test_clus))
  softpower<-readline("what softpower?")
  softpower<-as.integer(softpower)
  #run all of sams functions
  test_clus_tom<-ClusterTOM(datExpr = t(test_clus), softPower = softpower)
  test_clus_mods<-CutTOMTree(datExpr = t(test_clus), geneTree = test_clus_tom$geneTree, dissTOM = test_clus_tom$dissTOM, minModuleSize = 10 )
  test_clus_merge_mods<-MergeSimilarModules(datExpr=t(test_clus), dynamicColors = test_clus_mods$dyanmicColors, geneTree = test_clus_tom$geneTree, MEDissThres = 0.5)
  print(test_clus_merge_mods$merged_modules)                                          
  test_clus_merge_mods.isSig = sapply(test_clus_merge_mods$merged_modules, function(module){
    TestModuleSignificance(mod = module, dissTOM = test_clus_tom$dissTOM, expr.data = test_clus,
                           n_perm = 10000, pval = 0.05, n.bin = 10)
  })
  test_clus_merge_mods.isSig = test_clus_merge_mods$merged_modules[test_clus_merge_mods.isSig]
  print(test_clus_merge_mods.isSig)
  names(test_clus_merge_mods.isSig)<-paste(prefix,names(test_clus_merge_mods.isSig ), sep="_")
  return(test_clus_merge_mods.isSig)
}

####generate gene-gene expression pearson's correlation coefficient matrices####
#gene-gene correlation matrix per condition (viremia, suppressed, uninfected)
corCombined<-corSparse(t(Combined@assays$RNA@data))
corSuppressed<-corSparse(t(Suppressed@assays$RNA@data))
corViremic<-corSparse(t(Viremic@assays$RNA@data))
corUninfected<-corSparse(t(Uninfected@assays$RNA@data))
rownames(corCombined)<-colnames(corCombined)<-rownames(Combined)
rownames(corSuppressed)<-colnames(corSuppressed)<-rownames(Suppressed)
rownames(corViremic)<-colnames(corViremic)<-rownames(Viremic)
rownames(corUninfected)<-colnames(corUninfected)<-rownames(Uninfected)

corViremic_corHealthy <- list(corViremic, corUninfected)
corViremic_corSuppressed <- list(corViremic, corSuppressed)
corSuppressed_corHealthy <- list(corSuppressed, corUninfected)

#WGCNA modules detected in Viremia
DefaultAssay(Combined) <- "RNA"; Idents(Combined) <- "Condition"
modid(Combined, "Viremic", prefix="modules")
#pc = 7, power = 3

module_viremic_WGCNA_proliferation <- c("MKI67", "RRM2", "ASPM", "STMN1", "PCLAF", "CENPF",
                                        "HIST1H3B", "TYMS", "NUSAP1", "GTSE1", "POLQ", "NCAPG", 
                                        "KIF11", "HIST1H1B", "TPX2", "KNL1", "TOP2A", "CLSPN", 
                                        "DTL", "SHCBP1", "BIRC5", "DIAPH3", "MCM4", "HELLS", "TK1", "TSHZ2", "HJURP",   
                                        "KIF2C", "CIT", "KIF15", "NCAPG2", "SELL", "LIMS1", "ICOS")

module_viremic_WGCNA_cytotoxic <- c("NKG7", "CCL5", "ZEB2", "GZMH", "FGFBP2", "GZMA", "SYNE1", "CCL4", "ADGRG1",
                                    "GZMB", "PRF1",  "C1orf21", "PLEK", "ENC1", "CST7", "IFNG", "GNLY", "HLA-DRB5", 
                                    "HLA-DPB1", "HLA-DRB1", "CD74", "HLA-DRA", "CTSW", "SPON2", "HLA-DPA1", "MIAT", "TGFBR3", "MAF")

#WGCNA modules detected in Suppressed
DefaultAssay(Combined) <- "RNA"; Idents(Combined) <- "Condition"
modid(Combined, "Suppressed", prefix="modules")
#pc = 7, power = 3

module_suppressed_WGCNA_AP1 <- c("FOS", "FOSB", "JUN", "JUNB", "TNFAIP3", "ZFP36L2", 
                                 "ZFP36", "NFKBIA", "SLC2A3", "IER2", "RGCC", "NR4A2", "IFRD1",
                                 "IER5", "GADD45B", "LMNA", "DUSP2", "ARIH1", "KLF6", "VIM", 
                                 "FTH1", "TXNIP", "ITGB1", "NAMPT")


#gene-gene correlation matrix by HIV-1-infection in viremia (HIV-1 RNA+ DNA+ cells, HIV-1 RNA- DNA+ cells, HIV-1- cells)
corRNA_doublePositive_viremic <- corSparse(t(Viremic_HIV_RNA_DNA_pos@assays$RNA@data))
corDNApositive_viremic <- corSparse(t(Viremic_HIV_DNA_pos@assays$RNA@data))
corHIVnegative_viremic <- corSparse(t(Viremic_HIV_neg@assays$RNA@data))
rownames(corRNA_doublePositive_viremic)<-rownames(Viremic_HIV_RNA_DNA_pos)
colnames(corRNA_doublePositive_viremic) <- rownames(Viremic_HIV_RNA_DNA_pos)
rownames(corDNApositive_viremic)<-rownames(Viremic_HIV_DNA_pos)
colnames(corDNApositive_viremic) <- rownames(Viremic_HIV_DNA_pos)
rownames(corHIVnegative_viremic)<-rownames(Viremic_HIV_neg)
colnames(corHIVnegative_viremic) <- rownames(Viremic_HIV_neg)

corRNA_doublePositive_viremic_corDNApositive_viremic <- list(corRNA_doublePositive_viremic, corDNApositive_viremic)
corRNA_doublePositive_viremic_corHIVnegative_viremic <- list(corRNA_doublePositive_viremic, corHIVnegative_viremic)
corDNApositive_viremic_corHIVnegative_viremic <- list(corDNApositive_viremic, corHIVnegative_viremic)

#WGCNA modules detected in Viremic HIV-1 RNA+ cells
Idents(Viremic) <- "HIV_RNA_DNA2"; DefaultAssay(Viremic) <- "RNA"

modid(Viremic, "RNA and double positive", prefix="modules")
#pc = 8, power = 3

module_RNApos_cytotoxic <- c("CCL5","HLA-A","NKG7", "ZEB2", "FGFBP2", "GZMH", "GZMA", "IFNG", "HSPA5","SYNE1")

module_RNApos_proliferation <- c("MKI67", "ASPM", "CDK1", "HJURP", "BIRC5", "BUB1", "CENPF", "TPX2", "ANLN",  
                                  "GTSE1", "NCAPG", "PRR11", "KIF23", "CENPU", "GTF2H2C", "CIT", "RRM2", "BRCA2", "ECT2",
                                  "TOP2A", "PCLAF", "ESPL1", "IGHM", "CD84",  "ARHGEF35", "HIST1H4H", "IL21",
                                  "GNB4", "SGO2", "GATD1", "SPG11", "ANAPC2", "IKZF3", "YOD1", "GEN1", "BAK1", "MIR181A1HG",
                                  "SOX4", "MINPP1", "TSEN34", "GLRX3", "MRPL51", "KIF11", "C21orf58", "RC3H2", "MICB", "CD74",
                                  "HMGB3", "CCNB1", "SPATS2L", "NCSTN")


#gene-gene correlation matrix in HIV-1-infected cells
Idents(HIVposCells) <- "seurat_clusters"; DefaultAssay(HIVposCells) <- "RNA"
HIVpos_AP1Mod <- subset(HIVposCells, idents = c("AP1"))
HIVpos_IRFMod <- subset(HIVposCells, idents = c("IRF"))
HIVpos_cytotoxicMod <- subset(HIVposCells, idents = c("cytotoxic"))
HIVpos_MTMod <- subset(HIVposCells, idents = c("MT"))

corHIVpos_AP1Mod <- corSparse(t(HIVpos_AP1Mod@assays$RNA@data))
corHIVpos_IRFMod <- corSparse(t(HIVpos_IRFMod@assays$RNA@data))
corHIVpos_cytotoxicMod <- corSparse(t(HIVpos_cytotoxicMod@assays$RNA@data))
corHIVpos_MTMod <- corSparse(t(HIVpos_MTMod@assays$RNA@data))

rownames(corHIVpos_AP1Mod)<-rownames(HIVpos_AP1Mod)
colnames(corHIVpos_AP1Mod) <- rownames(HIVpos_AP1Mod)
rownames(corHIVpos_IRFMod)<-rownames(HIVpos_IRFMod)
colnames(corHIVpos_IRFMod) <- rownames(HIVpos_IRFMod)
rownames(corHIVpos_cytotoxicMod)<-rownames(HIVpos_cytotoxicMod)
colnames(corHIVpos_cytotoxicMod) <- rownames(HIVpos_cytotoxicMod)
rownames(corHIVpos_MTMod)<-rownames(HIVpos_MTMod)
colnames(corHIVpos_MTMod) <- rownames(HIVpos_MTMod)

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

corHIVpos_AP1Mod[is.nan(corHIVpos_AP1Mod)] <- 0
corHIVpos_IRFMod[is.nan(corHIVpos_IRFMod)] <- 0
corHIVpos_cytotoxicMod[is.nan(corHIVpos_cytotoxicMod)] <- 0
corHIVpos_MTMod[is.nan(corHIVpos_MTMod)] <- 0

corHIVpos_AP1Mod_corHIVpos_MTMod <- list(corHIVpos_AP1Mod, corHIVpos_MTMod)
corHIVpos_cytotoxicMod_corHIVpos_apoptoticMod <- list(corHIVpos_cytotoxicMod, corHIVpos_MTMod)
corHIVpos_IRFMod_corHIVpos_MTMod <- list(corHIVpos_IRFMod, corHIVpos_MTMod)

#WGCNA modules detected in all HIV-1+ cells
DefaultAssay(Combined) <- "RNA"; 
Combined$HIVpos <- Combined$HIV_RNA_DNA
Idents(Combined) <- "HIVpos"
Combined <-RenameIdents(Combined, 'posivie_positive' = 'positive', 'positive_negative' = 'positive', 'negative_positive' = 'positive', 'negative_negative' = 'negative') 
Combined$HIVpos <- Idents(Combined)

modid(Combined, "positive", prefix="modules")
#pc = 8, power = 3

HIVpos_AP1 <- c("FOSB"   ,    "JUNB"    ,   "FOS"   ,     "JUN"  ,      "NFKBIA",     "TNFAIP3"   , "AC016831.1", "AC020916.1" ,"IER2"  ,     "VIM"  ,      "PTMA"   ,   
                "AC016831.7", "IRF1"    ,   "NAMPT"   ,   "EIF4A3")

HIVpos_cytotoxic <- c("SYNE1" ,  "NKG7" ,   "ZEB2" ,   "GZMH"   , "FGFBP2" , "IFNG"  ,  "ENC1"  ,  "GZMB"   , "MRPL27",  "SLC16A6" ,"DUSP4" ,  "RGS1" ,   "ACOT13")


HIVpos_IRF <- c( "GTSE1"  ,    "MKI67"   ,   "TPX2"     ,  "ASPM" ,      "CDK1"    ,   "CENPF",      "BUB1"      , "KIF23",      "HJURP"     , "BRCA2",      "KIF15"     ,
                           "ECT2"    ,   "CKAP2L"   ,  "STIL"      , "DIAPH3",     "ANLN"   ,    "BIRC5"   ,   "RRM2"      , "TOP2A" ,     "NCAPG"    ,  "CENPU" ,     "DTL"       ,
                           "PRR11"   ,   "PCLAF"    ,  "KIF14"   ,   "PLK1"  ,     "MTFR2"     , "CDC20",     "KIF2C"     , "UBE2C"  ,    "ACTG1"      ,"ARPC1B" ,    "LIMS1"     ,
                           "PPIB"     ,  "ICOS"     ,  "SLC39A10"  , "POLQ"   ,    "MCM10"    ,  "MCM4"  ,     "HIST1H2BH"  ,"PBK"    ,    "HIST1H1B" ,  "BRCA1"  ,    "ESCO2"     ,
                           "BRIP1"     , "FBXO5"     , "ORC1"     ,  "CLSPN"   ,   "UHRF1"   ,   "HIST1H3B",   "FANCI"     , "UST"     ,   "CENPP"   ,   "SLBP"    ,   "ADGRF3"    ,
                           "DCTPP1"  ,   "HIST1H2BF" , "ATAD5"    ,  "TTK"      ,  "IL21"   ,    "CCNB2"    ,  "IGHM"     ,  "CCNB1"    ,  "CD70"  ,     "GNB4"     ,  "AL050341.2",
                           "LINC00996",  "SPC24"     , "HIST1H4H",   "WDR62"     , "NEMP1" ,     "DOT1L"     , "ARHGEF35",   "CCDC150"   , "BMERB1" ,    "TROAP"     , "TSEN34" ,   
                           "CMTM7"     , "AC004889.1", "CD59"  ,     "CCNA2"      ,"MAGOH",      "RMI2"       ,"ESPL1",      "PTTG1"      ,"KIF4A",      "TYMS"       ,"DLGAP5",    
                           "SHCBP1"    , "CES2"  )  

####function for split condition pearson's correlation coefficient heatmap matrix for WGCNA modules####
#modified from Collora et al. Immunity, 2022
splitcorheatmap<-function(cor, genes, prefix, n_ext=1,exp=TRUE){
  samples<-names(cor)
  plot<-c()
  if(exp==TRUE){for(i in 1:length(genes)){
    plot<-unique(c(plot,extract.cor(cor,genes[[i]], ngene = n_ext)))
  }
    plot2<-c()
    for(i in 1:length(plot)){
      plot2<-c(plot2, plot[[i]])
    }
    plot2<-unique(plot2)}else{plot2<-genes}
  
  corx<-cor[[1]][plot2,plot2]
  cory<-cor[[2]][plot2,plot2]
  corx[lower.tri(corx)]<-cory[lower.tri(cory)]
  x<-setNames(melt(corx), c("x","y","weight"))
  x[is.na(x)]<-0
  x$x<-factor(x$x, levels = plot2)
  x$y<-factor(x$y, levels = plot2)
  p<-ggplot(x, aes(x=x,y, fill=weight))+geom_tile()+
    scale_fill_gradient2(low = "blue",high = "red", limits=c(-1, 1), oob=squish)+
    theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.4))+theme(axis.ticks = element_blank())
  savename<-paste(prefix,"_",samples[1],"_",samples[2],"split.pdf", sep = "")
  ggsave(filename = savename, plot = p, device = "pdf", width = 10, height = 9)
}

####WGCNA modules add modulescore####
DefaultAssay(Combined) <- "RNA"
Combined <-AddModuleScore(Combined, features = list(module_viremic_WGCNA_proliferation), name = 'module_viremic_WGCNA_proliferation')
Combined <-AddModuleScore(Combined, features = list(module_viremic_WGCNA_cytotoxic), name = 'module_viremic_WGCNA_cytotoxic')
Combined <-AddModuleScore(Combined, features = list(module_suppressed_WGCNA_AP1), name = 'module_suppressed_WGCNA_AP1')

DefaultAssay(Viremic) <- "RNA"
Viremic <-AddModuleScore(Viremic1, features = list(module_RNApos_proliferation), name = 'module_RNApos_proliferation')

