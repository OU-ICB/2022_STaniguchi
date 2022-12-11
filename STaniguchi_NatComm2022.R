#Here we show de novo data processing of public human lung single-cell data using Seurat, R analytical tools for single cell genomics. 
#Data were downloaded from the NCBI Gene Expression Omnibus (GEO) with accession code GSE154826 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154826).

#Install packages
install.packages('Seurat')
library(Seurat)
install.packages('Matrix')
library(Matrix)
install.packages("pheatmap")
install.packages("RColorBrewer")
install.packages("ggplot2")
install.packages("gplots")
install.packages("calibrate")
install.packages("maptools")
install.packages("FactoMineR")
install.packages("rgl")
install.packages("caret")
install.packages("BBmisc")
install.packages("survival")
install.packages("survminer")
install.packages("tidyverse")
library(tidyverse)
install.packages("dplyr")
library(dplyr)
install.packages("data.table")
install.packages("reshape2")
library(reshape2)
install.packages("beeswarm")
install.packages("devtools")
install.packages('cowplot')

#Load the public data & check data quality
#Sample data we tested are as follows

#Normal group
#1, ID:1_370
#2, ID:7_378
#3, ID:10_403
#4, ID:14_406
#5, ID:16_408
#6, ID:18_410
#7, ID:28_514
#8, ID:37_564
#9, ID:41_570
#10, ID:44_572
#11, ID:46_578
#12, ID:49_581
#13, ID:307_190329

#Tumor group
#1, ID:2_370
#2, ID:8_378
#3, ID:11_403
#4, ID:15_406
#5, ID:17_408
#6, ID:19_410
#7, ID:29_514
#8, ID:38_564
#9, ID:42_570
#10, ID:45_572
#11, ID:47_578
#12, ID:50_581
#13, ID:308_190329

#Verified these data in turn
matrix_dir = "./Desktop/Folder/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv")
features.path <- paste0(matrix_dir, "features.tsv")
matrix.path <- paste0(matrix_dir, "matrix.mtx")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
cellid <- paste("Cell",1:length(barcode.names$V1))
colnames(mat) = cellid
rownames(mat) = feature.names$V2
data <- CreateSeuratObject(counts = mat,  project = "Quality Check", min.cells = 5)

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

#Remove low quality data
data.a <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 100)

#Extract high quality data
data.b <- subset(data.a, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)

#Count the number of high quality cells
#Calculate %(High Quality/All data)
length(data2$percent.mt)/length(data$percent.mt)*100

#Resluts of quality check
#Normal group
#1, #high quality cell -> 3445, %(high/all) -> 0.467258
#2, #high quality cell -> 22538, %(high/all) -> 3.056912
#3, #high quality cell -> 6184,  %(high/all) -> 0.8387587***
#4, #high quality cell -> 2982,  %(high/all) -> 0.4044596
#5, #high quality cell -> 4206,  %(high/all) -> 0.5704753
#6, #high quality cell -> 5355,  %(high/all) -> 0.7263184**
#7, #high quality cell -> 3363,  %(high/all) -> 0.4561361
#8, #high quality cell -> 8106,  %(high/all) -> 1.099447***
#9, #high quality cell -> 5265,  %(high/all) -> 0.7141113**
#10, #high quality cell -> 4962,  %(high/all) -> 0.6730143*
#11, #high quality cell -> 3940,  %(high/all) -> 0.5343967
#12, #high quality cell -> 6258,  %(high/all) -> 0.8487956***
#13, #high quality cell -> 4787,  %(high/all) -> 0.6492784*

#Tumor group
#1, #high quality cell -> 2213, %(high/all) -> 0.3001573
#2, #high quality cell -> 5922, %(high/all) -> 0.8032227**
#3, #high quality cell -> 3928, %(high/all) -> 0.5327691
#4, #high quality cell -> 4238, %(high/all) -> 0.5748155
#5, #high quality cell -> 5817, %(high/all) -> 0.7889811*
#6, #high quality cell -> 6250, %(high/all) -> 0.8477105***
#7, #high quality cell -> 3425, %(high/all) -> 0.4645454
#8, #high quality cell -> 12282, %(high/all) -> 1.665853
#9, #high quality cell -> 2746, %(high/all) -> 0.3724501
#10, #high quality cell -> 5925, %(high/all) -> 0.8036296**
#11, #high quality cell -> 6965, %(high/all) -> 0.9446886***
#12, #high quality cell -> 5375, %(high/all) -> 0.729031*
#13, #high quality cell -> 5958, %(high/all) -> 0.8081055***


#Select five data used in subsequent analysis
#Convert data file to SeuratObject
#Normal
normal1 <- CreateSeuratObject(counts = normal.mat, project = "Normal#3", min.cells = 5)
normal2 <- CreateSeuratObject(counts = normal.mat, project = "Normal#8", min.cells = 5)
normal3 <- CreateSeuratObject(counts = normal.mat, project = "Normal#12", min.cells = 5)
normal4 <- CreateSeuratObject(counts = normal.mat, project = "Normal#9", min.cells = 5)
normal5 <- CreateSeuratObject(counts = normal.mat, project = "Normal#6", min.cells = 5)

#LUAD(Tumor)
Tumor1 <- CreateSeuratObject(counts = Tumor.mat, project = "Tumor#6", min.cells = 5)
Tumor2 <- CreateSeuratObject(counts = Tumor.mat, project = "Tumor#11", min.cells = 5)
Tumor3 <- CreateSeuratObject(counts = Tumor.mat, project = "Tumor#13", min.cells = 5)
Tumor4 <- CreateSeuratObject(counts = Tumor.mat, project = "Tumor#10", min.cells = 5)
Tumor5 <- CreateSeuratObject(counts = Tumor.mat, project = "Tumor#2", min.cells = 5)


#Make meta.data
normal1@meta.data$Tissue <- "Normal"
normal2@meta.data$Tissue <- "Normal"
normal3@meta.data$Tissue <- "Normal"
normal4@meta.data$Tissue <- "Normal"
normal5@meta.data$Tissue <- "Normal"
Tumor1@meta.data$Tissue <- "LUAD"
Tumor2@meta.data$Tissue <- "LUAD"
Tumor3@meta.data$Tissue <- "LUAD"
Tumor4@meta.data$Tissue <- "LUAD"
Tumor5@meta.data$Tissue <- "LUAD"

normal1@meta.data$PatientID <- "10_403"
normal2@meta.data$PatientID <- "37_564"
normal3@meta.data$PatientID <- "49_581"
normal4@meta.data$PatientID <- "41_570"
normal5@meta.data$PatientID <- "18_410"
Tumor1@meta.data$PatientID <- "19_410"
Tumor2@meta.data$PatientID <- "47_578"
Tumor3@meta.data$PatientID <- "308_190329"
Tumor4@meta.data$PatientID <- "45_572"
Tumor5@meta.data$PatientID <- "8_378"

normal1@meta.data$Sample <- "Normal1"
normal2@meta.data$Sample <- "Normal2"
normal3@meta.data$Sample <- "Normal3"
normal4@meta.data$Sample <- "Normal4"
normal5@meta.data$Sample <- "Normal5"
Tumor1@meta.data$Sample <- "LUAD1"
Tumor2@meta.data$Sample <- "LUAD2"
Tumor3@meta.data$Sample <- "LUAD3"
Tumor4@meta.data$Sample <- "LUAD4"
Tumor5@meta.data$Sample <- "LUAD5"


#Quality check and Normalize data
normal1[["percent.mt"]] <- PercentageFeatureSet(normal1, pattern = "^MT-")
normal2[["percent.mt"]] <- PercentageFeatureSet(normal2, pattern = "^MT-")
normal3[["percent.mt"]] <- PercentageFeatureSet(normal3, pattern = "^MT-")
normal4[["percent.mt"]] <- PercentageFeatureSet(normal4, pattern = "^MT-")
normal5[["percent.mt"]] <- PercentageFeatureSet(normal5, pattern = "^MT-")
Tumor1[["percent.mt"]] <- PercentageFeatureSet(Tumor1, pattern = "^MT-")
Tumor2[["percent.mt"]] <- PercentageFeatureSet(Tumor2, pattern = "^MT-")
Tumor3[["percent.mt"]] <- PercentageFeatureSet(Tumor3, pattern = "^MT-")
Tumor4[["percent.mt"]] <- PercentageFeatureSet(Tumor4, pattern = "^MT-")
Tumor5[["percent.mt"]] <- PercentageFeatureSet(Tumor5, pattern = "^MT-")

#Extract good cell data
normal1 <- subset(normal1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20) 
normal2 <- subset(normal2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20) 
normal3 <- subset(normal3, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
normal4 <- subset(normal4, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
normal5 <- subset(normal5, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20) 
Tumor1 <- subset(Tumor1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20) 
Tumor2 <- subset(Tumor2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20) 
Tumor3 <- subset(Tumor3, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20) 
Tumor4 <- subset(Tumor4, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20) 
Tumor5 <- subset(Tumor5, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)

normal1 <- NormalizeData(normal1, normalization.method = "LogNormalize", scale.factor = 10000)
normal2 <- NormalizeData(normal2, normalization.method = "LogNormalize", scale.factor = 10000)
normal3 <- NormalizeData(normal3, normalization.method = "LogNormalize", scale.factor = 10000)
normal4 <- NormalizeData(normal4, normalization.method = "LogNormalize", scale.factor = 10000)
normal5 <- NormalizeData(normal5, normalization.method = "LogNormalize", scale.factor = 10000)
Tumor1 <- NormalizeData(Tumor1, normalization.method = "LogNormalize", scale.factor = 10000)
Tumor2 <- NormalizeData(Tumor2, normalization.method = "LogNormalize", scale.factor = 10000)
Tumor3 <- NormalizeData(Tumor3, normalization.method = "LogNormalize", scale.factor = 10000)
Tumor4 <- NormalizeData(Tumor4, normalization.method = "LogNormalize", scale.factor = 10000)
Tumor5 <- NormalizeData(Tumor5, normalization.method = "LogNormalize", scale.factor = 10000)

#Integrate datasets
all.list <- NULL
all.list <- list(normal1, normal2, normal3, normal4, normal5, Tumor1, Tumor2, Tumor3, Tumor4, Tumor5)
names(all.list) <- c("Normal1","Normal2", "Normal3", "Normal4", "Normal5", "Tumor1","Tumor2","Tumor3","Tumor4","Tumor5")
all.list 

#Normalize and identify variable features for each dataset independently
all.list <- lapply(X = all.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 500) 
})

#Select features that are repeatedly variable across datasets for integration run PCA on each
#Data set using these features
features <- SelectIntegrationFeatures(object.list = all.list)
all.list <- lapply(X = all.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
all.anchors <- FindIntegrationAnchors(object.list = all.list, anchor.features = features, reduction = "rpca")
all.integrated <- IntegrateData(anchorset = all.anchors)
DefaultAssay(all.integrated) <- "integrated"

#Perform clustring
all.integrated <- ScaleData(all.integrated, verbose = FALSE)
all.integrated <- RunPCA(all.integrated , npcs = 30, verbose = FALSE)
all.integrated <- RunUMAP(all.integrated , reduction = "pca", dims = 1:30)
all.integrated  <- FindNeighbors(all.integrated , reduction = "pca", dims = 1:30)
all.integrated  <- FindClusters(all.integrated , resolution = 0.3) 

#Visualization
DimPlot(all.integrated, reduction = "umap", label = T, label.size = 6, pt.size = 0.1) + NoLegend()
DimPlot(all.integrated, reduction = "umap", label = T, label.size = 0, pt.size = 0.1) + NoLegend()

DimPlot(all.integrated, reduction = "umap", group.by = "Tissue" , pt.size = 0.1)
DimPlot(all.integrated, reduction = "umap", group.by = "Tissue" , split.by = "Tissue", pt.size = 0.1)
DimPlot(all.integrated, reduction = "umap", group.by = "Tissue" , split.by = "PatientID" , pt.size = 0.1)

#Find markers of each cluster
markers <- FindAllMarkers(all.integrated, only.pos = TRUE)
top100.markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n=100, order_by=avg_log2FC)
write.csv(top100.markers, "./human lung_top100 markers.csv")

#Gene expressions
#11
FeaturePlot(all.integrated, features = "SPON2", pt.size = 0.3, min.cutoff = 0) #NK
FeaturePlot(all.integrated, features = "FCGR3A", pt.size = 0.3, min.cutoff = 0) #NK
FeaturePlot(all.integrated, features = "GZMB", pt.size = 0.3, min.cutoff = 0) #NK

#0,1,9
FeaturePlot(all.integrated, features = "CD3D", pt.size = 0.3, min.cutoff = 0) #T
FeaturePlot(all.integrated, features = "TRAC", pt.size = 0.3, min.cutoff = 0) #T

#8
FeaturePlot(all.integrated, features = "MZB1", pt.size = 0.3, min.cutoff = 0) #Plasma
FeaturePlot(all.integrated, features = "XBP1", pt.size = 0.3, min.cutoff = 0) #Plasma

#5
FeaturePlot(all.integrated, features = "MS4A1", pt.size = 0.3, min.cutoff = 0) #B
FeaturePlot(all.integrated, features = "CD79A", pt.size = 0.3, min.cutoff = 0) #B

#7
FeaturePlot(all.integrated, features = "SFTPC", pt.size = 0.3, min.cutoff = 0) #Epithelial cell
FeaturePlot(all.integrated, features = "MUC1", pt.size = 0.3, min.cutoff = 0) #Epithelial cell

#12
FeaturePlot(all.integrated, features = "TPSAB1", pt.size = 0.3, min.cutoff = 0) #Mast cell
FeaturePlot(all.integrated, features = "TPSB2", pt.size = 0.3, min.cutoff = 0) #Mast cell

FeaturePlot(all.integrated, features = "C1orf54", pt.size = 0.3, min.cutoff = 0) #cDC
FeaturePlot(all.integrated, features = "CLEC9A", pt.size = 0.3, min.cutoff = 0) #cDC

FeaturePlot(all.integrated, features = "IRF8", pt.size = 0.3, min.cutoff = 0) #pDC
FeaturePlot(all.integrated, features = "LILRA4", pt.size = 0.3, min.cutoff = 0) #pDC


#13
FeaturePlot(all.integrated, features = "EDN1", pt.size = 0.3, min.cutoff = 0) #Endothelial cell
FeaturePlot(all.integrated, features = "EPAS1", pt.size = 0.3, min.cutoff = 0) #Endothelial cell

#14
FeaturePlot(all.integrated, features = "CD68", pt.size = 0.3, min.cutoff = 0) #Dividing Mac
FeaturePlot(all.integrated, features = "MARCO", pt.size = 0.3, min.cutoff = 0) #Dividing Mac

#Remove unindetifiable clusters
new.cluster.id <- c("T1","T2","AM1","AM2","Mac1","B","Mac2","Epithelial","Plasma","T3","Mac3","NK","Mast","Endothelial","Dividing Mac","None1","None2","DC","None3","None4")
names(new.cluster.id) <- levels(all.integrated)
all.integrated <- RenameIdents(all.integrated, new.cluster.id)
human.lung.sorted <- subset(all.integrated, idents = c("T1","T2","AM1","AM2","Mac1","B","Mac2","Epithelial","Plasma","T3","Mac3","NK","Mast","Endothelial","Dividing Mac","DC"))

#Visualization
DimPlot(human.lung.sorted, reduction = "umap", label = T, label.size = 6, pt.size = 0.1) + NoLegend()
DimPlot(human.lung.sorted, reduction = "umap", label = T, label.size = 0, pt.size = 0.1) + NoLegend()

DimPlot(human.lung.sorted, reduction = "umap", group.by = "Tissue" , pt.size = 0.1)
DimPlot(human.lung.sorted, reduction = "umap", group.by = "Tissue" , split.by = "Tissue", pt.size = 0.1)
DimPlot(human.lung.sorted, reduction = "umap", group.by = "Tissue" , split.by = "PatientID" , pt.size = 0.1)
DimPlot(human.lung.sorted, reduction = "umap", split.by = "Tissue", pt.size = 0.1)

#Gene expressions
FeaturePlot(human.lung.sorted, features = "CD68", pt.size = 0.3, min.cutoff = 0) #Alveolar MP marker
FeaturePlot(human.lung.sorted, features = "MRC1", pt.size = 0.3, min.cutoff = 0) #Alveolar MP marker
FeaturePlot(human.lung.sorted, features = "CD163", pt.size = 0.3, min.cutoff = 0) #Alveolar MP marker
FeaturePlot(human.lung.sorted, features = "SIGLEC1", pt.size = 0.3, min.cutoff = 0) #Alveolar MP marker
FeaturePlot(human.lung.sorted, features = "MARCO", pt.size = 0.3, min.cutoff = 0) #Alveolar MP marker
FeaturePlot(human.lung.sorted, features = "IL10", pt.size = 0.3, min.cutoff = 0) #Alveolar MP marker
FeaturePlot(human.lung.sorted, features = "LY6E", pt.size = 0.3, min.cutoff = 0, split.by = c("Tissue")) #Alveolar MP marker

FeaturePlot(human.lung.sorted, features = "INHBA", pt.size = 0.3, min.cutoff = 0)
FeaturePlot(human.lung.sorted, features = "INHBA", pt.size = 0.3, min.cutoff = 0, split.by = c("Tissue"))

#Extract Alveolar Mac cluster2(AM2)
AM2<- subset(human.lung.sorted, idents = c("AM2"))
VlnPlot(AM2, feature = "INHBA", split.by = "Tissue", pt.size = 0) + geom_boxplot(outlier.shape = NA, width=0.1, alpha = 0.2, color = "black", position = position_dodge(0.9))
VlnPlot(AM2, feature = "INHBA", split.by = "Tissue", pt.size = 0) + stat_summary(fun = mean, geom='point', size = 25, colour = "black", shape = 95, position = position_dodge(0.9))

#Differentially expression genes, Normal vs Tumor in AM2
AM2.a <- AM2
AM2.a$Cluster.condition <- paste(AM2.a$seurat_clusters, AM2.a$Tissue, sep = "_")
Idents(AM2.a) <- AM2.a$Cluster.condition
DimPlot(AM2.a, reduction = "umap", label = T, label.size = 6, pt.size = 0.1) + NoLegend()
AM2.a.DEG <- FindMarkers(AM2.a, ident.1 = "3_LUAD", ident.2 = "3_Normal", logfc.threshold = 0, min.pct = 0)
write.csv(AM2.a.DEG, "./AM2_Tumor vs Normal.csv")

