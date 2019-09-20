library(scater)
library(Seurat)
load("data/PRJEB11202_SCE_norm.RData")

y_genes <- rowData(sce)$symbol[rowData(sce)$chr == "Y"]
mito_genes <- rowData(sce)$symbol[rowData(sce)$chr == "MT"]

# Initial Seurat analysis -------------------------------------------------

sce <- CreateSeuratObject(counts = counts(sce),
                          meta.data = as.data.frame(colData(sce)))

sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(sce, features = rownames(sce))

# PCA
sce <- RunPCA(sce, features = VariableFeatures(object = sce), npcs = 30)
DimPlot(sce, reduction = "pca")
ElbowPlot(sce)

# KNN clustering using 10 dims of PCA space
sce <- FindNeighbors(sce, dims = 1:5)
sce <- FindClusters(sce, resolution = 0.2)

# UMAP plots of clusters and GATA3 expression
sce <- RunUMAP(sce, dims = 1:5, reduction = "pca")
DimPlot(sce, reduction = "umap")
FeaturePlot(sce, features = c("GATA3", "VCX"), reduction = "umap")
FeaturePlot(sce, features = c("DDX3Y", "RPS4Y1", "TTTY15", "KDM5D"), reduction = "umap")


# Removal of genes from chromosome Y --------------------------------------

load("data/PRJEB11202_SCE_norm.RData")

# Remove genes in the Y chromosome
sce <- sce[!rowData(sce)$chr %in% c("Y"), ]

sce <- CreateSeuratObject(counts = counts(sce),
                          meta.data = as.data.frame(colData(sce)))
#sce <- PercentageFeatureSet(sce, features = y_genes, col.name = "percent.y")

# Find var features and scale data
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(sce, features = rownames(sce))
#sce <- ScaleData(sce, features = rownames(sce), vars.to.regress = "percent.y")

# PCA
sce <- RunPCA(sce, features = VariableFeatures(object = sce), npcs = 30)
DimPlot(sce, reduction = "pca")
ElbowPlot(sce)

# KNN clustering using 10 dims of PCA space
sce <- FindNeighbors(sce, dims = 1:5)
sce <- FindClusters(sce, resolution = 0.2)

# UMAP plots of clusters and GATA3 expression
sce <- RunUMAP(sce, dims = 1:5, reduction = "pca")
DimPlot(sce, reduction = "umap")
FeaturePlot(sce, features = c("GATA3"), reduction = "umap")

# Differential gene expression analysis
dgea <- FindMarkers(sce, ident.1 = 1, ident.2 = NULL)

ggplot(dgea, aes(avg_logFC)) + geom_histogram(bins = 20) + 
  labs(x = "Avg. log fold-change", y = "Counts") + 
  theme_classic(base_size = 15) 

# Select genes that are actve in each population
gata3_pos <- rownames(dgea)[dgea$avg_logFC > 0.25 & dgea$p_val_adj < 0.001]
gata3_neg <- rownames(dgea)[dgea$avg_logFC < -0.25 & dgea$p_val_adj < 0.001]

DoHeatmap(sce, features = c(gata3_pos, gata3_neg))

save(sce, file = "results/Seurat_Morula.RData")
save(dgea, gata3_pos, gata3_neg, file = "results/DGEA.RData")
