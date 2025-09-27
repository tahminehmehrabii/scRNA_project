if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")

library(Seurat)
library(dplyr)
library(patchwork)

# Load example CSV
counts <- read.csv("my_scRNA_counts.csv", row.names = 1)
seurat_obj <- CreateSeuratObject(counts = counts, project = "scRNA_example")

# QC (simulated dataset, minimal QC)
seurat_obj[["percent.mt"]] <- 0  # no mitochondrial genes in this toy data

# Normalize data
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seurat_obj), 10)
VariableFeaturePlot(seurat_obj) + LabelPoints(points = top10, repel = TRUE)

# Scaling
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

# PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
DimPlot(seurat_obj, reduction = "pca")

# Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:5)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:5)
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()

# Marker genes (toy example)
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

# Save Seurat object
saveRDS(seurat_obj, file = "scRNA_example_seurat.rds")
