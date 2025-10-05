# ========== 1. Load libraries ==========
library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(SingleR)
library(celldex)
library(infercnv)
library(DoubletFinder)
library(Matrix)
library(scater)
library(scran)
library(readr)

# ========== 2. Load raw count data ==========
# Replace the path below with your dataset folder (e.g., 10X matrix)
counts <- Read10X(data.dir = "D:/scRNA/raw_data/")  
seurat <- CreateSeuratObject(counts = counts, project = "HCC_scRNA", min.cells = 3, min.features = 200)

# ========== 3. Basic QC and filtering ==========
# Compute % mitochondrial genes
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC thresholds
seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & 
                 nCount_RNA > 500 & percent.mt < 20)

# ========== 4. Remove potential doublets ==========
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, nfeatures = 2000)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, npcs = 30)

sweep.res <- paramSweep_v3(seurat, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), "pK"]))
nExp_poi <- round(0.06 * ncol(seurat))

seurat <- doubletFinder_v3(seurat, PCs = 1:10, pN = 0.25, pK = pK, 
                           nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
df_col <- grep("DF.classifications", colnames(seurat@meta.data), value = TRUE)[1]
seurat <- subset(seurat, subset = (!!as.name(df_col)) == "Singlet")

# ========== 5. Normalization & scaling (SCTransform) ==========
seurat <- SCTransform(seurat, vars.to.regress = "percent.mt", verbose = FALSE)

# ========== 6. PCA, batch correction, and clustering ==========
seurat <- RunPCA(seurat, npcs = 30, verbose = FALSE)
# If batch column exists, run Harmony
if ("batch" %in% colnames(seurat@meta.data)) {
  seurat <- RunHarmony(seurat, group.by.vars = "batch", reduction = "pca", dims.use = 1:20)
  seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:20)
  seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:20)
} else {
  seurat <- RunUMAP(seurat, dims = 1:20)
  seurat <- FindNeighbors(seurat, dims = 1:20)
}
seurat <- FindClusters(seurat, resolution = 0.4)
DimPlot(seurat, reduction = "umap", label = TRUE) + NoLegend()

# ========== 7. Marker gene detection ==========
markers_all <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers_all, "Results/sc_markers_all_clusters.csv", row.names = FALSE)

# ========== 8. Cell type annotation (SingleR) ==========
ref <- celldex::HumanPrimaryCellAtlasData()
sce <- as.SingleCellExperiment(seurat)
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
seurat$SingleR_label <- pred$labels
DimPlot(seurat, group.by = "SingleR_label", label = TRUE, repel = TRUE)

# ========== 9. CNV detection (inferCNV) ==========
# Prepare annotation file (cell_name <tab> cluster)
cell_annot <- data.frame(cell = colnames(seurat), cluster = seurat$seurat_clusters)
write.table(cell_annot, "infercnv_cell_annotations.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Raw counts (non-normalized)
raw_counts <- GetAssayData(seurat, slot = "counts", assay = "RNA")

# Create inferCNV object
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = as.matrix(raw_counts),
  annotations_file = "infercnv_cell_annotations.txt",
  delim = "\t",
  gene_order_file = "gene_position_file.txt",  # must be prepared externally
  ref_group_names = c("normal")                # replace with your normal cell cluster name
)

# Run inferCNV
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 1,
                              out_dir = "Results/infercnv_out",
                              cluster_by_groups = TRUE,
                              denoise = TRUE,
                              HMM = TRUE)

# ========== 10. Save outputs ==========
saveRDS(seurat, file = "Results/seurat_processed.rds")
write.csv(seurat@meta.data, "Results/seurat_metadata.csv", row.names = TRUE)
save.image("Results/scRNA_pipeline_workspace.RData")



