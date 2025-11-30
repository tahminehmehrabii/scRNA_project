############################################################
# LSCC scRNA-seq pipeline (GSE150321)
############################################################

suppressPackageStartupMessages({
  ## scRNA
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
  library(tools)
  library(SingleCellExperiment)
  library(scDblFinder)})

## Main seed for all analyses (scRNA + Bulk + ML)
set.seed(1234)

############################################################
# Paths
############################################################

data_dir <- "D:/Single Cell/GSE150321_RAW"
out_dir  <- "D:/Single Cell/GSE150321_RAW"
setwd(data_dir)

project_id <- "LSCC_GSE150321"

############################################################
# LSCC01 
############################################################

file_path1 <- file.path(data_dir, "GSM4546857_LSCC01_DBEC_UMI.csv.gz")
sample_id1 <- "LSCC01"
dt1 <- fread(file = file_path1)
print(dt1[1:3, 1:4])
any(duplicated(dt1$Cell_Index))

# Build expression matrix
expr_mat1 <- as.matrix(dt1[, -1, with = FALSE])
any(duplicated(dt1[[1]]))
rownames(expr_mat1) <- dt1[[1]]
print(dim(expr_mat1))

# Main fix: transpose → gene x cell
expr_mat1 <- t(expr_mat1)
print(dim(expr_mat1))
all(expr_mat1 == floor(expr_mat1))
expr_df1 <- as.data.frame(expr_mat1)
print(expr_df1[1:2, 1:4])

# Build Seurat object
seu1 <- CreateSeuratObject(
  counts       = expr_mat1,
  project      = sample_id1,
  min.cells    = 1,
  min.features = 1
)
seu1$sample_id <- sample_id1
meta_df1 <- seu1@meta.data
print(head(meta_df1, 5))

############################################################
# LSCC02 
############################################################

file_path2 <- file.path(data_dir, "GSM4546858_LSCC02_DBEC_UMI.csv.gz")
sample_id2 <- "LSCC02"
dt2 <- fread(file = file_path2)
any(duplicated(colnames(dt2)[-1]))
print(dt2[1:3, 1:4])

# Build expression matrix
expr_mat2 <- as.matrix(dt2[, -1, with = FALSE])
any(duplicated(dt2[[1]]))
rownames(expr_mat2) <- dt2[[1]]
print(dim(expr_mat2))
all(expr_mat2 == floor(expr_mat2))
expr_df2 <- as.data.frame(expr_mat2)
print(expr_df2[1:2, 1:4])

# Build Seurat object
seu2 <- CreateSeuratObject(
  counts       = expr_mat2,
  project      = sample_id2,
  min.cells    = 1,
  min.features = 1
)
seu2$sample_id <- sample_id2
meta_df2 <- seu2@meta.data
print(head(meta_df2, 5))

############################################################
# STEP 1 — cell counts before QC
############################################################

n_cells1 <- ncol(seu1)
n_cells2 <- ncol(seu2)

cat("\n[Step 1] Number of cells before QC:\n")
cat("Sample", sample_id1, ":", n_cells1, "\n")
cat("Sample", sample_id2, ":", n_cells2, "\n")

stage1_qc <- data.frame(
  sample_id = c(sample_id1, sample_id2),
  n_cells   = c(n_cells1,   n_cells2),
  stage     = c("raw_before_QC", "raw_before_QC")
)
cat("\n[Step 1] stage1_qc:\n")
print(stage1_qc)

## Save stats before QC
write.csv(
  stage1_qc,
  file = file.path(data_dir, paste0(project_id, "_Stage1_preQC_cell_counts.csv")),
  row.names = FALSE
)

############################################################
# STEP 2 — QC per sample
############################################################

minFeature <- 200
maxFeature <- 7500
minCount   <- 400
maxCount   <- 40000
maxMT      <- 20

## LSCC01
seu1[["percent.mt"]] <- PercentageFeatureSet(seu1, pattern = "^MT-")
cat("\n[Step 2] LSCC01 meta before QC – first five rows and six columns:\n")
print(head(seu1@meta.data, 5))

seu1_qc <- subset(
  seu1,
  subset = nFeature_RNA > minFeature &
    nFeature_RNA < maxFeature &
    nCount_RNA    > minCount &
    nCount_RNA    < maxCount &
    percent.mt    < maxMT
)

cat("\n[Step 2] LSCC01 after QC – dimensions and some metadata rows:\n")
cat("Dimensions:", nrow(seu1_qc), "genes x", ncol(seu1_qc), "cells\n")
print(head(seu1_qc@meta.data, 5))

## LSCC02
seu2[["percent.mt"]] <- PercentageFeatureSet(seu2, pattern = "^MT-")
cat("\n[Step 2] LSCC02 meta before QC – first five rows and six columns:\n")
print(head(seu2@meta.data, 5))

seu2_qc <- subset(
  seu2,
  subset = nFeature_RNA > minFeature &
    nFeature_RNA < maxFeature &
    nCount_RNA    > minCount &
    nCount_RNA    < maxCount &
    percent.mt    < maxMT
)

cat("\n[Step 2] LSCC02 after QC – dimensions and some metadata rows:\n")
cat("Dimensions:", nrow(seu2_qc), "genes x", ncol(seu2_qc), "cells\n")
print(head(seu2_qc@meta.data, 5))

## Number of cells after QC
n_cells_qc1 <- ncol(seu1_qc)
n_cells_qc2 <- ncol(seu2_qc)

cat("\n[Step 2] Number of cells after QC:\n")
cat("Sample", sample_id1, ":", n_cells_qc1, "\n")
cat("Sample", sample_id2, ":", n_cells_qc2, "\n")

stage2_qc <- data.frame(
  sample_id = c(sample_id1, sample_id2),
  n_cells   = c(n_cells_qc1, n_cells_qc2),
  stage     = c("post_QC", "post_QC")
)
cat("\n[Step 2] stage2_qc:\n")
print(stage2_qc)

write.csv(
  stage2_qc,
  file = file.path(data_dir, paste0(project_id, "_Stage2_postQC_cell_counts.csv")),
  row.names = FALSE
)

############################################################
# filter-wise removed cell percentages and QC violin plots
############################################################

# Function to compute number and percentage of cells removed by each filter
qc_filter_stats <- function(seu_obj,
                            sample_label,
                            minFeature, maxFeature,
                            minCount,   maxCount,
                            maxMT) {
  md <- seu_obj@meta.data
  if (!all(c("nFeature_RNA", "nCount_RNA", "percent.mt") %in% colnames(md))) {
    warning("QC columns not found in meta.data for sample: ", sample_label)
    return(NULL)
  }
  
  total_cells <- nrow(md)
  
  f_min_feat <- md$nFeature_RNA > minFeature
  f_max_feat <- md$nFeature_RNA < maxFeature
  f_min_cnt  <- md$nCount_RNA   > minCount
  f_max_cnt  <- md$nCount_RNA   < maxCount
  f_mt       <- md$percent.mt   < maxMT
  
  # All filters combined (the same condition we use in subset)
  f_all <- f_min_feat & f_max_feat & f_min_cnt & f_max_cnt & f_mt
  
  kept_vec <- c(
    sum(f_min_feat),
    sum(f_max_feat),
    sum(f_min_cnt),
    sum(f_max_cnt),
    sum(f_mt),
    sum(f_all)
  )
  
  removed_vec <- total_cells - kept_vec
  removed_pct <- round(100 * removed_vec / total_cells, 1)
  
  data.frame(
    sample_id     = sample_label,
    filter        = c(
      paste0("nFeature_RNA > ", minFeature),
      paste0("nFeature_RNA < ", maxFeature),
      paste0("nCount_RNA > ", minCount),
      paste0("nCount_RNA < ", maxCount),
      paste0("percent.mt < ", maxMT),
      "All filters combined"
    ),
    total_cells   = total_cells,
    kept_cells    = kept_vec,
    removed_cells = removed_vec,
    removed_pct   = removed_pct,
    stringsAsFactors = FALSE
  )
}

## Compute per-filter tables for each sample
qc_stats1 <- qc_filter_stats(
  seu_obj      = seu1,
  sample_label = sample_id1,
  minFeature   = minFeature,
  maxFeature   = maxFeature,
  minCount     = minCount,
  maxCount     = maxCount,
  maxMT        = maxMT
)

qc_stats2 <- qc_filter_stats(
  seu_obj      = seu2,
  sample_label = sample_id2,
  minFeature   = minFeature,
  maxFeature   = maxFeature,
  minCount     = minCount,
  maxCount     = maxCount,
  maxMT        = maxMT
)

qc_stats_all <- rbind(qc_stats1, qc_stats2)

cat("\n[Step 2] QC filter statistics (per filter and combined):\n")
print(qc_stats_all)

qc_stats_file <- file.path(
  data_dir,
  paste0(project_id, "_Step2_QC_filter_stats_per_sample.csv")
)
write.csv(qc_stats_all, qc_stats_file, row.names = FALSE)
cat("\n[Step 2] QC filter statistics saved to:\n", qc_stats_file, "\n")

## QC violin plots for each sample (on all cells before QC)

# LSCC01
p_nf_1 <- VlnPlot(seu1, features = "nFeature_RNA", pt.size = 0.1) +
  geom_hline(yintercept = minFeature, linetype = "dashed", color = "red") +
  geom_hline(yintercept = maxFeature, linetype = "dashed", color = "red") +
  ggtitle(paste0(sample_id1, " – nFeature_RNA with QC thresholds"))

p_nc_1 <- VlnPlot(seu1, features = "nCount_RNA", pt.size = 0.1) +
  geom_hline(yintercept = minCount, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = maxCount, linetype = "dashed", color = "blue") +
  ggtitle(paste0(sample_id1, " – nCount_RNA with QC thresholds"))

p_mt_1 <- VlnPlot(seu1, features = "percent.mt", pt.size = 0.1) +
  geom_hline(yintercept = maxMT, linetype = "dashed", color = "purple") +
  ggtitle(paste0(sample_id1, " – percent.mt with QC threshold"))

ggsave(
  filename = file.path(data_dir, paste0(project_id, "_Step2_QC_", sample_id1, "_nFeature.png")),
  plot     = p_nf_1,
  width    = 5, height = 4, dpi = 300
)
ggsave(
  filename = file.path(data_dir, paste0(project_id, "_Step2_QC_", sample_id1, "_nCount.png")),
  plot     = p_nc_1,
  width    = 5, height = 4, dpi = 300
)
ggsave(
  filename = file.path(data_dir, paste0(project_id, "_Step2_QC_", sample_id1, "_percentMT.png")),
  plot     = p_mt_1,
  width    = 5, height = 4, dpi = 300
)

# LSCC02
p_nf_2 <- VlnPlot(seu2, features = "nFeature_RNA", pt.size = 0.1) +
  geom_hline(yintercept = minFeature, linetype = "dashed", color = "red") +
  geom_hline(yintercept = maxFeature, linetype = "dashed", color = "red") +
  ggtitle(paste0(sample_id2, " – nFeature_RNA with QC thresholds"))

p_nc_2 <- VlnPlot(seu2, features = "nCount_RNA", pt.size = 0.1) +
  geom_hline(yintercept = minCount, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = maxCount, linetype = "dashed", color = "blue") +
  ggtitle(paste0(sample_id2, " – nCount_RNA with QC thresholds"))

p_mt_2 <- VlnPlot(seu2, features = "percent.mt", pt.size = 0.1) +
  geom_hline(yintercept = maxMT, linetype = "dashed", color = "purple") +
  ggtitle(paste0(sample_id2, " – percent.mt with QC threshold"))

ggsave(
  filename = file.path(data_dir, paste0(project_id, "_Step2_QC_", sample_id2, "_nFeature.png")),
  plot     = p_nf_2,
  width    = 5, height = 4, dpi = 300
)
ggsave(
  filename = file.path(data_dir, paste0(project_id, "_Step2_QC_", sample_id2, "_nCount.png")),
  plot     = p_nc_2,
  width    = 5, height = 4, dpi = 300
)
ggsave(
  filename = file.path(data_dir, paste0(project_id, "_Step2_QC_", sample_id2, "_percentMT.png")),
  plot     = p_mt_2,
  width    = 5, height = 4, dpi = 300
)

############################################################
# EXTRA2 — summary of cell counts before and after QC
############################################################

# Number of cells before QC (from original Seurat objects)
n_cells_raw1 <- ncol(seu1)
n_cells_raw2 <- ncol(seu2)

# Summary before/after for each sample
qc_summary <- data.frame(
  sample_id       = c(sample_id1, sample_id2),
  n_cells_before  = c(n_cells_raw1, n_cells_raw2),
  n_cells_after   = c(n_cells_qc1,  n_cells_qc2),
  n_removed       = c(n_cells_raw1 - n_cells_qc1,
                      n_cells_raw2 - n_cells_qc2),
  removed_pct     = round(
    100 * c(n_cells_raw1 - n_cells_qc1,
            n_cells_raw2 - n_cells_qc2) /
      c(n_cells_raw1, n_cells_raw2),
    1
  ),
  stringsAsFactors = FALSE
)

cat("\n[Step 2] QC summary (before/after and % removed):\n")
print(qc_summary)

# If the difference in % removed between samples is large, print a note
diff_removed <- abs(qc_summary$removed_pct[1] - qc_summary$removed_pct[2])
if (diff_removed > 10) {
  cat("\n[Step 2] NOTE: Difference in % removed between samples is ",
      diff_removed, "% — please mention this potential bias in the manuscript.\n", sep = "")
}

write.csv(
  qc_summary,
  file = file.path(data_dir, paste0(project_id, "_Step2_QC_summary_before_after_removed.csv")),
  row.names = FALSE
)

############################################################
# STEP 3 — Merge
############################################################

lscc <- merge(
  x            = seu1_qc,
  y            = seu2_qc,
  add.cell.ids = c(sample_id1, sample_id2),
  project      = project_id
)

cat("\n[Step 3] Merged Seurat object:\n")
print(lscc)

cat("\n[Step 3] Cells per sample_id:\n")
print(table(lscc$sample_id))

cat("\n[Step 3] merged meta – first five rows and seven columns:\n")
print(head(lscc@meta.data, 5))

meta_stage3 <- lscc@meta.data
write.csv(
  meta_stage3,
  file = file.path(data_dir, paste0(project_id, "_Stage3_meta_merged.csv")),
  row.names = TRUE
)

lscc <- JoinLayers(lscc)
DefaultAssay(lscc) <- "RNA"

saveRDS(
  lscc,
  file = file.path(data_dir, paste0(project_id, "_Stage3_merged_seurat.rds"))
)

############################################################
# STEP 3b — Doublet detection
############################################################

cat("\n[Step 3b] Doublet detection with scDblFinder ...\n")

sce <- as.SingleCellExperiment(lscc)
colData(sce)$sample_id <- lscc$sample_id

set.seed(1234)  ## For reproducibility of scDblFinder
sce <- scDblFinder(sce, samples = "sample_id")

lscc$doublet_class <- sce$scDblFinder.class
lscc$doublet_score <- sce$scDblFinder.score

cat("\n[Step 3b] Doublet classes (before filtering):\n")
print(table(lscc$doublet_class))

cat("\n[Step 3b] meta with doublet_class – first five rows and eight columns:\n")
print(head(lscc@meta.data, 5))

lscc <- subset(lscc, subset = doublet_class == "singlet")

cat("\n[Step 3b] Number of cells after doublet removal:", ncol(lscc), "\n")
cat("\n[Step 3b] meta after doublet removal – first five rows and eight columns:\n")
print(head(lscc@meta.data, 5))

write.csv(
  lscc@meta.data,
  file = file.path(data_dir, paste0(project_id, "_Stage3b_meta_after_doubletRemoval.csv")),
  row.names = TRUE
)

############################################################
# STEP 4 — Normalize, HVG, Scale, PCA
############################################################

cat("\n[Step 4] NormalizeData ...\n")
lscc <- NormalizeData(lscc)

cat("\n[Step 4] FindVariableFeatures ...\n")
lscc <- FindVariableFeatures(
  lscc,
  selection.method = "vst",
  nfeatures        = 3000
)

cat("\n[Step 4] Number of variable features:",
    length(VariableFeatures(lscc)), "\n")

cat("\n[Step 4] First 10 HVGs:\n")
print(head(VariableFeatures(lscc), 10))

hvg_df <- data.frame(gene = VariableFeatures(lscc))
write.csv(
  hvg_df,
  file = file.path(data_dir, paste0(project_id, "_Stage4_HVG_genes.csv")),
  row.names = FALSE
)

cat("\n[Step 4] ScaleData on HVGs ...\n")
lscc <- ScaleData(
  lscc,
  features = VariableFeatures(lscc)
)

cat("\n[Step 4] RunPCA ...\n")
lscc <- RunPCA(
  lscc,
  features = VariableFeatures(lscc)
)

cat("\n[Step 4] First 5 cells and first 3 PCs from PCA embeddings:\n")
print(head(Embeddings(lscc, "pca")[, 1:3], 5))

p_elbow <- ElbowPlot(lscc, ndims = 50) +
  ggtitle("Stage4 - PCA elbow plot")

ggsave(
  file.path(data_dir, paste0(project_id, "_Stage4_PCA_Elbow.png")),
  plot   = p_elbow,
  width  = 6,
  height = 5,
  dpi    = 1200
)

saveRDS(
  lscc,
  file = file.path(data_dir, paste0(project_id, "_Stage4_postPCA_seurat.rds"))
)

############################################################
# STEP 5 — Neighbors, clustering, UMAP
############################################################

dims_use <- 1:30

cat("\n[Step 5] FindNeighbors ...\n")
lscc <- FindNeighbors(
  lscc,
  dims = dims_use
)

cat("\n[Step 5] FindClusters ...\n")
lscc <- FindClusters(
  lscc,
  resolution = 0.5
)

cat("\n[Step 5] Cell counts per cluster:\n")
print(table(lscc$seurat_clusters))

set.seed(1234)  ## For reproducible UMAP
cat("\n[Step 5] RunUMAP ...\n")
lscc <- RunUMAP(
  lscc,
  dims = dims_use
)

cat("\n[Step 5] First 5 cells in UMAP coordinates:\n")
print(head(Embeddings(lscc, "umap"), 5))

lscc$celltype_main <- as.character(lscc$seurat_clusters)

p_umap_cluster <- DimPlot(
  lscc,
  reduction = "umap",
  group.by  = "seurat_clusters",
  label     = TRUE,
  repel     = TRUE
) +
  ggtitle("Stage5 - Global UMAP - clusters") +
  NoLegend()

p_umap_sample <- DimPlot(
  lscc,
  reduction = "umap",
  group.by  = "sample_id"
) +
  ggtitle("Stage5 - Global UMAP - samples")

ggsave(
  file.path(data_dir, paste0(project_id, "_Stage5_UMAP_clusters.png")),
  plot   = p_umap_cluster,
  width  = 7,
  height = 6,
  dpi    = 1200
)

ggsave(
  file.path(data_dir, paste0(project_id, "_Stage5_UMAP_samples.png")),
  plot   = p_umap_sample,
  width  = 7,
  height = 6,
  dpi    = 1200
)

write.csv(
  lscc@meta.data,
  file = file.path(data_dir, paste0(project_id, "_Stage5_meta_with_clusters.csv")),
  row.names = TRUE
)

saveRDS(
  lscc,
  file = file.path(data_dir, paste0(project_id, "_Stage5_postUMAP_seurat.rds"))
)

############################################################
# EXTRA (Step 5) — Batch effect check between LSCC01 and LSCC02
############################################################

cat("\n[Step 5 – EXTRA] Batch check based on sample_id vs clusters:\n")

# Table: cluster × sample
tab_cs <- table(lscc$seurat_clusters, lscc$sample_id)
cat("\n[Step 5 – EXTRA] Cell counts per (cluster, sample):\n")
print(tab_cs)

# Proportion of each sample inside each cluster (rows = cluster, sum=1 per row)
prop_by_cluster <- prop.table(tab_cs, margin = 1)
cat("\n[Step 5 – EXTRA] Proportion of each sample inside each cluster (rows sum to 1):\n")
print(round(prop_by_cluster, 3))

# Proportion of each cluster inside each sample (columns = sample, sum=1 per column)
prop_by_sample <- prop.table(tab_cs, margin = 2)
cat("\n[Step 5 – EXTRA] Proportion of each cluster inside each sample (columns sum to 1):\n")
print(round(prop_by_sample, 3))

# Simple chi-square test for cluster vs sample association
if (all(dim(tab_cs) > 1)) {
  chi_res <- chisq.test(tab_cs)
  cat("\n[Step 5 – EXTRA] Chi-square test for association between clusters and sample_id:\n")
  print(chi_res)
} else {
  cat("\n[Step 5 – EXTRA] Not enough levels in cluster/sample to run chi-square test.\n")
}

# Save tables to files for paper / Supplement
write.csv(
  as.data.frame(tab_cs),
  file = file.path(data_dir,
                   paste0(project_id, "_Step5_table_cells_per_cluster_sample.csv")),
  row.names = FALSE
)

write.csv(
  as.data.frame(prop_by_cluster),
  file = file.path(data_dir,
                   paste0(project_id, "_Step5_prop_samples_within_each_cluster.csv")),
  row.names = FALSE
)

write.csv(
  as.data.frame(prop_by_sample),
  file = file.path(data_dir,
                   paste0(project_id, "_Step5_prop_clusters_within_each_sample.csv")),
  row.names = FALSE
)

cat("\n[Step 5 – EXTRA] Batch-effect tables saved.\n")

############################################################
# STEP 6 (scRNA) — Global cluster markers (up + down)
############################################################

cat("\n[Step 6] FindAllMarkers for all clusters ...\n")

all_markers <- FindAllMarkers(
  lscc,
  only.pos        = FALSE,
  min.pct         = 0.25,
  logfc.threshold = 0.25
)

cat("\n[Step 6] Dimensions of marker table (rows = gene-cluster pairs):\n")
print(dim(all_markers))

cat("\n[Step 6] First 10 marker rows:\n")
print(head(all_markers, 10))

# Direction of expression (up / down in each cluster)
all_markers$direction <- ifelse(all_markers$avg_log2FC > 0, "up", "down")

write.csv(
  all_markers,
  file = file.path(data_dir,
                   paste0(project_id, "_Stage6_cluster_markers_global_up_down.csv")),
  row.names = FALSE
)

cat("\n===== Script finished up to Step 6. Continuing with further steps ... =====\n")

############################################################
# STEP 7 (scRNA) — Canonical lineage marker DotPlot
############################################################

cat("\n[Step 7] Building marker panel and drawing DotPlot for major lineages ...\n")

t_markers <- c("CD3D", "CD3E", "CD2", "TRAC", "CD4", "CD8A", "CD8B")
b_markers <- c("MS4A1", "CD79A", "CD79B", "CD19", "CD74", "MZB1", "SDC1")
myeloid_markers <- c("LYZ", "CD68", "CD163", "LST1", "CSF1R", "ITGAM", "FCGR3A")
epith_tumor_markers <- c("EPCAM", "KRT8", "KRT18", "KRT5", "KRT14", "KRT17", "TP63")
epith_hnscc_extra <- c("KRT4", "KRT6A", "KRT6B", "KRT16", "KRT1", "KRT10", "KRT13")
fibro_markers <- c("COL1A1", "COL1A2", "COL3A1", "TAGLN", "ACTA2", "PDGFRA", "DCN")
endo_markers <- c("PECAM1", "CDH5", "VWF", "KDR", "ENG", "RGS5")
prolif_markers <- c("MKI67", "TOP2A", "PCNA", "BIRC5")
mast_markers <- c("TPSAB1", "TPSB2", "CPA3")
nk_markers <- c("NKG7", "GNLY", "PRF1", "GZMB", "CTSW")
tam_markers <- c("MRC1", "CD163", "SPP1", "APOE")

marker_panel <- unique(c(
  t_markers,
  b_markers,
  myeloid_markers,
  epith_tumor_markers,
  epith_hnscc_extra,
  fibro_markers,
  endo_markers,
  prolif_markers,
  mast_markers,
  nk_markers,
  tam_markers
))

## Keep only genes that actually exist in lscc
marker_panel <- marker_panel[marker_panel %in% rownames(lscc)]

cat("\n[Step 7] Number of genes in marker_panel present in lscc:", length(marker_panel), "\n")

write.csv(
  data.frame(marker = marker_panel),
  file = file.path(data_dir, paste0(project_id, "_Stage7_marker_panel_used.csv")),
  row.names = FALSE
)

if (length(marker_panel) > 0) {
  p_dot <- DotPlot(lscc, features = marker_panel) +
    RotatedAxis() +
    ggtitle("Stage7 - Canonical & LSCC lineage markers across clusters")
  
  ggsave(
    file.path(data_dir, paste0(project_id, "_Stage7_lineage_DotPlot.png")),
    plot   = p_dot,
    width  = 10,
    height = 6,
    dpi    = 1200
  )
  
  cat("\n[Step 7] DotPlot saved: Stage7_lineage_DotPlot.png\n")
} else {
  cat("\n[Step 7] None of the marker_panel genes were found in lscc.\n")
}

############################################################
# STEP 8 (scRNA) — Manual mapping: cluster -> major cell type
############################################################

cat("\n[Step 8] Manual mapping: cluster -> major cell type ...\n")

cat("\n[Step 8] Current cluster levels:\n")
print(levels(lscc))  # Idents(lscc) = seurat_clusters here

## Mapping based on mean expression of markers (lineage_scores)
cluster_to_celltype <- c(
  "T cells",           # cluster 0
  "Tumor/Epithelial",  # cluster 1
  "T cells",           # cluster 2
  "Tumor/Epithelial",  # cluster 3
  "T cells",           # cluster 4
  "B/Plasma cells",    # cluster 5
  "Tumor/Epithelial",  # cluster 6
  "Tumor/Epithelial",  # cluster 7
  "B/Plasma cells",    # cluster 8
  "Fibroblasts",       # cluster 9
  "Myeloid",           # cluster 10
  "Endothelial",       # cluster 11
  "B/Plasma cells",    # cluster 12
  "Fibroblasts",       # cluster 13
  "B/Plasma cells",    # cluster 14
  "Tumor/Epithelial",  # cluster 15
  "B/Plasma cells",    # cluster 16
  "T cells"            # cluster 17
)

if (length(cluster_to_celltype) != length(levels(lscc))) {
  warning("[Step 8] length of cluster_to_celltype does not match number of clusters. celltype_main remains as cluster IDs.")
} else {
  names(cluster_to_celltype) <- levels(lscc)
  lscc <- RenameIdents(lscc, cluster_to_celltype)
  lscc$celltype_main <- Idents(lscc)
}

p_umap_ct <- DimPlot(
  lscc,
  reduction = "umap",
  group.by  = "celltype_main",
  label     = TRUE,
  repel     = TRUE
) +
  ggtitle("Stage8 - Global UMAP - main cell types")

ggsave(
  file.path(data_dir, paste0(project_id, "_Stage8_UMAP_celltypes.png")),
  plot   = p_umap_ct,
  width  = 7,
  height = 6,
  dpi    = 1200
)

write.csv(
  lscc@meta.data,
  file = file.path(data_dir, paste0(project_id, "_Stage8_meta_with_celltypes.csv")),
  row.names = TRUE
)

############################################################
# STEP 9 (scRNA) — Per-lineage subsetting & reclustering
############################################################

cat("\n[Step 9] Per-lineage subsetting & reclustering ...\n")

run_lineage_pipeline <- function(obj,
                                 subset_ids,
                                 lineage_name,
                                 max_dims  = 20,
                                 res       = 0.5,
                                 project_id,
                                 data_dir) {
  message("\nProcessing lineage: ", lineage_name)
  
  ## Keep only identities that actually exist in Idents(obj)
  valid_ids <- intersect(subset_ids, levels(Idents(obj)))
  if (length(valid_ids) == 0) {
    message("  No matching identities for ", lineage_name, " (skipping).")
    return(NULL)
  }
  
  sub <- subset(obj, idents = valid_ids)
  if (ncol(sub) < 50) {
    message("  Too few cells for ", lineage_name, " (", ncol(sub), "), skipping.")
    return(NULL)
  }
  
  sub <- NormalizeData(sub)
  sub <- FindVariableFeatures(sub, selection.method = "vst", nfeatures = 3000)
  sub <- ScaleData(sub, features = VariableFeatures(sub))
  
  npcs_use <- min(max_dims, 50, ncol(sub), nrow(sub))
  sub <- RunPCA(sub, features = VariableFeatures(sub), npcs = npcs_use)
  
  n_pcs_computed <- ncol(Embeddings(sub, "pca"))
  dims_use <- 1:min(max_dims, n_pcs_computed)
  message("  Using ", length(dims_use), " PCs (1:", max(dims_use), ") for ", lineage_name)
  
  sub <- FindNeighbors(sub, dims = dims_use)
  sub <- FindClusters(sub, resolution = res)
  
  set.seed(1234) ## For UMAP reproducibility in each lineage
  sub <- RunUMAP(sub, dims = dims_use)
  
  p1 <- DimPlot(
    sub,
    reduction = "umap",
    group.by  = "seurat_clusters",
    label     = TRUE,
    repel     = TRUE
  ) +
    ggtitle(paste0("Stage9 - ", lineage_name, " UMAP - clusters")) +
    NoLegend()
  
  p2 <- DimPlot(
    sub,
    reduction = "umap",
    group.by  = "sample_id"
  ) +
    ggtitle(paste0("Stage9 - ", lineage_name, " UMAP - samples"))
  
  ggsave(
    file.path(data_dir, paste0(project_id, "_Stage9_", lineage_name, "_UMAP_clusters.png")),
    plot   = p1,
    width  = 6,
    height = 5,
    dpi    = 1200
  )
  
  ggsave(
    file.path(data_dir, paste0(project_id, "_Stage9_", lineage_name, "_UMAP_samples.png")),
    plot   = p2,
    width  = 6,
    height = 5,
    dpi    = 1200
  )
  
  lineage_markers <- FindAllMarkers(
    sub,
    only.pos        = TRUE,
    min.pct         = 0.25,
    logfc.threshold = 0.25
  )
  
  write.csv(
    lineage_markers,
    file = file.path(data_dir, paste0(project_id, "_Stage9_", lineage_name, "_markers.csv")),
    row.names = FALSE
  )
  
  write.csv(
    sub@meta.data,
    file = file.path(data_dir, paste0(project_id, "_Stage9_", lineage_name, "_meta.csv")),
    row.names = TRUE
  )
  
  saveRDS(
    sub,
    file = file.path(data_dir, paste0(project_id, "_Stage9_", lineage_name, "_seurat.rds"))
  )
  
  return(sub)
}

## Define lineages based on celltype_main set in Step 8
t_lineage_ids       <- c("T cells")
b_lineage_ids       <- c("B/Plasma cells")
myeloid_lineage_ids <- c("Myeloid")
epith_lineage_ids   <- c("Tumor/Epithelial")
fibro_lineage_ids   <- c("Fibroblasts")
endo_lineage_ids    <- c("Endothelial")

t_obj <- run_lineage_pipeline(
  obj          = lscc,
  subset_ids   = t_lineage_ids,
  lineage_name = "T_cells",
  max_dims     = 20,
  res          = 0.5,
  project_id   = project_id,
  data_dir     = data_dir
)

b_obj <- run_lineage_pipeline(
  obj          = lscc,
  subset_ids   = b_lineage_ids,
  lineage_name = "B_cells",
  max_dims     = 20,
  res          = 0.5,
  project_id   = project_id,
  data_dir     = data_dir
)

myeloid_obj <- run_lineage_pipeline(
  obj          = lscc,
  subset_ids   = myeloid_lineage_ids,
  lineage_name = "Myeloid",
  max_dims     = 20,
  res          = 0.5,
  project_id   = project_id,
  data_dir     = data_dir
)

epith_obj <- run_lineage_pipeline(
  obj          = lscc,
  subset_ids   = epith_lineage_ids,
  lineage_name = "EpithelialTumor",
  max_dims     = 20,
  res          = 0.5,
  project_id   = project_id,
  data_dir     = data_dir
)

fibro_obj <- run_lineage_pipeline(
  obj          = lscc,
  subset_ids   = fibro_lineage_ids,
  lineage_name = "Fibroblasts",
  max_dims     = 20,
  res          = 0.5,
  project_id   = project_id,
  data_dir     = data_dir
)

endo_obj <- run_lineage_pipeline(
  obj          = lscc,
  subset_ids   = endo_lineage_ids,
  lineage_name = "Endothelial",
  max_dims     = 20,
  res          = 0.5,
  project_id   = project_id,
  data_dir     = data_dir
)

############################################################
# STEP 10 (scRNA) — Save final global Seurat object
############################################################

cat("\n[Step 10] Saving global annotated Seurat object ...\n")

saveRDS(
  lscc,
  file = file.path(data_dir, paste0(project_id, "_Stage10_Global_annotated_seurat.rds"))
)

cat("\n[Step 10] Global annotated Seurat object saved at:\n",
    file.path(data_dir, paste0(project_id, "_Stage10_Global_annotated_seurat.rds")), "\n")

############################################################
# STEP 11 (scRNA) — Epithelial labels (Malignant vs Keratinocyte_like)
############################################################

cat("\n[Step 11] Epithelial labels — Malignant vs Keratinocyte_like ...\n")

delta_thr <- 0.10

markers_epith <- list(
  A_name = "Malignant",
  B_name = "Keratinocyte_like",
  A_up   = c("MKI67","TOP2A","UBE2C","PCNA","TYMS","BIRC5","CCNB1","AURKA","AURKB","KIF11","TK1","H2AFZ","EPCAM","VIM"),
  B_up   = c("KRT5","KRT14","KRT15","KRT1","KRT10","DSG3","IVL","SPRR1A","SPRR1B","TP63","KRT6A","KRT6B")
)

say_msg <- function(x) {
  message("\n", x, "\n")
}

norm_cluster_name <- function(x) {
  x <- as.character(x)
  sub("^[^0-9-]*", "", x)
}

avg_by_cluster <- function(obj, assay = "RNA", slot = "data") {
  if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
    stop("seurat_clusters is not present in metadata.")
  }
  ae <- suppressWarnings(
    try(
      AverageExpression(obj, group.by = "seurat_clusters", assays = assay, slot = slot),
      silent = TRUE
    )
  )
  if (inherits(ae, "try-error") || is.null(ae[[assay]])) {
    stop("AverageExpression failed.")
  }
  as.matrix(ae[[assay]])
}

match_genes_to_ref <- function(genes, ref) {
  idx <- match(toupper(genes), toupper(ref))
  ref[na.omit(idx)]
}

module_score_avg <- function(avg_expr, gene_set) {
  gs <- match_genes_to_ref(gene_set, rownames(avg_expr))
  if (length(gs) == 0) {
    return(rep(NA_real_, ncol(avg_expr)))
  }
  colMeans(avg_expr[gs, , drop = FALSE])
}

auto_pick_clusters <- function(obj, markers, delta = 0.0, top_k_each = NULL) {
  avg <- avg_by_cluster(obj)
  cl_raw <- colnames(avg)
  cl     <- norm_cluster_name(cl_raw)
  
  sA <- module_score_avg(avg, markers$A_up)
  sB <- module_score_avg(avg, markers$B_up)
  
  d <- sA - sB
  names(d) <- cl
  
  A <- names(d[d >=  delta])
  B <- names(d[d <= -delta])
  
  ## Remove overlaps
  A <- setdiff(A, B)
  B <- setdiff(B, A)
  
  if (length(A) == 0 || length(B) == 0) {
    ord <- sort(d, decreasing = TRUE)
    k <- if (is.null(top_k_each)) max(1, floor(length(ord) / 2)) else top_k_each
    A <- names(head(ord, k))
    B <- setdiff(names(tail(ord, k)), A)
  }
  
  A_final <- intersect(norm_cluster_name(levels(obj$seurat_clusters)), A)
  B_final <- intersect(norm_cluster_name(levels(obj$seurat_clusters)), B)
  
  list(
    A      = A_final,
    B      = B_final,
    scores = data.frame(
      cluster = cl,
      score_A = sA,
      score_B = sB,
      diff    = as.numeric(d[cl]),
      check.names = FALSE
    )
  )
}

write_labels_epith <- function(obj, picks, nameA, nameB, suffix) {
  grp <- ifelse(
    norm_cluster_name(obj$seurat_clusters) %in% picks$A, nameA,
    ifelse(norm_cluster_name(obj$seurat_clusters) %in% picks$B, nameB, NA)
  )
  
  df <- data.frame(
    cell_barcode   = colnames(obj),
    seurat_cluster = as.character(obj$seurat_clusters),
    assigned_group = grp,
    check.names    = FALSE
  )
  
  df <- df[!is.na(df$assigned_group), , drop = FALSE]
  
  fwrite(
    df,
    file.path(data_dir, paste0(project_id, "_", suffix, "_cells.csv"))
  )
  
  scores_df <- picks$scores
  scores_df$prelabel <- ifelse(
    norm_cluster_name(scores_df$cluster) %in% picks$A, nameA,
    ifelse(norm_cluster_name(scores_df$cluster) %in% picks$B, nameB, "none")
  )
  
  fwrite(
    scores_df,
    file.path(data_dir, paste0(project_id, "_", suffix, "_autoSelection.csv"))
  )
  
  invisible(df)
}

say_msg("##\n[Step 11] Epithelial labels — Malignant vs Keratinocyte_like\n##")

if (!exists("epith_obj") || is.null(epith_obj)) {
  stop("epith_obj does not exist. Cannot run Step 11.")
}

message("  seurat_clusters levels (epith_obj):")
print(levels(epith_obj$seurat_clusters))

picks_ep <- tryCatch(
  auto_pick_clusters(epith_obj, markers_epith, delta = delta_thr),
  error = function(e) {
    stop("Automatic selection of epithelial clusters failed: ", e$message)
  }
)

message("  Malignant clusters (auto): {", paste(picks_ep$A, collapse = ", "), "}")
message("  Keratinocyte_like clusters (auto): {", paste(picks_ep$B, collapse = ", "), "}")

write_labels_epith(
  epith_obj,
  picks_ep,
  "Malignant",
  "Keratinocyte_like",
  "Stage11_Epithelial_Malignant_vs_KeratinocyteLike"
)

say_msg("Step 11 finished — epithelial cell labels & autoSelection tables saved.")

############################################################
# STEP 12 (scRNA) — Epithelial DEGs (Malignant vs Keratinocyte_like)
############################################################

cat("\n[Step 12] Epithelial DEGs — Malignant vs Keratinocyte_like (strict) ...\n")

min_cells_per_group <- 50
logfc_thr_strict    <- 0.50
minpct_strict       <- 0.20
padj_thr_strict     <- 0.05
test_use_method     <- "wilcox"
only_pos_genes      <- FALSE

apply_labels_from_csv <- function(obj,
                                  labels_csv,
                                  id_col    = "cell_barcode",
                                  group_col = "assigned_group",
                                  ident_name = ".__label__") {
  if (!file.exists(labels_csv)) {
    message("  File does not exist: ", labels_csv)
    return(NULL)
  }
  
  lab <- tryCatch(fread(labels_csv), error = function(e) NULL)
  if (is.null(lab) || !(id_col %in% names(lab)) || !(group_col %in% names(lab))) {
    message("  Required columns not found in ", labels_csv)
    return(NULL)
  }
  
  setnames(lab, c(id_col, group_col), c("cell_barcode", "assigned_group"), skip_absent = TRUE)
  lab <- lab[lab$cell_barcode %in% colnames(obj), ]
  if (nrow(lab) == 0) {
    message("  No overlapping cells between labels and Seurat object.")
    return(NULL)
  }
  
  md <- obj@meta.data
  md$.__cb__ <- rownames(md)
  md <- dplyr::left_join(md, as.data.frame(lab), by = c(".__cb__" = "cell_barcode"))
  rownames(md) <- md$.__cb__
  md$.__cb__   <- NULL
  obj@meta.data <- md
  
  obj[[ident_name]] <- obj@meta.data$assigned_group
  
  keep_cells <- rownames(obj@meta.data)[!is.na(obj@meta.data[[ident_name]])]
  sub <- subset(obj, cells = keep_cells)
  if (ncol(sub) == 0) {
    message("  Empty subset after label assignment.")
    return(NULL)
  }
  
  Idents(sub) <- ident_name
  sub
}

run_deg_strict <- function(obj,
                           g1     = "Malignant",
                           g2     = "Keratinocyte_like",
                           suffix = "Stage12_Epithelial_Malignant_vs_KeratinocyteLike_DEGs") {
  n1 <- sum(Idents(obj) == g1)
  n2 <- sum(Idents(obj) == g2)
  
  message("  Number of cells — ", g1, ": ", n1, " | ", g2, ": ", n2)
  
  if (n1 < min_cells_per_group || n2 < min_cells_per_group) {
    message("  Insufficient cells (minimum ", min_cells_per_group, ") for ", suffix, ". Skipping.")
    return(NULL)
  }
  
  fwrite(
    data.table(group = c(g1, g2), n_cells = c(n1, n2)),
    file.path(data_dir, paste0(project_id, "_", suffix, "_group_counts.csv"))
  )
  
  message("  Running FindMarkers (", g1, " vs ", g2, ") ...")
  
  deg <- tryCatch(
    {
      FindMarkers(
        obj,
        ident.1         = g1,
        ident.2         = g2,
        only.pos        = only_pos_genes,
        logfc.threshold = logfc_thr_strict,
        min.pct         = minpct_strict,
        test.use        = test_use_method,
        return.thresh   = 1
      )
    },
    error = function(e) {
      message("  Error in FindMarkers: ", e$message)
      NULL
    }
  )
  
  if (is.null(deg) || nrow(deg) == 0) {
    message("  DEG result is empty for ", suffix)
    return(NULL)
  }
  
  tbl <- data.frame(
    gene = rownames(deg),
    deg,
    row.names   = NULL,
    check.names = FALSE
  )
  
  ## Direction based on logFC/log2FC
  fc_col <- intersect(c("avg_log2FC", "avg_logFC", "avg_log2fc"), colnames(tbl))[1]
  if (!is.na(fc_col)) {
    tbl$direction <- ifelse(tbl[[fc_col]] > 0, "up", "down")
  }
  
  ## Filter by p_val_adj < padj_thr_strict
  padj_col <- intersect(c("p_val_adj", "p_val"), colnames(tbl))[1]
  if (!is.na(padj_col)) {
    tbl_filt <- tbl[which(tbl[[padj_col]] < padj_thr_strict), , drop = FALSE]
  } else {
    tbl_filt <- tbl
  }
  
  if (nrow(tbl_filt) == 0) {
    message("  No genes passed padj <", padj_thr_strict, " for ", suffix)
    return(NULL)
  }
  
  out_path <- file.path(data_dir, paste0(project_id, "_", suffix, ".csv"))
  fwrite(tbl_filt, out_path)
  message("  Saved epithelial DEGs (padj <", padj_thr_strict, ") to: ", out_path)
  
  tbl_filt
}

say_msg("##\n[Step 12] Epithelial DEGs — Malignant vs Keratinocyte_like (strict)\n##")

lab_ep <- file.path(
  data_dir,
  paste0(project_id, "_Stage11_Epithelial_Malignant_vs_KeratinocyteLike_cells.csv")
)

epi_sub <- apply_labels_from_csv(epith_obj, lab_ep)

if (!is.null(epi_sub)) {
  groups <- levels(Idents(epi_sub))
  g1 <- if ("Malignant" %in% groups) "Malignant" else groups[1]
  g2 <- if ("Keratinocyte_like" %in% groups) "Keratinocyte_like" else setdiff(groups, g1)[1]
  run_deg_strict(
    epi_sub,
    g1,
    g2,
    "Stage12_Epithelial_Malignant_vs_KeratinocyteLike_DEGs"
  )
}

say_msg("Step 12 finished — epithelial DEG table and group counts saved.")
