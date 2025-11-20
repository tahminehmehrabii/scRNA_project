
############################################################
# Full single-cell pipeline for LSCC (GSE150321) - STEPWISE
# Each stage saves CSV / image / RDS outputs
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

############################################################
# 0. Paths & basic settings
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  data_dir <- args[1]
} else {
  data_dir <- "D:/Single Cell/GSE150321_RAW"
}

if (length(args) >= 2) {
  project_id <- args[2]
} else {
  project_id <- "LSCC_GSE150321"
}

cat("Data directory:", data_dir, "\n")
cat("Project ID    :", project_id, "\n\n")

setwd(data_dir)

############################################################
# Stage 1 — Read raw CSV.gz files (per sample)
############################################################

raw_files <- list.files(data_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)
cat("Found files:\n")
print(raw_files)

if (length(raw_files) == 0) {
  stop("No .csv.gz files found in the folder. Please check the data_dir path.")
}

read_one_sample <- function(file_path, sample_id = NULL) {
  if (is.null(sample_id)) {
    sample_id <- gsub("\\.csv\\.gz$", "", basename(file_path))
  }
  message("Reading: ", file_path)
  dt <- fread(file_path)
  expr_mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(expr_mat) <- dt[[1]]
  if (any(duplicated(rownames(expr_mat)))) {
    rownames(expr_mat) <- make.unique(rownames(expr_mat))
  }
  seu <- CreateSeuratObject(
    counts       = expr_mat,
    project      = sample_id,
    min.cells    = 1,
    min.features = 1
  )
  seu$sample_id <- sample_id
  return(seu)
}

seu_list <- list()
for (f in raw_files) {
  sid <- gsub("\\.csv\\.gz$", "", basename(f))
  seu_list[[sid]] <- read_one_sample(f, sample_id = sid)
}

cat("\n[Stage1] Number of cells before QC:\n")
n_cells_raw <- sapply(seu_list, ncol)
print(n_cells_raw)

# Save raw summary table
stage1_qc <- data.frame(
  sample_id = names(n_cells_raw),
  n_cells   = as.numeric(n_cells_raw),
  stage     = "raw_before_QC"
)
write.csv(stage1_qc,
          file = file.path(data_dir, paste0(project_id, "_Stage1_raw_cell_counts.csv")),
          row.names = FALSE)

############################################################
# Stage 2 — QC per sample (nFeature, nCount, percent.mt)
############################################################

minFeature <- 200
maxFeature <- 7500
minCount   <- 400
maxCount   <- 40000
maxMT      <- 20

seu_list <- lapply(seu_list, function(seu) {
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  subset(
    seu,
    subset = nFeature_RNA > minFeature &
      nFeature_RNA < maxFeature &
      nCount_RNA    > minCount &
      nCount_RNA    < maxCount &
      percent.mt    < maxMT
  )
})

cat("\n[Stage2] Number of cells after QC:\n")
n_cells_qc <- sapply(seu_list, ncol)
print(n_cells_qc)

stage2_qc <- data.frame(
  sample_id = names(n_cells_qc),
  n_cells   = as.numeric(n_cells_qc),
  stage     = "post_QC"
)
write.csv(stage2_qc,
          file = file.path(data_dir, paste0(project_id, "_Stage2_postQC_cell_counts.csv")),
          row.names = FALSE)

############################################################
# Stage 3 — Merge samples
############################################################

first_name <- names(seu_list)[1]
lscc <- seu_list[[first_name]]
if (length(seu_list) > 1) {
  other_names <- names(seu_list)[-1]
  lscc <- merge(
    x            = lscc,
    y            = seu_list[other_names],
    add.cell.ids = names(seu_list),
    project      = project_id
  )
}

cat("\n[Stage3] Merged object:\n")
print(lscc)
cat("\nCells per sample_id:\n")
print(table(lscc$sample_id))

# Save metadata after merge
meta_stage3 <- lscc@meta.data
write.csv(meta_stage3,
          file = file.path(data_dir, paste0(project_id, "_Stage3_meta_merged.csv")),
          row.names = TRUE)

# Seurat v5: join layers
lscc <- JoinLayers(lscc)
DefaultAssay(lscc) <- "RNA"

# Save merged Seurat object as RDS
saveRDS(lscc, file = file.path(data_dir, paste0(project_id, "_Stage3_merged_seurat.rds")))

############################################################
# Stage 4 — Normalize, HVG, Scale, PCA
############################################################

lscc <- NormalizeData(lscc)
lscc <- FindVariableFeatures(lscc, selection.method = "vst", nfeatures = 3000)
cat("\n[Stage4] Number of variable features:", length(VariableFeatures(lscc)), "\n")

# Save list of variable genes
write.csv(
  data.frame(gene = VariableFeatures(lscc)),
  file = file.path(data_dir, paste0(project_id, "_Stage4_HVG_genes.csv")),
  row.names = FALSE
)

lscc <- ScaleData(lscc, features = VariableFeatures(lscc))
lscc <- RunPCA(lscc, features = VariableFeatures(lscc))

# Elbow plot
p_elbow <- ElbowPlot(lscc, ndims = 50) +
  ggtitle("Stage4 - PCA elbow plot")
ggsave(file.path(data_dir, paste0(project_id, "_Stage4_PCA_Elbow.png")),
       plot = p_elbow, width = 6, height = 5, dpi = 1200)

# Save Seurat object after PCA
saveRDS(lscc, file = file.path(data_dir, paste0(project_id, "_Stage4_postPCA_seurat.rds")))

############################################################
# Stage 5 — Neighbor graph, Clustering, UMAP
############################################################

dims_use <- 1:30
lscc <- FindNeighbors(lscc, dims = dims_use)
lscc <- FindClusters(lscc, resolution = 0.5)
lscc <- RunUMAP(lscc, dims = dims_use)

# celltype_main initially equals cluster IDs
lscc$celltype_main <- as.character(lscc$seurat_clusters)

# UMAP by clusters
p_umap_cluster <- DimPlot(lscc, reduction = "umap",
                          group.by = "seurat_clusters",
                          label = TRUE, repel = TRUE) +
  ggtitle("Stage5 - Global UMAP - clusters") +
  NoLegend()

# UMAP by sample
p_umap_sample <- DimPlot(lscc, reduction = "umap",
                         group.by = "sample_id") +
  ggtitle("Stage5 - Global UMAP - samples")

ggsave(file.path(data_dir, paste0(project_id, "_Stage5_UMAP_clusters.png")),
       plot = p_umap_cluster, width = 7, height = 6, dpi = 1200)
ggsave(file.path(data_dir, paste0(project_id, "_Stage5_UMAP_samples.png")),
       plot = p_umap_sample, width = 7, height = 6, dpi = 1200)

# Save metadata and Seurat object
write.csv(lscc@meta.data,
          file = file.path(data_dir, paste0(project_id, "_Stage5_meta_with_clusters.csv")),
          row.names = TRUE)
saveRDS(lscc, file = file.path(data_dir, paste0(project_id, "_Stage5_postUMAP_seurat.rds")))

############################################################
# Stage 6 — Global markers per cluster (UP + DOWN)
############################################################

all_markers <- FindAllMarkers(
  lscc,
  only.pos        = FALSE,   # both up and down
  min.pct         = 0.25,
  logfc.threshold = 0.25
)

# Add expression direction column
all_markers$direction <- ifelse(all_markers$avg_log2FC > 0, "up", "down")

# Save full marker table
write.csv(
  all_markers,
  file = file.path(data_dir, paste0(project_id, "_Stage6_cluster_markers_global_up_down.csv")),
  row.names = FALSE
)

############################################################
# Stage 7 — Canonical + LSCC-relevant lineage DotPlot
############################################################

## T cells
t_markers <- c("CD3D", "CD3E", "CD2", "TRAC", "CD4", "CD8A", "CD8B")

## B / Plasma cells
b_markers <- c("MS4A1", "CD79A", "CD79B", "CD19", "CD74", "MZB1", "SDC1")

## Myeloid (general)
myeloid_markers <- c("LYZ", "CD68", "CD163", "LST1", "CSF1R", "ITGAM", "FCGR3A")

## Tumor / epithelial (base)
epith_tumor_markers <- c("EPCAM", "KRT8", "KRT18", "KRT5", "KRT14", "KRT17", "TP63")

## Extra squamous / HNSCC epithelial markers
epith_hnscc_extra <- c("KRT4", "KRT6A", "KRT6B", "KRT16", "KRT1", "KRT10", "KRT13")

## Fibroblasts / CAF
fibro_markers <- c("COL1A1", "COL1A2", "COL3A1", "TAGLN", "ACTA2", "PDGFRA", "DCN")

## Endothelial
endo_markers <- c("PECAM1", "CDH5", "VWF", "KDR", "ENG", "RGS5")

## Proliferation
prolif_markers <- c("MKI67", "TOP2A", "PCNA", "BIRC5")

## Mast cells (very relevant in HNSCC)
mast_markers <- c("TPSAB1", "TPSB2", "CPA3")

## NK / cytotoxic T
nk_markers <- c("NKG7", "GNLY", "PRF1", "GZMB", "CTSW")

## TAM / SPP1+ macrophages (optional but useful)
tam_markers <- c("MRC1", "CD163", "SPP1", "APOE")

## Finalize marker panel
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

## Keep only genes that actually exist in the Seurat object
marker_panel <- marker_panel[marker_panel %in% rownames(lscc)]

## Save list of markers used
write.csv(
  data.frame(marker = marker_panel),
  file = file.path(data_dir, paste0(project_id, "_Stage7_marker_panel_used.csv")),
  row.names = FALSE
)

## Draw DotPlot
if (length(marker_panel) > 0) {
  p_dot <- DotPlot(lscc, features = marker_panel) +
    RotatedAxis() +
    ggtitle("Stage7 - Canonical & LSCC lineage markers across clusters")
  
  ggsave(
    file.path(data_dir, paste0(project_id, "_Stage7_lineage_DotPlot.png")),
    plot = p_dot, width = 10, height = 6, dpi = 1200
  )
}



############################################################
# Stage 8 — Manual mapping: cluster -> major cell type
############################################################

cat("\n[Stage8] Current cluster levels:\n")
print(levels(lscc))

# This is just an example; you can change it later if needed
cluster_to_celltype <- c(
  "Tumor/Epithelial",  # 0
  "T cells",           # 1
  "Myeloid",           # 2
  "B/Plasma cells",    # 3
  "Fibroblasts",       # 4
  "Endothelial",       # 5
  "Tumor/Epithelial",  # 6
  "T cells",           # 7
  "Myeloid",           # 8
  "Tumor/Epithelial",  # 9
  "B/Plasma cells",    # 10
  "Fibroblasts",       # 11
  "Myeloid",           # 12
  "Tumor/Epithelial",  # 13
  "T cells"            # 14
)

if (length(cluster_to_celltype) != length(levels(lscc))) {
  warning("Stage8: Length of cluster_to_celltype does not match number of clusters. celltype_main will remain as cluster IDs for now.")
} else {
  names(cluster_to_celltype) <- levels(lscc)
  lscc <- RenameIdents(lscc, cluster_to_celltype)
  lscc$celltype_main <- Idents(lscc)
}

# UMAP with celltype_main
p_umap_ct <- DimPlot(lscc, reduction = "umap",
                     group.by = "celltype_main", label = TRUE, repel = TRUE) +
  ggtitle("Stage8 - Global UMAP - main cell types")
ggsave(file.path(data_dir, paste0(project_id, "_Stage8_UMAP_celltypes.png")),
       plot = p_umap_ct, width = 7, height = 6, dpi = 1200)

# Save metadata after annotation
write.csv(lscc@meta.data,
          file = file.path(data_dir, paste0(project_id, "_Stage8_meta_with_celltypes.csv")),
          row.names = TRUE)

############################################################
# Stage 9 — Per-lineage subsetting & reclustering
############################################################

run_lineage_pipeline <- function(obj, subset_ids, lineage_name,
                                 dims_use = 1:20, res = 0.5, project_id, data_dir) {
  message("Processing lineage: ", lineage_name)
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
  sub <- RunPCA(sub, features = VariableFeatures(sub))
  sub <- FindNeighbors(sub, dims = dims_use)
  sub <- FindClusters(sub, resolution = res)
  sub <- RunUMAP(sub, dims = dims_use)
  
  # UMAP plots
  p1 <- DimPlot(sub, reduction = "umap", group.by = "seurat_clusters",
                label = TRUE, repel = TRUE) +
    ggtitle(paste0("Stage9 - ", lineage_name, " UMAP - clusters")) +
    NoLegend()
  p2 <- DimPlot(sub, reduction = "umap", group.by = "sample_id") +
    ggtitle(paste0("Stage9 - ", lineage_name, " UMAP - samples"))
  
  ggsave(file.path(data_dir, paste0(project_id, "_Stage9_", lineage_name, "_UMAP_clusters.png")),
         plot = p1, width = 6, height = 5, dpi = 1200)
  ggsave(file.path(data_dir, paste0(project_id, "_Stage9_", lineage_name, "_UMAP_samples.png")),
         plot = p2, width = 6, height = 5, dpi = 1200)
  
  # markers
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
  
  # Save metadata and RDS
  write.csv(
    sub@meta.data,
    file = file.path(data_dir, paste0(project_id, "_Stage9_", lineage_name, "_meta.csv")),
    row.names = TRUE
  )
  saveRDS(sub,
          file = file.path(data_dir, paste0(project_id, "_Stage9_", lineage_name, "_seurat.rds")))
  return(sub)
}

# These names must match the final celltype_main values
t_lineage_ids       <- c("T cells")
b_lineage_ids       <- c("B/Plasma cells")
myeloid_lineage_ids <- c("Myeloid")
epith_lineage_ids   <- c("Tumor/Epithelial")
fibro_lineage_ids   <- c("Fibroblasts")
endo_lineage_ids    <- c("Endothelial")

t_obj       <- run_lineage_pipeline(lscc, t_lineage_ids,       "T_cells",      project_id = project_id, data_dir = data_dir)
b_obj       <- run_lineage_pipeline(lscc, b_lineage_ids,       "B_cells",      project_id = project_id, data_dir = data_dir)
myeloid_obj <- run_lineage_pipeline(lscc, myeloid_lineage_ids, "Myeloid",      project_id = project_id, data_dir = data_dir)
epith_obj   <- run_lineage_pipeline(lscc, epith_lineage_ids,   "EpithelialTumor", project_id = project_id, data_dir = data_dir)
fibro_obj   <- run_lineage_pipeline(lscc, fibro_lineage_ids,   "Fibroblasts",  project_id = project_id, data_dir = data_dir)
endo_obj    <- run_lineage_pipeline(lscc, endo_lineage_ids,    "Endothelial",  project_id = project_id, data_dir = data_dir)

############################################################
# Stage 10 — Save final global object
############################################################

saveRDS(lscc,
        file = file.path(data_dir, paste0(project_id, "_Stage10_Global_annotated_seurat.rds")))
cat("\nGlobal annotated Seurat object saved to:",
    file.path(data_dir, paste0(project_id, "_Stage10_Global_annotated_seurat.rds")), "\n")

cat("\nPipeline finished.\n")


############################################################
# Stage 11 — Epithelial labels (Malignant vs Keratinocyte_like, NO DEG)
# Stage 12 — Epithelial DEGs (Malignant vs Keratinocyte_like, stringent)
# Prerequisites: epith_obj, data_dir, and project_id must exist from previous stages
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(data.table)
})

## ====================== Parameters ======================
# Stage 11 (automatic cluster selection using module scores)
delta_thr <- 0.10

# Markers to distinguish malignant vs keratinocyte-like epithelial cells
markers_epith <- list(
  A_name = "Malignant",
  B_name = "Keratinocyte_like",
  A_up   = c("MKI67","TOP2A","UBE2C","PCNA","TYMS","BIRC5","CCNB1","AURKA","AURKB","KIF11","TK1","H2AFZ","EPCAM","VIM"),
  B_up   = c("KRT5","KRT14","KRT15","KRT1","KRT10","DSG3","IVL","SPRR1A","SPRR1B","TP63","KRT6A","KRT6B")
)

# Stage 12 (strict DEG analysis)
min_cells_per_group <- 50
logfc_thr_strict    <- 0.50
minpct_strict       <- 0.20
test_use_method     <- "wilcox"
only_pos_genes      <- FALSE   # both up and down

## ====================== Helper functions ======================
.say <- function(x) message("\n", x, "\n")
.norm <- function(x){ x <- as.character(x); sub("^[^0-9-]*","",x) }

.avg_by_cluster <- function(obj, assay="RNA", slot="data"){
  if (!"seurat_clusters" %in% colnames(obj@meta.data))
    stop("seurat_clusters is not present in metadata.")
  ae <- suppressWarnings(try(AverageExpression(obj, group.by="seurat_clusters", assays=assay, slot=slot), silent=TRUE))
  if (inherits(ae,"try-error") || is.null(ae[[assay]]))
    stop("AverageExpression failed.")
  as.matrix(ae[[assay]])
}

.match_genes <- function(genes, ref){
  idx <- match(toupper(genes), toupper(ref))
  ref[na.omit(idx)]
}

.module_score <- function(avg_expr, gene_set){
  gs <- .match_genes(gene_set, rownames(avg_expr))
  if (length(gs)==0) return(rep(NA_real_, ncol(avg_expr)))
  colMeans(avg_expr[gs, , drop=FALSE])
}

.auto_pick_clusters <- function(obj, markers, delta=0.0, top_k_each=NULL){
  avg <- .avg_by_cluster(obj)
  cl_raw <- colnames(avg); cl <- .norm(cl_raw)
  sA <- .module_score(avg, markers$A_up)
  sB <- .module_score(avg, markers$B_up)
  d  <- sA - sB; names(d) <- cl
  
  A <- names(d[d >=  delta])
  B <- names(d[d <= -delta])
  A <- setdiff(A,B); B <- setdiff(B,A)
  
  if (length(A)==0 || length(B)==0){
    ord <- sort(d, decreasing=TRUE)
    k <- if (is.null(top_k_each)) max(1, floor(length(ord)/2)) else top_k_each
    A <- names(head(ord, k))
    B <- setdiff(names(tail(ord, k)), A)
  }
  
  list(
    A = intersect(.norm(levels(obj$seurat_clusters)), A),
    B = intersect(.norm(levels(obj$seurat_clusters)), B),
    scores = data.frame(cluster=cl, score_A=sA, score_B=sB, diff=as.numeric(d[cl]), check.names=FALSE)
  )
}

.write_labels_epith <- function(obj, picks, nameA, nameB, suffix){
  grp <- ifelse(.norm(obj$seurat_clusters) %in% picks$A, nameA,
                ifelse(.norm(obj$seurat_clusters) %in% picks$B, nameB, NA))
  df <- data.frame(
    cell_barcode   = colnames(obj),
    seurat_cluster = as.character(obj$seurat_clusters),
    assigned_group = grp,
    check.names = FALSE
  )
  df <- df[!is.na(df$assigned_group), , drop=FALSE]
  # Save
  fwrite(df, file.path(data_dir, paste0(project_id, "_", suffix, "_cells.csv")))
  fwrite(picks$scores, file.path(data_dir, paste0(project_id, "_", suffix, "_autoSelection.csv")))
  invisible(df)
}

.apply_labels_from_csv <- function(obj, labels_csv,
                                   id_col="cell_barcode", group_col="assigned_group",
                                   ident_name=".__label__"){
  if (!file.exists(labels_csv)) { message(" File not found: ", labels_csv); return(NULL) }
  lab <- tryCatch(fread(labels_csv), error=function(e) NULL)
  if (is.null(lab) || !(id_col %in% names(lab)) || !(group_col %in% names(lab))) {
    message(" Required columns were not found in ", labels_csv, "."); return(NULL)
  }
  setnames(lab, c(id_col, group_col), c("cell_barcode","assigned_group"), skip_absent=TRUE)
  lab <- lab[lab$cell_barcode %in% colnames(obj), ]
  if (nrow(lab)==0) { message(" No cell overlap with object for: ", labels_csv); return(NULL) }
  
  md <- obj@meta.data
  md$.__cb__ <- rownames(md)
  md <- md |> dplyr::left_join(as.data.frame(lab), by=c(".__cb__"="cell_barcode"))
  rownames(md) <- md$.__cb__; md$.__cb__ <- NULL
  obj@meta.data <- md
  
  obj[[ident_name]] <- obj@meta.data$assigned_group
  sub <- subset(obj, subset = !is.na(!!as.name(ident_name)))
  if (ncol(sub)==0) { message(" Subset is empty after labeling."); return(NULL) }
  Idents(sub) <- ident_name
  sub
}

.run_deg_strict <- function(obj, g1="Malignant", g2="Keratinocyte_like",
                            suffix="Stage12_Epithelial_Malignant_vs_KeratinocyteLike_DEGs"){
  n1 <- sum(Idents(obj)==g1); n2 <- sum(Idents(obj)==g2)
  message("  Number of cells — ", g1, ": ", n1, " | ", g2, ": ", n2)
  if (n1 < min_cells_per_group || n2 < min_cells_per_group) {
    message(" Not enough cells (minimum ", min_cells_per_group, "). Skipping: ", suffix)
    return(NULL)
  }
  # Save group counts
  fwrite(data.table(group=c(g1,g2), n_cells=c(n1,n2)),
         file.path(data_dir, paste0(project_id, "_", suffix, "_group_counts.csv")))
  
  message("  Running FindMarkers (", g1, " vs ", g2, ") ...")
  deg <- tryCatch({
    FindMarkers(obj,
                ident.1         = g1,
                ident.2         = g2,
                only.pos        = only_pos_genes,
                logfc.threshold = logfc_thr_strict,
                min.pct         = minpct_strict,
                test.use        = test_use_method
    )
  }, error=function(e){ message(" Error in FindMarkers: ", e$message); NULL })
  
  if (is.null(deg) || nrow(deg)==0) { message(" DEG result is empty: ", suffix); return(NULL) }
  
  tbl <- data.frame(gene=rownames(deg), deg, row.names=NULL, check.names=FALSE)
  fc_col <- intersect(c("avg_log2FC","avg_logFC","avg_log2fc"), colnames(tbl))[1]
  if (!is.na(fc_col)) tbl$direction <- ifelse(tbl[[fc_col]] > 0, "up", "down")
  
  out_path <- file.path(data_dir, paste0(project_id, "_", suffix, ".csv"))
  fwrite(tbl, out_path)
  message(" Saved: ", out_path)
  tbl
}

## ====================== Stage 11 ======================
.say("############################################################\n[Stage11] Epithelial labels — Malignant vs Keratinocyte_like\n############################################################")
if (!exists("epith_obj") || is.null(epith_obj)) {
  stop("epith_obj does not exist. Cannot run Stage11/12.")
}
message("  seurat_clusters levels (epith_obj):"); print(levels(epith_obj$seurat_clusters))

picks_ep <- tryCatch(
  .auto_pick_clusters(epith_obj, markers_epith, delta=delta_thr),
  error=function(e){ stop("Automatic selection of epithelial clusters failed: ", e$message) }
)
message("  Malignant clusters: {", paste(picks_ep$A, collapse=", "), "}")
message("  Keratinocyte_like clusters: {", paste(picks_ep$B, collapse=", "), "}")

# Save cell labels and cluster score table
.write_labels_epith(epith_obj, picks_ep, "Malignant","Keratinocyte_like",
                    "Stage11_Epithelial_Malignant_vs_KeratinocyteLike")

.say("Stage11 finished — the following files were saved:\n - ..._Stage11_Epithelial_Malignant_vs_KeratinocyteLike_cells.csv\n - ..._Stage11_Epithelial_Malignant_vs_KeratinocyteLike_autoSelection.csv")

## ====================== Stage 12 ======================
.say("############################################################\n[Stage12] Epithelial DEGs — Malignant vs Keratinocyte_like (strict)\n############################################################")

# Apply Stage 11 labels and run strict DEG
lab_ep <- file.path(data_dir, paste0(project_id, "_Stage11_Epithelial_Malignant_vs_KeratinocyteLike_cells.csv"))
epi_sub <- .apply_labels_from_csv(epith_obj, lab_ep)
if (!is.null(epi_sub)) {
  # If exact group names are missing, fall back to closest available levels
  groups <- levels(Idents(epi_sub))
  g1 <- if ("Malignant" %in% groups) "Malignant" else groups[1]
  g2 <- if ("Keratinocyte_like" %in% groups) "Keratinocyte_like" else setdiff(groups, g1)[1]
  .run_deg_strict(epi_sub, g1, g2, "Stage12_Epithelial_Malignant_vs_KeratinocyteLike_DEGs")
}

.say("Stage12 finished — DEG table and group counts were saved.")
