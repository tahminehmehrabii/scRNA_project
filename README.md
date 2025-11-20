# Single-cell RNA-seq  Project
## LSCC scRNA-seq Pipeline Overview (Seurat, R)

This repository contains a fully scripted single-cell RNA-seq pipeline for LSCC (GSE150321) implemented in R/Seurat.  
Each stage saves intermediate CSV/PNG/RDS outputs for reproducibility.

### Stage 0 – Setup
- Parse command-line arguments (`data_dir`, `project_id`) and set working directory.
- Load required libraries (Seurat, data.table, dplyr, ggplot2).

### Stage 1 – Load raw data
- Scan `data_dir` for per-sample `.csv.gz` files.
- For each sample: read count matrix, create `Seurat` object, attach `sample_id`.
- Save raw cell counts per sample (`*_Stage1_raw_cell_counts.csv`).

### Stage 2 – Quality control (per sample)
- Compute `nFeature_RNA`, `nCount_RNA`, and mitochondrial percentage (`percent.mt`).
- Filter cells using hard thresholds (e.g. 200–7500 genes, 400–40000 UMIs, <20% MT).
- Save post-QC cell counts per sample (`*_Stage2_postQC_cell_counts.csv`).

### Stage 3 – Merge samples
- Merge all QC-filtered Seurat objects into one multi-sample object.
- Join layers (Seurat v5), set `RNA` as default assay.
- Save merged metadata and Seurat object:
  - `*_Stage3_meta_merged.csv`
  - `*_Stage3_merged_seurat.rds`

### Stage 4 – Normalization & PCA
- Normalize data and select highly variable genes (HVGs, `vst`, 3000 genes).
- Save HVG list (`*_Stage4_HVG_genes.csv`).
- Scale data on HVGs and run PCA.
- Export PCA elbow plot (`*_Stage4_PCA_Elbow.png`) and post-PCA object (`*_Stage4_postPCA_seurat.rds`).

### Stage 5 – Graph clustering & UMAP
- Build kNN graph (`FindNeighbors`) and perform clustering (`FindClusters`, res = 0.5).
- Run UMAP on selected PCs.
- Generate UMAPs:
  - By cluster (`*_Stage5_UMAP_clusters.png`)
  - By sample (`*_Stage5_UMAP_samples.png`)
- Save metadata with clusters and Seurat object:
  - `*_Stage5_meta_with_clusters.csv`
  - `*_Stage5_postUMAP_seurat.rds`

### Stage 6 – Cluster marker genes
- Run `FindAllMarkers` across clusters (up- and down-regulated genes).
- Add direction column based on log2FC (up vs down).
- Save global marker table (`*_Stage6_cluster_markers_global_up_down.csv`).

### Stage 7 – Canonical lineage markers (DotPlot)
- Define curated marker panels for:
  - T cells, B/plasma cells, myeloid cells, tumor/epithelial, fibroblasts/CAF, endothelial, proliferation, mast, NK/cytotoxic, TAM/SPP1+ macrophages.
- Filter to markers present in the dataset and save panel (`*_Stage7_marker_panel_used.csv`).
- Plot DotPlot of marker expression across clusters (`*_Stage7_lineage_DotPlot.png`).

### Stage 8 – Map clusters to major cell types
- Manually map `seurat_clusters` to major lineages:
  - Tumor/Epithelial, T cells, B/Plasma cells, Myeloid, Fibroblasts, Endothelial.
- Store mapping in `celltype_main` and visualize UMAP by cell type:
  - `*_Stage8_UMAP_celltypes.png`
- Save updated metadata (`*_Stage8_meta_with_celltypes.csv`).

### Stage 9 – Lineage-specific sub-pipelines
- For each major lineage (T cells, B cells, myeloid, epithelial tumor, fibroblasts, endothelial):
  - Subset cells of that lineage.
  - Re-run normalization, HVG selection, scaling, PCA, neighbors, clustering, and UMAP.
  - Export:
    - UMAP by cluster and sample (`*_Stage9_<Lineage>_UMAP_*.png`)
    - Lineage-specific marker genes (`*_Stage9_<Lineage>_markers.csv`)
    - Lineage metadata (`*_Stage9_<Lineage>_meta.csv`)
    - Lineage Seurat object (`*_Stage9_<Lineage>_seurat.rds`).

### Stage 10 – Global annotated object
- Save final global Seurat object with cluster and cell-type annotations:
  - `*_Stage10_Global_annotated_seurat.rds`.

### Stage 11 – Epithelial labeling (Malignant vs Keratinocyte-like)
- Use predefined marker sets for malignant and keratinocyte-like programs.
- Compute module scores per cluster and automatically assign clusters to:
  - `Malignant` vs `Keratinocyte_like`.
- Save:
  - Per-cell labels (`*_Stage11_Epithelial_Malignant_vs_KeratinocyteLike_cells.csv`)
  - Cluster-level scores (`*_Stage11_Epithelial_Malignant_vs_KeratinocyteLike_autoSelection.csv`).

### Stage 12 – Epithelial DEGs (Malignant vs Keratinocyte-like)
- Apply labels from Stage 11 and subset epithelial tumor object.
- Run strict differential expression (Wilcoxon, log2FC and detection filters).
- Save:
  - Group cell counts (`*_Stage12_Epithelial_Malignant_vs_KeratinocyteLike_DEGs_group_counts.csv`)
  - DEG table with direction (`*_Stage12_Epithelial_Malignant_vs_KeratinocyteLike_DEGs.csv`).
