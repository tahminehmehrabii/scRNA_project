# scRNA_project
# Single-cell RNA-seq Example Project

This is a toy example of single-cell RNA-seq analysis using **Seurat** in R.

## Files
- `my_scRNA_counts.csv`: Example count matrix (genes x cells)
- `scRNA_project.R`: Full Seurat analysis pipeline
- `scRNA_example_seurat.rds`: Saved Seurat object after analysis

## Steps
1. Load CSV into Seurat
2. Normalize data
3. Identify variable features
4. Scale data
5. PCA and clustering
6. UMAP visualization
7. Marker gene identification
8. Save Seurat object

## Requirements
- R >= 4.0
- Packages: Seurat, dplyr, patchwork

## How to Run
```R
source("scRNA_project.R")
