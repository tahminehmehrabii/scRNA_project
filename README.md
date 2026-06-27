# LSCC Single-Cell RNA-seq Analysis and Refined inferCNV Pipeline

This repository contains an R-based single-cell RNA-seq workflow for the analysis of laryngeal squamous cell carcinoma (LSCC) samples from **GSE206332**. The pipeline integrates quality control, clustering, cell-type annotation, CNV-based malignant-cell identification, malignant subclustering, refined inferCNV analysis, and High-CNV marker discovery.

## **Dataset**

* **Dataset:** GSE206332
* **Samples:** GSM6251294, GSM6251297, and GSM6251300
* **Disease:** Laryngeal squamous cell carcinoma (LSCC)
* **Data type:** 10x Genomics single-cell RNA-seq data

## **Required R Packages**

The pipeline requires the following R packages:

```r
Seurat
SeuratObject
Matrix
data.table
dplyr
ggplot2
harmony
infercnv
patchwork
pheatmap
grid
png
ggrepel
```

## **Analysis Workflow**



### **1. Environment Setup**

The script clears the R environment, defines global options, sets a random seed, loads required packages, and creates output directories for figures, RDS objects, and inferCNV results.

### **2. Input Data Preparation**

Raw 10x Genomics files, including the matrix, barcode, and feature files, are loaded for each LSCC sample. Gene-order information based on the hg38 reference genome is used for inferCNV analysis.

### **3. Quality Control**

Cells are filtered according to the following criteria:

* Retain cells with at least **200 detected genes**.
* Remove cells with mitochondrial transcript percentages of **20% or higher**.
* Generate a QC summary table showing cell counts before and after filtering.

### **4. Data Integration and Dimensionality Reduction**

Filtered samples are merged and processed using the following steps:

* Log-normalization of gene-expression data.
* Identification of 2,000 highly variable genes.
* Scaling of gene-expression values.
* Principal component analysis (PCA).
* Harmony batch correction based on sample identity.
* Graph-based clustering and t-SNE visualization.

### **5. Global Cluster Marker Detection**

Differential-expression analysis is performed across global cell clusters using `FindAllMarkers`. A dedicated differential-expression assay is generated to avoid layer-related issues in Seurat v5.

Marker tables are exported for:

* All detected cluster markers.
* Statistically significant global cluster markers.

### **6. Cell-Type Annotation**

Major cell populations are annotated using canonical marker genes and module-score analysis.

The annotated cell types include:

* B cells
* Endothelial cells
* Epithelial cells
* Fibroblasts
* Myeloid cells
* NK cells
* T cells

t-SNE plots and marker DotPlots are generated to visualize cell-type assignments.

### **7. Initial inferCNV Analysis**

Initial inferCNV analysis is performed separately for each sample.

* **Epithelial cells** are evaluated as potential malignant populations.
* **Myeloid cells** are used as the reference population.
* CNV scores are calculated from inferCNV expression profiles.
* Epithelial cells with CNV scores above the 95th percentile of reference myeloid cells are classified as malignant-like cells.

### **8. Malignant Cell Re-Clustering**

Malignant epithelial cells are extracted and re-clustered to investigate intratumoral heterogeneity.

The workflow includes:

* Normalization and variable-feature selection.
* PCA and graph-based clustering.
* t-SNE visualization of malignant epithelial subclusters.
* Calculation of initial CNV-score summaries for each malignant subcluster.

### **9. Low-CNV Reference Cluster Selection**

Weak or low-CNV malignant subclusters are identified using the initial CNV-score distribution.

The script supports two approaches:

* **Automatic selection:** Clusters in the lowest CNV-score quartile are selected.
* **Manual selection:** Users can specify low-CNV clusters after inspecting the initial CNV violin plot.

Selected low-CNV clusters are used as internal reference populations for refined inferCNV analysis.

### **10. Refined inferCNV Analysis**

Refined inferCNV is performed using the selected low-CNV malignant subclusters as reference groups.

This step includes:

* Generation of refined CNV scores.
* Classification of malignant subclusters as Low-CNV or High-CNV.
* Visualization of refined CNV-score distributions.
* Export of inferCNV heatmaps and CNV summary tables.

### **11. High-CNV Malignant Cell Selection**

Low-CNV malignant subclusters are removed after refined inferCNV analysis.

The remaining cells represent High-CNV malignant epithelial populations and are used for downstream marker identification.

### **12. High-CNV Marker Gene Identification**

Cluster-specific marker genes are identified among the remaining High-CNV malignant subclusters.

Markers are retained using the following criteria:

* `avg_log2FC ≥ 1`
* `min.pct ≥ 0.25`
* Adjusted `P-value < 0.05`

The pipeline exports:

* Full High-CNV marker tables.
* Filtered High-CNV marker tables.
* Final unique marker-gene lists in CSV and TXT formats.

### **13. Figure Generation**

The workflow generates publication-ready figures, including:

* Global t-SNE clusters.
* Annotated cell-type t-SNE plot.
* Canonical marker DotPlot.
* Malignant epithelial subcluster t-SNE plot.
* Initial CNV-score violin plot.
* Refined inferCNV heatmap.
* Refined CNV-score violin plot.
* Combined multi-panel manuscript figure.

All figures are exported at high resolution.

## **Figure Reconstruction Scripts**

Additional scripts are included for reconstructing publication-ready figures without rerunning the complete single-cell analysis.

These scripts:

* Load existing Seurat and malignant-cell RDS objects.
* Reuse existing inferCNV output images.
* Do not rerun QC, normalization, Harmony, clustering, inferCNV, or marker detection.
* Generate standardized high-resolution PNG, TIFF, and PDF figures.
* Improve label placement in t-SNE plots using non-overlapping labels.
* Create enlarged and cropped inferCNV panels for manuscript preparation.

## **Main Output Files**

### **Figures**

```text
Figure_SC_01A_tSNE_global_clusters.png
Figure_SC_01B_tSNE_final_celltypes_with_malignant.png
Figure_SC_01C_DotPlot_celltype_markers.png
Figure_SC_01D_tSNE_malignant_subclusters.png
Figure_SC_01E_Initial_CNV_score_violin_for_reference_selection.png
Figure_SC_01F_inferCNV.png
Figure_SC_01G_Refined_CNV_score_violin_malignant_subclusters.png
```

### **Final High-CNV Marker Files**

```text
Final_High_CNV_Malignant_Genes.csv
Final_High_CNV_Malignant_Gene_Names.csv
Final_High_CNV_Malignant_Gene_Names.txt
High_CNV_Malignant_Cluster_FindAllMarkers_after_lowCNV_removed_logFC1_all.csv
```

### **RDS Objects**

```text
FINAL_seurat_object_with_malignant_annotation.rds
MALIGNANT_all_object_with_refined_CNV.rds
HIGH_CNV_MALIGNANT_object_after_lowCNV_removal.rds
HIGH_CNV_MALIGNANT_object_DE_fixed.rds
```

## **Notes**

* Update all directory paths before running the script.
* Ensure that the hg38 gene-order file is available before starting inferCNV analysis.
* The figure-reconstruction scripts require previously generated RDS objects and inferCNV output files.
* High-CNV marker genes are identified only after low-CNV malignant subclusters are removed.
