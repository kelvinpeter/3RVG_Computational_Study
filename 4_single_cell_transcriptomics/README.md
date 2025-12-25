\# 4. Single-cell Transcriptomic Analysis



This directory contains the complete single-cell RNA-seq analysis pipeline used in this study to generate \*\*Figure 5\*\*.



\## Purpose

To validate a conserved multi-gene inflammatory signature across brain cell populations in Alzheimer’s disease using independent single-nucleus RNA-seq datasets.



\## Software Requirements

\- R (≥ 4.2)

\- Seurat (v4.3)

\- ggplot2

\- dplyr

\- pheatmap



\## Directory Structure

4\_single\_cell\_transcriptomics/

├── data/

│   └── single\_cell/        # Input count matrices (not redistributed)

├── code/

│   └── single\_cell/        # Analysis scripts

├── results/

│   └── figures/            # Output figures

└── README.md



\## Analysis Workflow

1\. Quality control and filtering of single-cell data

2\. Normalization, dimensionality reduction, clustering, and UMAP

3\. Gene-level visualization using feature, violin, and bubble plots



\## Prioritized Genes

TNF, IL1B, CLEC7A, CSNK1E, ARAF, S100A8, APOE, CXCL1, APP, ABL1



\## Data Sources

Publicly available Alzheimer’s disease single-cell atlases, including HuMicA and HumMigA33.



\## Reproducibility

All scripts required to reproduce Figure 5 (panels A–F) are provided in this directory.





