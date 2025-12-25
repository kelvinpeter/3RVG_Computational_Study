library(Seurat)
library(pheatmap)

seu <- readRDS("../../results/seurat_processed.rds")

genes <- c("TNF","IL1B","CLEC7A","CSNK1E","ARAF",
           "S100A8","APOE","CXCL1","APP","ABL1")

VlnPlot(seu, features = genes, stack = TRUE, flip = TRUE)

avg.exp <- AverageExpression(seu, features = genes)$RNA
pheatmap(avg.exp, scale = "row")
