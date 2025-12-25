library(Seurat)
library(ggplot2)

seu <- readRDS("../../results/seurat_processed.rds")

genes <- c("TNF","IL1B","CLEC7A","CSNK1E","ARAF",
           "S100A8","APOE","CXCL1","APP","ABL1")

pdf("../../results/figures/Figure5_UMAP_FeaturePlots.pdf")
for (g in genes) {
  print(FeaturePlot(seu, features = g))
}
dev.off()
