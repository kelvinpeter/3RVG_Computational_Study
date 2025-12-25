library(Seurat)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)

seu <- RunPCA(seu, npcs = 30)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.6)

seu <- RunUMAP(seu, dims = 1:20)
saveRDS(seu, "../../results/seurat_processed.rds")
