library(Seurat)
library(dplyr)

sc <- Read10X(data.dir = "../../data/single_cell/")
seu <- CreateSeuratObject(counts = sc, min.cells = 3, min.features = 200)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

seu <- subset(seu,
              subset = nFeature_RNA > 500 &
                       nFeature_RNA < 6000 &
                       percent.mt < 10)