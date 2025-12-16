# DEG Analysis for GSE16161 - Psoriasis Dataset
# Author: [Saleem Iqbal]
# Date: [Dec-17-2025]

# Load required libraries
library(GEOquery)
library(limma)
library(umap)
library(maptools)

# Load series and platform data from GEO
gset <- getGEO("GSE16161", GSEMatrix = TRUE, AnnotGPL = TRUE)

if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Make proper column names
fvarLabels(gset) <- make.names(fvarLabels(gset))

# Group membership for samples (Lesional vs Normal)
gsms <- "000000000111111111"
sml <- strsplit(gsms, split = "")[[1]]

# Log2 transformation if needed
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

if (LogC) {
  ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex)
}

# Assign groups
gs <- factor(sml)
groups <- make.names(c("Lesional", "Normal"))
levels(gs) <- groups
gset$group <- gs

# Design matrix
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

# Remove missing values
gset <- gset[complete.cases(exprs(gset)), ]

# Fit linear model
fit <- lmFit(gset, design)

# Contrasts
cts <- c(paste(groups[1], "-", groups[2], sep = ""))
cont.matrix <- makeContrasts(contrasts = cts, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)

# Empirical Bayes
fit2 <- eBayes(fit2, 0.01)

# Top table
tT <- topTable(fit2, adjust = "fdr", sort.by = "B", number = 250)
tT <- subset(tT, select = c("ID", "adj.P.Val", "P.Value", "t", "B", "logFC", 
                            "GB_ACC", "SPOT_ID", "Gene.Symbol", "Gene.symbol", 
                            "Gene.title"))

# Save results
write.csv(tT, "GSE16161_DEG_results.csv", row.names = FALSE)

# Create output directory for plots
dir.create("Plots", showWarnings = FALSE)

# Visualization
# 1. Histogram of adjusted p-values
pdf("Plots/Padj_distribution.pdf")
hist(tT$adj.P.Val, col = "grey", border = "white", 
     xlab = "P-adj", ylab = "Number of genes", 
     main = "P-adj value distribution")
dev.off()

# 2. Volcano plot
pdf("Plots/Volcano_plot.pdf")
volcanoplot(fit2, coef = 1, main = "Volcano Plot: Lesional vs Normal", 
            pch = 20, highlight = length(which(decideTests(fit2)[, 1] != 0)))
dev.off()

# 3. Boxplot
pdf("Plots/Expression_boxplot.pdf")
ord <- order(gs)
palette(c("#1B9E77", "#7570B3"))
boxplot(exprs(gset)[, ord], boxwex = 0.6, notch = TRUE, 
        main = "GSE16161 Expression Distribution", 
        outline = FALSE, las = 2, col = gs[ord])
legend("topleft", groups, fill = palette(), bty = "n")
dev.off()

# 4. UMAP
ex_clean <- na.omit(exprs(gset))
ex_clean <- ex_clean[!duplicated(ex_clean), ]
ump <- umap(t(ex_clean), n_neighbors = 8, random_state = 123)

pdf("Plots/UMAP_plot.pdf", width = 8, height = 6)
par(mar = c(3, 3, 2, 6), xpd = TRUE)
plot(ump$layout, main = "UMAP Plot", xlab = "", ylab = "", 
     col = as.numeric(gs), pch = 20, cex = 1.5)
legend("topright", inset = c(-0.15, 0), legend = levels(gs), 
       pch = 20, col = 1:nlevels(gs), title = "Group", pt.cex = 1.5)
dev.off()

# Session info
sink("session_info.txt")
sessionInfo()
sink()