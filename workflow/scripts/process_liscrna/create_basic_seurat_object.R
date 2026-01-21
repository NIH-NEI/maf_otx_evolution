library(Seurat)

# Folder that contains: barcodes.tsv.gz, features.tsv.gz, matrix.mtx (or matrix.mtx.gz)
data_dir <- snakemake@params[["data_dir"]]
seurat_file <- snakemake@output[["seu"]]
pdf_file <- snakemake@output[["pdf"]]

print(data_dir)
print(seurat_file)
print(pdf_file)


#data_dir <- "imports/liscrna/pig/pig_li_1/"
#data_dir <- "imports/liscrna/crabMacaque/crabMacaque_li_1/"

pdf(pdf_file, height = 10, width = 10)

# Read 10x Genomics-formatted data
counts <- Read10X(data.dir = data_dir, gene.column = 1)

# Standard Seurat workflow (mostly defaults; UMAP/Neighbors need dims specified)
seu <- CreateSeuratObject(counts = counts, project = "sample1")

cat("Total cells :", ncol(seu), "\n")
cat("Empty cells (0 genes):", sum(seu$nFeature_RNA == 0), "\n")
cat("Very low quality (<200 genes):", sum(seu$nFeature_RNA < 200), "\n")
cat("Very low quality (<100 genes):", sum(seu$nFeature_RNA < 100), "\n")
cat("Very low quality (<50 genes):", sum(seu$nFeature_RNA < 50), "\n")

VlnPlot(seu,
        features = c("nCount_RNA", "nFeature_RNA"),
        log = TRUE,
        pt.size = 0.0)


#seu <- NormalizeData(seu)
#seu <- FindVariableFeatures(seu)
#seu <- ScaleData(seu)
#seu <- RunPCA(seu, npcs=50)
#
#ElbowPlot(seu, ndims = 50)
#
## Common default-style choice is to use the first 30 PCs
#seu <- RunUMAP(seu, dims = 1:30)
#DimPlot(seu, reduction = "umap")
#FeaturePlot(seu, reduction = "umap", features = c("RHO", "OPN1LW", "VSX2", "PAX6"))
#
#seu <- FindNeighbors(seu, dims = 1:20)
#seu <- FindClusters(seu)
#
#p1 <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters")
#p2 <- FeaturePlot(seu, reduction = "umap", features = c("RHO", "OPN1LW", "VSX2", "PAX6"))
#
#p1 + p2

dev.off()

#saveRDS(seu, "pig_seu.RDS")
#saveRDS(seu, "crabMacaque_seu.RDS")
saveRDS(seu, seurat_file)
