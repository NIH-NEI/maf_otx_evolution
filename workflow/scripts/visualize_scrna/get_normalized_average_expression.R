suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(purrr)
    library(stringr)
    library(tibble)
})

seurat_file <- snakemake@input[["seu"]]
tsv_file <- snakemake@output[[1]]
species <- snakemake@wildcards[["org"]]

print(seurat_file)
print(tsv_file)
print(species)

UpperCase_genes = function(object, integration = FALSE){

  rownames(object@assays$RNA@counts)<- toupper(rownames(object@assays$RNA@counts))
  rownames(object@assays$RNA@data)<- toupper(rownames(object@assays$RNA@data))
  rownames(object@assays$RNA@scale.data)<- toupper(rownames(object@assays$RNA@scale.data))
  object@assays$RNA@var.features <- toupper(object@assays$RNA@var.features)
  #rownames(object@assays$RNA@meta.features)<- toupper(rownames(object@assays$RNA@meta.features))

  if (integration){
    rownames(object@assays$integrated@counts)<- toupper(rownames(object@assays$integrated@counts))
    rownames(object@assays$integrated@data)<- toupper(rownames(object@assays$integrated@data))
    rownames(object@assays$integrated@scale.data)<- toupper(rownames(object@assays$integrated@scale.data))
    object@assays$integrated@var.features <- toupper(object@assays$integrated@var.features)
    #rownames(object@assays$integration@meta.features)<- toupper(rownames(object@assays$integration@meta.features))
  }

  return(object)
}

features <- c("NRL", "MAF", "MAFA", "MAFB", "MAFF", "MAFG", "MAFK", "CRX", "OTX1", "OTX2")
classes <- c("RGC", "BP", "GabaAC", "GlyAC", "HC", "Cone", "Rod", "MG")

obj <- readRDS(seurat_file)
if (species %in% c("deerMouse")) {
obj <- UpdateSeuratObject(obj)
}
obj <- UpperCase_genes(obj)
DefaultAssay(obj) <- "RNA"

Idents(obj) <- "cell_class"
obj <- subset(obj, idents = classes)

obj@meta.data$cell_class0 = factor(obj@meta.data$cell_class, levels = classes)
levels(obj@meta.data$cell_class0) = c("RGC","BC","AC","AC","HC","Cone","Rod", "MG")
Idents(obj) <- "cell_class0"

# Calculate average expression for all genes by class
features_in <- features[features %in% rownames(obj)]

avg_exp <- AverageExpression(obj, features = features_in, assays = "RNA", slot = "counts")$RNA
#colnames(avg_exp) <- paste0(species, "_", colnames(avg_exp))

# Normalize
log_exp <- log1p(avg_exp)
class_norm <- avg_exp / apply(avg_exp, 1, max)

print(log_exp)
print(class_norm)

df = (
    as.data.frame(class_norm)
    |> rownames_to_column("gene")
    |> mutate(OrganismShortName = species)
)
print(df)
readr::write_tsv(df, tsv_file)
