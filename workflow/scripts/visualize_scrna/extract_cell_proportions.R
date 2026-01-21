suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(purrr)
    library(stringr)
})

seurat_file <- snakemake@input[["seu"]]
tsv_file <- snakemake@output[[1]]
organism <- snakemake@wildcards[["org"]]
sample_id <- paste0(organism, "_hahn_1")

#organism <- "human"
#sample_id <- "human_hahn_1"
#seurat_file <- "/data/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Human_initial.rds"
#tsv_file <- "haha.tsv"

print(seurat_file)
print(tsv_file)
print(organism)

seurat_obj <- readRDS(seurat_file)

df <- (
    seurat_obj@meta.data
    |> mutate(count = 1)
    |> group_by(cell_class)
    |> summarize(count = sum(count))
    |> ungroup()
    |> mutate(OrganismShortName = organism)
    |> mutate(SampleID = sample_id)
    |> select(OrganismShortName, SampleID, cell_class, count)
)
print(df)
readr::write_tsv(df, tsv_file)
