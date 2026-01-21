# Load libraries
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(dplyr)

# ---- Inputs ----
expr_file <- "scratch/kallisto_normalized_quants/ortho~maf_focused/norm~tmm_abundance_edger/expr_samples__tmm_abundance_edger__maf_focused.tsv"   # rows = genes, columns = samples
expr_file <- "scratch/kallisto_normalized_interspecies/ortho~maf_focused/norm~tmm_abundance_edger/expr_matrix__tmm_abundance_edger__maf_focused.tsv"   # rows = genes, columns = samples
species_file <- "configs/species_table.tsv"      # first column should be sample IDs
genes_of_interest <- c("CMAF", "NRL", "MAFA", "MAFB")  # example genes
genes_of_interest <- c("CMAF", "NRL", "MAFA", "MAFB", "OTX1", "OTX2", "CRX")  # example genes

# ---- Load data ----
expr <- read.table(expr_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
species_order <- (
    read.delim(species_file)
    |> select(OrganismOrder177, OrganismID, SurrogateID, OrganismShortName, OrganismColor)
)
print(species_order)

dim(expr)
colnames(expr)

samples <- (
    data.frame(sample_id = colnames(expr))
    |> mutate(OrganismID = sub("_[0-9]+$", "", sample_id))
    |> left_join(species_order)
    |> arrange(OrganismOrder177, sample_id)
    |> mutate(short_name = str_replace(sample_id, OrganismID, OrganismShortName))
    |> mutate(ShortNameFactor = factor(OrganismShortName, levels = unique(OrganismShortName), ordered = TRUE))
)
print(samples)

# Ensure samples are ordered as in sample file
expr <- expr[, samples$sample_id]
colnames(expr) <- samples$short_name

# Subset for selected genes
expr_sel <- expr[rownames(expr) %in% genes_of_interest, ]
expr_sel <- expr_sel[genes_of_interest, ]

# Optional: scale per gene (row)
#expr_scaled <- t(scale(t(expr_sel)))
expr_scaled <- log(1+expr_sel)

colors <- (
    read.delim("configs/species_colors.txt")
    |> pull(HexCode2, OrganismColor)
)

# ---- Plot heatmap ----
Heatmap(t(expr_scaled),
        name = "log(1+CPM)",
        #col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        row_split = samples$ShortNameFactor,
        row_title_rot = 0,
        row_title_gp = gpar(fontsize = 4),
        row_gap = unit(0.05, "mm"),
        column_split = c("MAF", "MAF", "MAF", "MAF", "OTX", "OTX", "OTX"),
        column_gap = unit(1.0, "mm"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = TRUE,
        show_row_names = FALSE,
        right_annotation = rowAnnotation(
            df = data.frame(
                Group = samples$OrganismColor,
                Annotation = ifelse(samples$OrganismID == samples$SurrogateID, "Original", "Surrogate")
            ),
            col = list(
                Group = colors,
                Annotation = c("Original" = "green", "Surrogate" = "red")
            )
        )
)
