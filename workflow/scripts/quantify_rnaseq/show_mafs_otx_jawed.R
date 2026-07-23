# Load libraries
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(dplyr)

# rows = genes, columns = samples
expr_file <- snakemake@input[["norm"]]
species_file <- snakemake@input[["orgs"]]
pdf_file <- snakemake@output[["pdf"]]

genes_of_interest <- c(
    "CMAF" = "MAFL",
    "NRL" = "MAFL",
    "MAFA" = "MAFL",
    "MAFB" = "MAFL",
    "MAFF" = "MAFS",
    "MAFG" = "MAFS",
    "MAFK" = "MAFS",
    "OTX1" = "OTX",
    "OTX2" = "OTX",
    "CRX" = "OTX"
)

names(genes_of_interest)
unname(genes_of_interest)

# ---- Load data ----
expr <- read.table(expr_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
species_order <- (
    read.delim(species_file)
    |> select(OrganismOrder181, OrganismID, SurrogateID, OrganismShortName, OrganismColor)
)
print(species_order)

dim(expr)
colnames(expr)

samples <- (
    data.frame(sample_id = colnames(expr))
    |> mutate(OrganismID = sub("_[0-9]+$", "", sample_id))
    |> left_join(species_order)
    |> arrange(OrganismOrder181, sample_id)
    |> mutate(short_name = str_replace(sample_id, OrganismID, OrganismShortName))
    |> mutate(ShortNameFactor = factor(OrganismShortName, levels = unique(OrganismShortName), ordered = TRUE))
)
print(samples)

samples$pos <- 1:nrow(samples)

pos_df <- (
    samples
    |> group_by(OrganismShortName)
    |> summarize(idx = floor(mean(pos)))
    |> ungroup()
)
print(pos_df)


# Ensure samples are ordered as in sample file
expr <- expr[, samples$sample_id]
colnames(expr) <- samples$short_name

# Subset for selected genes
expr_sel <- expr[rownames(expr) %in% names(genes_of_interest), ]
expr_sel <- expr_sel[names(genes_of_interest), ]

# Optional: scale per gene (row)
#expr_scaled <- t(scale(t(expr_sel)))
expr_scaled <- log(1+expr_sel)

colors <- (
    read.delim("configs/species_colors.txt")
    |> pull(HexCode2, OrganismColor)
)

# IMPORTANT: use numeric indices in 'at' below (not rownames)
la <- rowAnnotation(
  group_label = anno_mark(
    at = pos_df$idx,            # <--- numeric indices (1..nrow(mat)), NOT row names
    labels = pos_df$OrganismShortName,             # shown labels
    labels_rot = 0,
    padding = unit(0.02, "mm"),
    which = "row",
    labels_gp = gpar(fontsize = 7)
  ),
  width = unit(3, "cm")
)

# ---- Plot heatmap ----
hm <- Heatmap(t(expr_scaled),
        name = "log(1+CPM)",
        col = colorRamp2(c(0, 5, 10), c("blue", "white", "red")),
        row_split = samples$ShortNameFactor,
        row_title_rot = 0,
        row_title_gp = gpar(fontsize = 4),
        row_title = NULL,
        row_gap = unit(0.05, "mm"),
        column_split = unname(genes_of_interest),
        column_gap = unit(1.0, "mm"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = TRUE,
        show_row_names = FALSE,
        right_annotation = la,
        left_annotation = rowAnnotation(
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
pdf(pdf_file, height = 7.5, width = 5.5)
draw(hm, merge_legends = TRUE)
dev.off()
