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
    "MAFL" = "MAFL",
    "MAFS" = "MAFS",
    "OTX" = "OTX"
)

names(genes_of_interest)
unname(genes_of_interest)

# ---- Load data ----
expr <- (
    read.table(expr_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
)



organisms <- (
    read.delim(species_file)
    |> select(OrganismOrder181, OrganismID, SurrogateID, OrganismShortName, OrganismColor)
    |> arrange(OrganismOrder181)
    |> mutate(OrganismColorFactor = factor(OrganismColor, levels = unique(OrganismColor), ordered = TRUE))
)
print(organisms)


# Subset for selected genes
expr_sel <- expr[rownames(expr) %in% names(genes_of_interest), ]
expr_sel <- expr_sel[names(genes_of_interest), ]

# Optional: scale per gene (row)
#expr_scaled <- t(scale(t(expr_sel)))
expr_scaled <- log(1+expr_sel)

expr_scaled <- (
    as.data.frame(expr_scaled)
    |> rownames_to_column("gene")
    |> gather(key = "sample_id", value = "expr", -gene)
    |> mutate(OrganismID = sub("_[0-9]+$", "", sample_id))
    |> group_by(gene, OrganismID)
    |> summarize(expr = median(expr, na.rm = TRUE))
    |> spread(OrganismID, expr)
    |> column_to_rownames("gene")
    |> as.matrix()
)
print(expr_scaled)

organisms <- filter(organisms, OrganismID %in% colnames(expr_scaled))

colors <- (
    read.delim("configs/species_colors.txt")
    |> pull(HexCode2, OrganismColor)
)


## Ensure rows and columns are ordered properly
expr_scaled <- expr_scaled[names(genes_of_interest), organisms$OrganismID]
colnames(expr_scaled) <- organisms$OrganismShortName

# ---- Plot heatmap ----
hm <- Heatmap(t(expr_scaled),
        name = "log(1+CPM)",
        col = colorRamp2(c(0, 5, 10), c("blue", "white", "red")),
        row_split = organisms$OrganismColorFactor,
        row_title_rot = 0,
        row_title_gp = gpar(fontsize = 4),
        row_title = NULL,
        row_gap = unit(1.0, "mm"),
        column_split = unname(genes_of_interest),
        column_gap = unit(1.0, "mm"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = TRUE,
        show_row_names = TRUE,
        left_annotation = rowAnnotation(
            df = data.frame(
                Group = organisms$OrganismColor,
                Annotation = ifelse(organisms$OrganismID == organisms$SurrogateID, "Original", "Surrogate")
            ),
            col = list(
                Group = colors,
                Annotation = c("Original" = "green", "Surrogate" = "red")
            )
        )
)
pdf(pdf_file, height = 7.5, width = 15.5)
draw(hm, merge_legends = TRUE)
dev.off()
