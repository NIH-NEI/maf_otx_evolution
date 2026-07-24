suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(purrr)
    library(stringr)
    library(patchwork)
    library(tibble)
    library(ComplexHeatmap)
    library(circlize)
})

species_file <- "configs/species181_table.tsv"      # first column should be sample IDs
genes_of_interest <- c("MAF", "NRL", "MAFA", "MAFB")  # example genes
genes_of_interest <- c("MAF", "NRL", "MAFA", "MAFB", "OTX1", "OTX2", "CRX")  # example genes


#tsv_file <- "scratch/scrna_dotplots/allorgs_avg_expr_normalized.tsv"

tsv_file <- snakemake@input[[1]]
pdf_file <- snakemake@output[[1]]
print(tsv_file)
print(pdf_file)

expr <- (
    read.delim(tsv_file)
    |> select(gene, OrganismShortName, Rod)
    |> filter(gene %in% genes_of_interest)
    |> filter(OrganismShortName != "seaLamprey")
    |> spread(key = OrganismShortName, value = Rod)
    |> column_to_rownames("gene")
)
print(expr)


species_order <- (
    read.delim(species_file)
    |> select(OrganismOrder181, OrganismID, SurrogateID, OrganismShortName, OrganismColor)
)
print(species_order)

dim(expr)
colnames(expr)

samples <- (
    data.frame(OrganismShortName = colnames(expr))
    |> left_join(species_order)
    |> arrange(OrganismOrder181)
    |> mutate(ShortNameFactor = factor(OrganismShortName, levels = unique(OrganismShortName), ordered = TRUE))
)

# Ensure samples are ordered as in sample file
expr <- expr[, samples$OrganismShortName]

# Subset for selected genes
expr_sel <- expr[rownames(expr) %in% genes_of_interest, ]
expr_sel <- expr_sel[genes_of_interest, ]

# Optional: scale per gene (row)
#expr_scaled <- t(scale(t(expr_sel)))
#expr_scaled <- log(1+expr_sel)

colors <- (
    read.delim("configs/species_colors.txt")
    |> pull(HexCode2, OrganismColor)
)

# ---- Plot heatmap ----
hm <- Heatmap(expr_sel,
        name = "NormAvgExpr",
        col = colorRamp2(c(0, 1), c("white", "#447241")),
        column_split = samples$ShortNameFactor,
        column_title_rot = 90,
        column_title_gp = gpar(fontsize = 9),
        column_gap = unit(0.05, "mm"),
        row_split = c("MAF", "MAF", "MAF", "MAF", "OTX", "OTX", "OTX"),
        row_gap = unit(1.0, "mm"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        top_annotation = columnAnnotation(
            df = data.frame(
                Group = samples$OrganismColor
                #, Annotation = ifelse(samples$OrganismID == samples$SurrogateID, "Original", "Surrogate")
            ),
            col = list(
                Group = colors
                #, Annotation = c("Original" = "green", "Surrogate" = "red")
            )
        )
)

pdf(pdf_file, height = 6, width = 8)
draw(hm, merge_legends = TRUE)
dev.off()
