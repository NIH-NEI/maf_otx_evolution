### Packages

suppressPackageStartupMessages({
  library(FactoMineR)
  library(factoextra)
  library(ggplot2)
  library(ggforce)
  library(optparse)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
  library(patchwork)
})

### Parser

arg_list <- list(
  make_option(opt_str = c("-m", "--matrix"), type = "character", default = "ordered",
              help = "Optional: Input similarity matrix from DALI. 
              Default: ordered"),
  make_option(opt_str = c("-p", "--protein_data"), type = "character", default = "proteindata.txt",
              help = "Optional: Input tab-delimited protein data file with columns [ProteinID  species/group gene]. No header.
              Default: proteindata.txt"),
  make_option(opt_str = c("-t", "--plot_title"), type = "character", default = NA,
              help = "Optional: Information to be appended to plot title.
              Default: Correspondence Analysis plot | Heatmap"),
  make_option(opt_str = c("-c", "--colors"), type = "character", default = NA,
              help = "Optional: Table of hexcode colors to be used with columns [species/group  color]. No header.
              Default: default ggplot colors")
)

parser <- parse_args(OptionParser(usage = "Rscript correspondence_analysis_facto.R [-i [ordered] | -p [proteindata.txt] | -t [title] | -c [color file]", option_list = arg_list, description = "R script to produce correspondence analysis plot and heatmap from DALI results. Output PDF file containing both plots is created in working directory."))

### Test cases

# parser <- parse_args(OptionParser(usage = "Rscript correspondence_analysis_facto.R [-i [ordered] | -p [proteindata.txt] | -t title]", option_list = arg_list, description = "R script to produce correspondence analysis plot and heatmap from DALI results. Output PDF file containing both plots is created in working directory."), args = c("--colors=color_table.txt", "--plot_title=title"))

### Script

org_colors_orig <- (
    read.delim("configs/species_colors.txt")
)
shortColor <- pull(org_colors_orig, OrganismShortColor, OrganismColor)

org_colors <- (
    org_colors_orig
    |> mutate(OrganismColor = factor(OrganismShortColor, levels = unique(OrganismShortColor)))
    |> pull(HexCode2, OrganismColor)
)
print(org_colors)

query2pid <- (
    read.delim("imports/FromLeo/query.list", header = FALSE, col.names = c("QueryID", "ProteinID"))
)
head(query2pid)

# Input data
proteindata_table <- (
    read.delim(file = parser$protein_data, sep = "\t")
    |> inner_join(query2pid)
    |> filter(GeneFamily %in% c("LargeMAF", "Outgroup"))
    |> filter(GeneGroup != "Ignore")
    #|> filter(GeneGroup != "CMAF")
    #|> filter(ProteinSource != "Denovo")
    |> select(ProteinID, OrganismColor, GeneGroup, GeneFamily, QueryID)
    |> mutate(OrganismColor = shortColor[OrganismColor])
)
head(proteindata_table)

similarity_matrix <- read.table(file = parser$matrix, skip = 1, sep = "\t", row.names = 1)
colnames(similarity_matrix) <- rownames(similarity_matrix)

similarity_matrix <- similarity_matrix[proteindata_table$QueryID, proteindata_table$QueryID]

similarity_matrix[1:10,1:10]

rownames(similarity_matrix) <- proteindata_table$ProteinID #query2pid[rownames(similarity_matrix)]
colnames(similarity_matrix) <- proteindata_table$ProteinID #query2pid[colnames(similarity_matrix)]


# Default title and colors
if(!is.na(parser$plot_title)) {
  plot_title <- c(paste("DALI Correspondence Analysis,", parser$plot_title), paste("DALI Similarity Heatmap,", parser$plot_title))
} else {
  plot_title <- c("DALI Correspondence Analysis", "DALI Similarity Heatmap")
}

if(!is.na(parser$colors)) {
  color_table <- read.table(file = parser$colors, sep = "\t", comment.char = "")$V2
  names(color_table) <- read.table(file = parser$colors, sep = "\t", comment.char = "")$V1
} else {
  color_table <- hue_pal()(length(unique(proteindata_table$OrganismColor)))
  names(color_table) <- unique(proteindata_table$OrganismColor)
}

# Run CA, export results
res.ca <- CA(similarity_matrix, graph = FALSE)
row_results <- res.ca$row$coord %>%
  as.data.frame() %>%
  rownames_to_column(var = "ProteinID") %>%
  full_join(., proteindata_table, by = "ProteinID")

# Generate Correspondence Analysis plot
CA_plot <-
  row_results %>%
  ggplot(aes(x = `Dim 1`, y = `Dim 2`, shape = GeneGroup, color = OrganismColor)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 1) + geom_vline(xintercept = 0, color = "black", linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = org_colors) +
    ggtitle(plot_title[1]) +
    xlab("Dim 1") + ylab("Dim 2") +
    ylim(c(-0.5,1)) +
    xlim(c(-0.5,1)) +
    labs(color = "Species", shape = "Gene") +
    theme(plot.title = element_text(hjust = 0.5))

# Generate Heatmap

heatmap_data <- data.matrix(similarity_matrix)
heatmap_cols <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)
heatmap_anno <- data.frame(Species = proteindata_table$OrganismColor, row.names = proteindata_table$ProteinID)
heatmap_plot <- pheatmap(heatmap_data, 
                         scale = "none",
                         cluster_cols = TRUE, cluster_rows = TRUE,
                         color = heatmap_cols,
                         # height = 5, width = 5,
                         legend = FALSE, legend_breaks = c(4, 10, 16),
                         show_rownames = F, show_colnames = F, angle_col = 90, main = plot_title[2],
                         annotation_col = heatmap_anno, annotation_row = heatmap_anno,
                         annotation_names_row = F, annotation_names_col = F,
                         annotation_colors = list(Species = color_table), silent = TRUE)

# ComplexHeatmap version of your pheatmap call
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)

# data + colors (same palette logic as your code)
heatmap_data <- as.matrix(similarity_matrix)
heatmap_cols <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)
col_fun <- colorRamp2(seq(min(heatmap_data, na.rm = TRUE),
                          max(heatmap_data, na.rm = TRUE),
                          length.out = length(heatmap_cols)),
                      heatmap_cols)

# annotation data (make sure rownames match matrix dimnames)
heatmap_anno <- data.frame(
  Species = proteindata_table$OrganismColor,
  row.names = proteindata_table$ProteinID
)

# keep only IDs present in the matrix (avoids mismatches)
common_ids <- intersect(rownames(heatmap_data), rownames(heatmap_anno))
heatmap_data <- heatmap_data[common_ids, common_ids, drop = FALSE]
heatmap_anno  <- heatmap_anno[common_ids, , drop = FALSE]

# annotations (column and row)
ha_col <- HeatmapAnnotation(
  df = heatmap_anno,
  col = list(Species = color_table),
  show_annotation_name = FALSE
)

ha_row <- rowAnnotation(
  df = heatmap_anno,
  col = list(Species = color_table),
  show_annotation_name = FALSE
)

# heatmap (cluster rows/cols TRUE, hide row/col names, no legend)
ht <- Heatmap(
  heatmap_data,
  name = NULL,                     # no legend title (and legend turned off below)
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = ha_col,
  left_annotation = ha_row,
  column_title = plot_title[2],
  heatmap_legend_param = list(at = c(4, 10, 16))  # like your legend_breaks (legend hidden anyway)
)

# Output files

if(!is.na(parser$plot_title)) {
  file_title <- c(paste(parser$plot_title, "-DALI_plots.pdf", sep = ""))
} else {
  file_title <- c("DALI_plots.pdf")
}

pdf(file = file_title, width = 8, height = 6)
CA_plot
grid::grid.newpage()
heatmap_plot
# draw (legend = FALSE like your pheatmap call; also hides annotation legends)
draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     show_heatmap_legend = FALSE,
     show_annotation_legend = FALSE)
dev.off()

