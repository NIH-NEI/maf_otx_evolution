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
              help = "Optional: Input tab-delimited protein data file with columns [identifier  species/group gene]. No header.
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

# Input data
proteindata_table <- read.table(file = parser$protein_data, col.names = c("identifier", "species_group", "gene_group"), sep = "\t")
similarity_matrix <- read.table(file = parser$matrix, skip = 1, sep = "\t", row.names = 1)
colnames(similarity_matrix) <- rownames(similarity_matrix)

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
  color_table <- hue_pal()(length(unique(proteindata_table$species_group)))
  names(color_table) <- unique(proteindata_table$species_group)
}

# Run CA, export results
res.ca <- CA(similarity_matrix, graph = FALSE)
row_results <- res.ca$row$coord %>%
  as.data.frame() %>%
  rownames_to_column(var = "identifier") %>%
  full_join(., proteindata_table, by = "identifier")

# Generate Correspondence Analysis plot
CA_plot <-
  row_results %>%
  ggplot(aes(x = `Dim 1`, y = `Dim 2`, shape = gene_group, color = species_group)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 1) + geom_vline(xintercept = 0, color = "black", linewidth = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = color_table) +
    ggtitle(plot_title[1]) +
    xlab("Dim 1") + ylab("Dim 2") +
    labs(color = "Species", shape = "Gene") +
    theme(plot.title = element_text(hjust = 0.5))

# Generate Heatmap

heatmap_data <- data.matrix(similarity_matrix)
heatmap_cols <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)
heatmap_anno <- data.frame(Species = proteindata_table$species_group, row.names = proteindata_table$identifier)
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
dev.off()

