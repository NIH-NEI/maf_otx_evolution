suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(purrr)
    library(stringr)
})

tsv_file <- snakemake@input[[1]]
pdf_file <- snakemake@output[[1]]

print(tsv_file)

df <- (
    read.delim(tsv_file)
    # optional: keep a stable ordering
    |> mutate(
        object = factor(object, levels = unique(object)),
        # row label with hierarchy
        row_label = paste(object, celltype, sep = " | ")
    )
    |> filter(celltype == "Rod")
)
print(df)

genes <- c("NRL", "MAF", "MAFA", "MAFB", "MAFF", "MAFG", "MAFK", "CRX", "OTX1", "OTX2")

# -----------------------------
# Plot: single ggplot
# -----------------------------
plot_mega_dotplot <- function(
        df,
        genes,
        dot_max_size = 6,
        x_text_angle = 60
) {
    # Keep columns in the supplied gene order
    df <- df %>%
        mutate(gene = factor(gene, levels = genes)) %>%
        # y order: object then celltype (grouped)
        arrange(object, celltype) %>%
        mutate(row_label = factor(row_label, levels = unique(row_label)))

    ggplot(df, aes(x = gene, y = row_label)) +
        geom_point(aes(size = pct_exp, color = avg_exp), alpha = 0.9) +
        scale_size(range = c(0, dot_max_size)) +
        theme_classic(base_size = 12) +
        theme(
            axis.title = element_blank(),
            axis.text.x = element_text(angle = x_text_angle, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 9),
            panel.grid.major = element_line(linewidth = 0.2, color = "grey92")
        )
}

p <- plot_mega_dotplot(df, genes)

pdf(pdf_file, width = 8, height = 6)
p
dev.off()
