suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(purrr)
    library(stringr)
})

tsv_file <- snakemake@input[[1]]
pdf_file <- snakemake@output[[1]]

print(tsv_file)

cell_types <- list(
    Rod = c("Rod"),
    Cone = c("Cone"),
    BP = c("BP"),
    HC = c("HC"),
    AC = c("GabaAC", "GlyAC"),
    RGC = c("RGC"),
    MG = c("MG"),
    Other = c("MicroG", "Other")
)
#print(cell_types)

cell_types_df <- stack(cell_types)
colnames(cell_types_df) <- c("cell_type", "cell_class")
print(cell_types_df)

type2class = pull(cell_types_df, cell_class, cell_type)
print(type2class)



df <- (
    read.delim(tsv_file)
    # optional: keep a stable ordering
    |> mutate(
        SampleID = factor(SampleID, levels = unique(SampleID)),
    )
    |> rename(cell_type := cell_class)
    |> mutate(cell_class = type2class[cell_type])
    |> group_by(SampleID, cell_class)
    |> summarize(count = sum(count))
)
print(df)

total_cells <- df |> group_by(SampleID) |> summarize(total_cells = sum(count)) |> ungroup()
print(total_cells)

df <- (
    df
    |> left_join(total_cells)
    |> mutate(celltype_proportion = count / total_cells)
)
print(df)

photoreceptors <- (
    df |> filter(cell_class %in% c("Rod", "Cone"))
)
print(photoreceptors)
total_photoreceptors <- photoreceptors |> group_by(SampleID) |> summarize(total_photoreceptors = sum(count)) |> ungroup()
print(total_photoreceptors)

photoreceptors <- (
    photoreceptors
    |> left_join(total_photoreceptors)
    |> mutate(pr_proportion = count / total_photoreceptors)
)
print(photoreceptors)

# -----------------------------
# Plot: single ggplot
# -----------------------------
plot_proportions <- function(
        df,
        col,
        dot_max_size = 6,
        x_text_angle = 60
) {
    ## Keep columns in the supplied gene order
    #df <- df %>%
    #    mutate(gene = factor(gene, levels = genes)) %>%
    #    # y order: object then celltype (grouped)
    #    arrange(object, celltype) %>%
    #    mutate(row_label = factor(row_label, levels = unique(row_label)))
    df$proportion = df[,col, drop = TRUE]
    print(df)

    ggplot(df, aes(x = SampleID, y = proportion, fill = cell_class)) +
        geom_col() +
        theme_classic(base_size = 12) +
        theme(
            axis.title = element_blank(),
            axis.text.x = element_text(angle = x_text_angle, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 9),
            panel.grid.major = element_line(linewidth = 0.2, color = "grey92")
        )
}


pdf(pdf_file, width = 8, height = 6)
plot_proportions(df, "celltype_proportion")
plot_proportions(photoreceptors, "pr_proportion")
dev.off()
