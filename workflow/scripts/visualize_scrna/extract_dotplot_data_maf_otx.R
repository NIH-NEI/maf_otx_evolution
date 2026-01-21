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

print(seurat_file)
print(tsv_file)
print(organism)

genes_to_plot <- (
    read.delim("configs/hahn_maf_otx.tsv")
    |> as.list()
)
#print(genes_to_plot)
print(genes_to_plot[[organism]])

seurat_obj <- readRDS(seurat_file)


# -----------------------------
# Helper: extract dotplot stats for ONE object, while skipping missing genes
# -----------------------------
dotplot_stats_one_object <- function(
        seu,
        genes,
        celltype_col = "celltype",
        object_name = "obj1",
        assay = NULL,
        slot = "data"   # "data" (log-normalized), or "counts" if you prefer
) {
    if (!is.null(assay)) DefaultAssay(seu) <- assay
    # --- Fix duplicate cell names (most common cause of your error)
    if (anyDuplicated(Cells(seu)) > 0) {
      seu <- RenameCells(seu, new.names = make.unique(Cells(seu)))
    }

    ## --- Guard against duplicate feature names
    #if (anyDuplicated(rownames(seu)) > 0) {
    #  #stop("This Seurat object has duplicate feature (gene) names in rownames(seu). Fix upstream before DotPlot().")
    #}
    ## genes that actually exist in this object
    genes_present <- intersect(genes, rownames(seu))
    if (length(genes_present) == 0) {
        return(tibble())  # nothing to contribute
    }

    print(genes_present)
    # This line is important because chicken seurat object has
    #  NA in idents but not in cell_class column
    Idents(seu) <- "cell_class"

    print(colnames(seu@meta.data))
    print(Idents(seu) %>% anyNA())
    print(seu@meta.data |> filter(is.na(cell_class)))
    print("Going inside DotPlot ....")

    # Use Seurat's DotPlot to compute avg.exp and pct.exp consistently
    p <- DotPlot(seu, features = genes_present, group.by = celltype_col, assay = assay)
    df <- p$data %>%
        as_tibble() %>%
        transmute(
            object = object_name,
            celltype = as.character(id),
            gene = as.character(features.plot),
            avg_exp = avg.exp,
            pct_exp = pct.exp
        )

    df
}

df <- dotplot_stats_one_object(
    seu = seurat_obj,
    genes = genes_to_plot[[organism]],
    celltype_col = "cell_class",
    object_name = organism,
    assay = NULL,
    slot = "data"
)

print(df)
readr::write_tsv(df, tsv_file)
quit()

# -----------------------------
# Main: build combined stats table for all objects
# -----------------------------
build_mega_dotplot_df <- function(
        seurat_list,          # named list of Seurat objects
        genes,
        celltype_col = "celltype",
        assay = NULL,
        slot = "data"
) {
    imap_dfr(seurat_list, ~dotplot_stats_one_object(
        seu = .x,
        genes = genes,
        celltype_col = celltype_col,
        object_name = .y,
        assay = assay,
        slot = slot
    )) %>%
        # optional: keep a stable ordering
        mutate(
            object = factor(object, levels = names(seurat_list)),
            # row label with hierarchy
            row_label = paste(object, celltype, sep = " | ")
        )
}

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

# -----------------------------
# Example usage
# -----------------------------
# seurat_list must be a named list:
seurat_list <- list(
    human = readRDS("/vf/users/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Human_initial.rds")
    , crabMacaque = readRDS("/vf/users/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Macaque_initial.rds")
    , marmoset = readRDS("/vf/users/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Marmoset_initial.rds")
    , labMouse = readRDS("/vf/users/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Mouse_initial.rds")
    , grassMouse = readRDS("/vf/users/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Rhabdomys_initial.rds")
    , deerMouse = readRDS("/vf/users/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Peromyscus_initial.rds")
    , groundSquirrel = readRDS("/vf/users/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Squirrel_initial.rds")
    , treeShrew = readRDS("/vf/users/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Tree_shrew_initial.rds")
    , ferret = readRDS("/vf/users/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Ferret_initial.rds")
    , sheep = readRDS("/vf/users/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Sheep_initial.rds")
    , cow = readRDS("/vf/users/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Cow_initial.rds")
    , opossum = readRDS("/vf/users/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Opossum_initial.rds")
    , chicken = readRDS("/vf/users/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Chicken_initial.rds")
    , brownAnole = readRDS("/vf/users/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Lizard_initial.rds")
    , zebrafish = readRDS("/vf/users/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Zebrafish_initial.rds")
    , lamprey = readRDS("/vf/users/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Lamprey_initial.rds")
)

genes <- c("MAFA", "MAFB", "MAF", "NRL", "MAFF", "MAFG", "MAFK", "CRX", "OTX1", "OTX2")

df <- build_mega_dotplot_df(
  seurat_list = seurat_list,
  genes = genes,
  celltype_col = "cell_class",  # change to your metadata column
  assay = "RNA"
)
p <- plot_mega_dotplot(df, genes)
p














# ==================



seurat_file <- "/vf/users/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Mouse_initial.rds"
#seurat_file <- "/vf/users/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Lamprey_initial.rds"
seurat_obj <- readRDS(seurat_file)

genes <- data.frame(
    gene = rownames(seurat_obj)
)

head(genes)

head(seurat_obj@meta.data)

meta_cols <- colnames(seurat_obj@meta.data)
meta_cols


mafs <- grep("maf", ignore.case = TRUE, value = TRUE, rownames(seurat_obj))
mafs

otxs <- grep("otx", ignore.case = TRUE, value = TRUE, rownames(seurat_obj))
otxs

grep("NRL", ignore.case = TRUE, value = TRUE, rownames(seurat_obj))

grep("CRX", ignore.case = TRUE, value = TRUE, rownames(seurat_obj))

# Define the list of genes you want to visualize
genes_to_plot <- list(
    human = c("MAFA", "MAFB", "MAF", "NRL", "MAFF", "MAFG", "MAFK", "CRX", "OTX1", "OTX2"),
#    mouse = c("MAFA", "MAFB", "MAF", "NRL", "MAFF", "MAFG", "MAFK", "CRX", "OTX1", "OTX2"),
    mouse = c("MAFA", "MAFB", "MAF", "NRL", "CRX", "OTX1", "OTX2"),
#    mouse = c("Mafa", "Mafb", "Maf", "Nrl", "Maff", "Mafg", "Mafk", "Crx", "Otx1", "Otx2"),
    lamprey = c("MAFA", "MAFB", "MAF", "NRL", "MAFF", "MAFG", "MAFK", "CRX", "OTX1", "OTX2", "OTX5"),
    extra = c()
)

dot_plot_org <- function(organism) {

    genes_to_plot <- genes_to_plot[[organism]]

    # Check if the genes are present in the Seurat object
    missing_genes <- setdiff(genes_to_plot, rownames(seurat_obj))
    if (length(missing_genes) > 0) {
        warning("The following genes are not present in the Seurat object: ", paste(missing_genes, collapse = ", "))
    }

    # Filter out missing genes
    genes_to_plot <- intersect(genes_to_plot, rownames(seurat_obj))

    # Generate the dot plot
    if (length(genes_to_plot) > 0) {
        dotplot <- DotPlot(seurat_obj, features = genes_to_plot, group.by = "cell_class") +
            scale_color_gradient(low = "blue", high = "red") +
            theme_minimal() +
            labs(title = "Dot Plot of Genes Across Cell Types",
                 y = "Cell Types",
                 x = "Genes")

        # Save the plot as a PDF

        # Display the plot
        print(dotplot)
    } else {
        message("No valid genes to plot.")
    }

}

pdf("dotplot_genes_celltypes.pdf", width = 8, height = 6)

dot_plot_org("lamprey")
dot_plot_org("mouse")
dot_plot_org("human")

dev.off()
