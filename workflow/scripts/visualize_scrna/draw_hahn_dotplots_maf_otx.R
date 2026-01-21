library(Seurat)
library(tidyverse)

seurat_file <- "/vf/users/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Mouse_initial.rds"
#seurat_file <- "/vf/users/pals2/celltype-evolution/data/ShekharData/SeuratObjects/Lamprey_initial.rds"

seurat_file <- snakemake@input[["seu"]]
pdf_file <- snakemake@output[[1]]
organism <- snakemake@wildcards[["org"]]

print(seurat_file)
print(pdf_file)
print(organism)

genes_to_plot <- (
    read.delim("configs/hahn_maf_otx.tsv")
    |> as.list()
)
print(genes_to_plot)

seurat_obj <- readRDS(seurat_file)

genes <- data.frame(
    gene = rownames(seurat_obj)
)
head(genes)

head(seurat_obj@meta.data)

mafs <- grep("maf", ignore.case = TRUE, value = TRUE, rownames(seurat_obj))
mafs

otxs <- grep("otx", ignore.case = TRUE, value = TRUE, rownames(seurat_obj))
otxs

grep("NRL", ignore.case = TRUE, value = TRUE, rownames(seurat_obj))
grep("CRX", ignore.case = TRUE, value = TRUE, rownames(seurat_obj))


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

pdf(pdf_file, width = 8, height = 6)
dot_plot_org(organism)
dev.off()
