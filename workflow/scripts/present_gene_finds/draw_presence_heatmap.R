library(ComplexHeatmap)
library(circlize)
library(tidyr)
library(dplyr)

infile <- "exports/curated_genes/draft_combined_maf_otx.tsv"
infile <- "scratch/grand_list/grand_list.tsv"

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

df <- (
    read.delim(infile)
    |> filter(GeneGroup != "Ignore")
    |> mutate(GeneGroup = ifelse(GeneGroup == "MAFL", "CMAF", GeneGroup))
    |> mutate(GeneGroup = ifelse(GeneGroup == "MAFS", "MAFF", GeneGroup))
    |> mutate(GeneGroup = ifelse(GeneGroup == "OTX", "OTX1", GeneGroup))
#    |> select(OrganismShortName, OrganismColor, ProteinSource, ProteinID, GeneGroup)
#    |> select(-OrganismID, -AnnotationSource, -bulkrna107)
#    |> select(-nMAFL, -nMAFS, -nOTX)
#    |> pivot_longer(
#        cols = c(-OrganismShortName, -OrganismColor, -AnnotationType),
#        names_to = "GeneGroup",
#        values_to = "Presence"
#    )
    |> group_by(OrganismShortName, OrganismColor, ProteinSource, GeneGroup)
    |> summarize(Presence = n_distinct(ProteinID))
#    |> mutate(Count = ifelse(AnnotationSource == "Annotated", 1,
#                    ifelse(AnnotationSource == "Denovo", 2, 0)))
    |> mutate(Count = ifelse(Presence == 0, 0,
                      ifelse(ProteinSource == "Annotated", 1,
                      ifelse(ProteinSource == "Denovo", 2, 0))))
    |> group_by(OrganismShortName, OrganismColor, GeneGroup)
    |> summarize(Presence = sum(Count))
    |> mutate(Presence = recode(Presence, `0` = "None", `1` = "Annotated", `2` = "Denovo", `3` = "Both"))
    |> pivot_wider(
        names_from = GeneGroup,
        values_from = Presence,
        values_fill = "None"
    )
    |> left_join(
        read.delim("configs/species173_table.tsv")
        |> select(OrganismShortName, OrganismOrder173, bulkrna107)
    )
    |> arrange(OrganismOrder173)
    |> ungroup()
    |> mutate(OrganismColor = shortColor[OrganismColor])
)
head(df)


gene_grouping = c(
    "CMAF" = "Large MAF",
    "MAFA" = "Large MAF",
    "MAFB" = "Large MAF",
    "NRL" = "Large MAF",
    "MAFF" = "Small MAF",
    "MAFG" = "Small MAF",
    "MAFK" = "Small MAF",
    "OTX1" = "OTX",
    "OTX2" = "OTX",
    "CRX" = "OTX"
)

draw_heatmap <- function(df, outfile) {

mat <- (
    df |> select(-OrganismShortName, -OrganismColor, -OrganismOrder173, -bulkrna107)
    |> as.matrix()
)
rownames(mat) <- pull(df, OrganismShortName)
mat <- mat[,names(gene_grouping)]
mat[1:5,]

organism_group <- (
    df
    |> mutate(OrganismGroup = factor(df$OrganismColor, levels = names(org_colors), ordered = TRUE))
    |> mutate(HasRNAseq = factor(df$bulkrna107, levels = c("Yes", "No"), ordered = TRUE))
    |> select(OrganismGroup, HasRNAseq)
)
row.names(organism_group) <- rownames(mat)

gene_group <- (
    data.frame(gene = colnames(mat))
    |> mutate(GeneGroup = gene_grouping[gene])
)
gene_group <- factor(gene_group$GeneGroup, levels = c("Large MAF", "Small MAF", "OTX"), ordered = TRUE)
names(gene_group) <- colnames(mat)

mat[is.na(mat)] <- "None"  # or keep NA if you prefer gray

# ===================================
# 2. Create annotations
# ===================================

# Row annotation: OrganismGroup
ra <- rowAnnotation(
  df = organism_group,
  col = list(
        OrganismGroup = org_colors,
        HasRNAseq = c("Yes" = "#4BCAAD", "No" = "#CE6633")
    ),
  annotation_name_gp = gpar(fontsize = 4),
  width = unit(4, "mm")
)

# Column annotation: GeneGroup
top_ann <- HeatmapAnnotation(
  GeneGroup = gene_group,
  col = list(GeneGroup = c(
    "Large MAF" = "#8dd3c7",
    "Small MAF" = "#ffffb3",
    "OTX" = "#bebada"
  )),
  annotation_name_gp = gpar(fontsize = 4),
  height = unit(4, "mm")
)

# ===================================
# 3. Make the heatmap
# ===================================

ht <- Heatmap(
  mat,
  name = "Protein presence",
  col = c(
    "Annotated" = "#DCEDC8", #"#99F9FF", # annotated
    "Denovo" = "#8BC34A", #"#FFCA99", # denovo
    "Both" = "#33691E", #"#B66DFF" # both
    "None" = "#FFCCBC" #"grey99", # absent
  ),

  # Row settings
  cluster_rows = FALSE,
  show_row = TRUE,
  row_split = organism_group$OrganismGroup,           # splits & orders by OrganismGroup
  row_title = "%s",                     # shows group name
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 2),

  # Column settings
  cluster_columns = FALSE,
  column_split = gene_group,            # group genes by functional category
  column_title = "Genes",
  column_title_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 0,
  column_names_centered = TRUE,
  top_annotation = top_ann,

  # General appearance
  show_row_names = TRUE,
  show_column_names = TRUE,            # set TRUE if <100 genes
#  heatmap_legend_param = list(
#    title = "Presence",
#    at = c(0, 1),
#    labels = c("Absent", "Present")
#  ),

  # Add row annotation on the left
  left_annotation = ra,

  column_gap = unit(0.2, "mm"),
  row_gap = unit(0.2, "mm")

#  # Size
#  width = unit(12, "cm"),
#  height = unit(16, "cm")
)

# Draw it
pdf(outfile, height = 8, width = 10.25)
draw(ht, merge_legends = TRUE)
dev.off()

}

draw_heatmap(df, "presence_heatmap.pdf")

df <- (
    df
    |> mutate(across(where(is.character), ~ ifelse(. %in% c("Denovo"), "None", .)))
    |> mutate(across(where(is.character), ~ ifelse(. %in% c("Both"), "Annotated", .)))
)

head(df)
draw_heatmap(df, "presence_heatmap_annotated.pdf")
