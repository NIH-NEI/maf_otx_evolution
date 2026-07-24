# --------------------------------------------------------------
# R script to reproduce the phylogenetic tree style in your image
# Packages: ape, ggtree, ggplot2, treeio (optional)
# --------------------------------------------------------------

# Install if needed (uncomment first time)
# install.packages(c("ape", "ggtree", "ggplot2"))
# BiocManager::install("ggtree")   # if ggtree is not installed yet

library(ape)
library(tidyverse)
library(ggtree)
library(rphylopic)
library(ggtreeExtra)
library(ggnewscale)


datasets <- c("Bulk" = "firebrick", "SC" = "green", "Annotated" = "blue", "No" = "grey95")

colors <- (
    read.delim("configs/species_colors.txt")
)
print(colors)

org_groups = unique(colors$OrganismColor)

colors <- (
  colors |> mutate(OrganismColor = factor(OrganismColor, levels = org_groups, ordered = TRUE))
)


tree <- read.tree("imports/trees/species181_topology_with_scientific_names.nwk")
tree$tip.label <- sub("\\|.*", "", tree$tip.label)


df <- (
  read.delim("configs/species1XX_table.tsv")
  #|> left_join(read.delim("configs/species173_annotations.tsv") |> select(OrganismShortName, AnnotationSource))
  |> mutate(
    taxa = OrganismShortName
  )
  |> mutate(bulk = ifelse(bulkrna107 == "Yes", "Bulk", "No"))
  |> mutate(sc = ifelse((liscrna24 == "Yes") | (hahnscrna17 == "Yes"), "SC", "No"))
  |> mutate(annotated = ifelse(AnnotationSource %in% c("GCF", "Ensembl", "Custom"), "Annotated", "No"))
  |> relocate(taxa)
)

# Step 5: Create base ggtree plot
p <- ggtree(tree, layout = "fan", open.angle=5, ladderize = FALSE)  # Change layout if desired (e.g., "rectangular,circular")
print(p)

# Step 6: Add tip labels colored by group
p <- (
  p %<+% df
)
print(p)

mylevels <- c("Bulk" = "firebrick", "SC" = "green", "Annotated" = "blue", "No" = "grey95")


dat <- (
  df[, c("taxa", "bulk", "sc", "annotated")]
  |> mutate(bulk = factor(bulk, levels = names(mylevels)))
  |> mutate(sc = factor(sc, levels = names(mylevels)))
  |> column_to_rownames("taxa")
)

dat_long <- (
  df[, c("taxa", "bulk", "sc", "annotated")]
  |> pivot_longer(
    cols = -taxa,              # all columns except label
    names_to = "variable",      # heatmap column name
    values_to = "value"         # cell value
  )
  |> mutate(value = factor(value, levels = names(mylevels)))
)

dat2 <- (
  df[, c("taxa", "OrganismColor")]
  |> pivot_longer(
    cols = -taxa,              # all columns except label
    names_to = "variable",      # heatmap column name
    values_to = "value"         # cell value
  )
  |> mutate(value = factor(value, levels = org_groups))
)

p3 <- (
  p
  + new_scale_fill()
  + geom_fruit(
    data    = dat2,
    geom    = geom_tile,
    mapping = aes(x = variable, y = taxa, fill = value),
    offset  = 0.01,   # distance from tree to heatmap
    pwidth  = 0.5    # relative width of the heatmap band
  )
  + scale_fill_manual(values=pull(colors, HexCode2, OrganismColor), name="Groups")
)

print(p3)

p4 <- (
  p3
  + geom_tiplab(aes(color = OrganismColor), size = 2, align = TRUE, offset = 1.0, show.legend = FALSE)
  + scale_color_manual(values = pull(colors, HexCode2, OrganismColor))
  #  + theme(legend.position = "none")
)

dat <- (
  df[, c("taxa", "bulk", "sc", "annotated")]
  |> pivot_longer(
    cols = -taxa,              # all columns except label
    names_to = "variable",      # heatmap column name
    values_to = "value"         # cell value
  )
  |> mutate(value = factor(value, levels = names(mylevels)))
)

p2 <- (
  p4
  + new_scale_fill()
  + geom_fruit(
    data    = dat_long,
    geom    = geom_tile,
    mapping = aes(x = variable, y = taxa, fill = value),
    offset  = 0.25,   # distance from tree to heatmap
    pwidth  = 0.05    # relative width of the heatmap band
  )
  + scale_fill_manual(values=mylevels, name="Dataset")
  + new_scale_fill()
)
print(p2)

pdf("species181.pdf", height=8, width=10)
print(p2)
dev.off()
