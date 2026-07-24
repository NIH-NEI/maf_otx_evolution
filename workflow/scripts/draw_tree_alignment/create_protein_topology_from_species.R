library(ape)
library(tidyverse)

# ---- INPUT FILE PATHS ----
subset_fasta_path <- "subset.fasta"
original_tree_path <- "original_tree.nwk"
output_tree_path <- "subset_tree.nwk"

subset_fasta_path <- snakemake@input[["aln"]]
original_tree_path <- snakemake@input[["tree"]]
output_tree_path <- snakemake@output[["tree"]]

print(subset_fasta_path)
print(original_tree_path)
print(output_tree_path)



pdf("Rplots.pdf", height=15, width=10)

# ---- READ INPUTS ----
cat("Reading original tree...\n")
tree <- read.tree(original_tree_path)
plot(tree)

cat("Reading full and subset FASTA files...\n")
subset_fasta <- read.FASTA(subset_fasta_path)

# ---- PROCESSING ----
subset_labels <- names(subset_fasta)
print(subset_labels)

df <- (
    tibble(taxa = subset_labels)
    |> separate(taxa, into = c("ShortName", "GeneSymbol"), sep = "__", remove = FALSE)
)
print(df)

# Identify tips to drop
tips_to_drop <- setdiff(tree$tip.label, df$ShortName)
cat("Dropping", length(tips_to_drop), "tips from the tree...\n")

# Prune tree
tree <- drop.tip(tree, tips_to_drop)
print(tree)
plot(tree)

# ---- Step 2: For each duplicated ShortName in tree, add tips ----
for (sn in unique(df$ShortName)) {
  n_taxa <- sum(df$ShortName == sn)
  has_node_labels <- (sum(is.na(tree$node.label)) == 0)
  if (n_taxa > 1 && sn %in% tree$tip.label) {
    # Find the node to split
    tree <- bind.tree(
      x = tree,
      y = stree(n = n_taxa, tip.label = df$taxa[df$ShortName == sn]),
      where = which(tree$tip.label == sn)
    )
    plot(tree)
    title(paste("Before dropping", sn))
    nodelabels()
    print(tree$node.label)
    if (has_node_labels) {
      tree$node.label[is.na(tree$node.label)] <- sn
    }
    print(tree$node.label)
    # Drop the old tip
    tree <- drop.tip(tree, sn)
    plot(tree)
    title(paste("After dropping", sn))
  } else {
      tree$tip.label[tree$tip.label == sn] <- df$taxa[df$ShortName == sn]
      plot(tree)
    title(paste("Renamed", sn))
  }
}
plot(tree)

if ((length(setdiff(df$taxa, tree$tip.label)) > 0) || 
    (length(setdiff(tree$tip.label, df$taxa)) > 0) ) {
    print("Error!!!! Not done")
    quit()
}

# ---- OUTPUT ----
cat("Writing pruned tree to", output_tree_path, "...\n")
write.tree(tree, output_tree_path)

cat("Done! Pruned tree has", length(tree$tip.label), "tips.\n")

