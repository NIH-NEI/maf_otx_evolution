#library(TDbook)
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(viridis)
library(ggnewscale)
library(ape)

snakemake@source("mymsaplot.R")

tree_file <- snakemake@input[["tree"]]
msa_fasta_file <- snakemake@input[["aln"]]
seq_info_file <- snakemake@input[["info"]]
#sp_info_file <- snakemake@input[["spgrp"]]
#outgroup_file <- snakemake@input[["outgroup"]]
pdf_file <- snakemake@output[["plot"]]
final_tree_file <- snakemake@output[["tree"]]
final_itol_file <- snakemake@output[["itol"]]
final_fasta_file <- snakemake@output[["fasta"]]
final_tip_file <- snakemake@output[["tip"]]

tree_file
seq_info_file
pdf_file

colors <- (
    read.delim("configs/species_colors.txt")
)
print(colors)

org_groups = unique(colors$OrganismColor)

colors <- (
  colors |> mutate(OrganismColor = factor(OrganismColor, levels = org_groups, ordered = TRUE))
)

gene_group_colors = c(
    "MAFL"="#7F0000", "CMAF"="#B22222", "NRL"="#E41A1C", "MAFA"= "#FF7F00", "MAFB"= "#FDB863",
    "MAFS"= "#08306B", "MAFF"= "#2171B5", "MAFG"= "#6BAED6", "MAFK"= "#9ECAE1",
    "OTX"= "#00441B", "OTX1"= "#238B45", "OTX2"= "#66C2A4", "CRX"= "#B2DF8A",
    "Outgroup"= "gray"
)


tree_seq_nwk <- read.tree(tree_file)
print(tree_seq_nwk)

seq_info <- (
    read.delim(seq_info_file)
#    |> left_join(read.delim(sp_info_file))
)

outgroups = seq_info |> filter(GeneGroup == 'Outgroup') |> pull(ProteinID) |> head(1)
print(outgroups)

if (length(outgroups) > 0) {
    tree_seq_nwk <- root(tree_seq_nwk, outgroup = outgroups[outgroups %in% tree_seq_nwk$tip.label])
}
tree_seq_nwk$tip.label

seq_info <- (
    data.frame(ProteinID = tree_seq_nwk$tip.lab)
    ## ALWAYS MAKE SURE THERE IS A COLUMN NAMED taxa
    |> mutate(taxa = ProteinID)
    |> left_join(seq_info)
)
seq_info |> select(taxa, ProteinID, OrganismColor, GeneGroup, Status)

gene_groups <- unique(seq_info$GeneGroup)
statuses <- unique(seq_info$Status)

fasta_obj <- read_validate_fasta(msa_fasta_file, tree_seq_nwk)

name2taxa <- seq_info |> pull(taxa, ProteinID)
names(fasta_obj) <- name2taxa[names(fasta_obj)]

tree_seq_nwk$tip.label <- unname(name2taxa[tree_seq_nwk$tip.label])
names(seq_info)

p <- (
    ggtree(tree_seq_nwk)
    %<+% seq_info
)

current_offset <- 0

current_offset <- current_offset + 1
p <- mymsaplot(p, fasta_obj, offset=current_offset, width=3)

dat1 <- (
  seq_info
  |> select(taxa, GeneGroup)
  |> pivot_longer(
    cols = -taxa,              # all columns except label
    names_to = "variable",      # heatmap column name
    values_to = "value"         # cell value
  )
  |> mutate(value = factor(value, levels = gene_groups))
)

p <- (
  p
  + new_scale_fill()
  + geom_fruit(
    data    = dat1,
    geom    = geom_tile,
    mapping = aes(x = variable, y = taxa, fill = value),
    offset  = 0.5,   # distance from tree to heatmap
    pwidth  = 0.3    # relative width of the heatmap band
  )
  + scale_fill_manual(values=gene_group_colors, name="Groups")
)

dat2 <- (
  seq_info
  |> select(taxa, Status)
  |> pivot_longer(
    cols = -taxa,              # all columns except label
    names_to = "variable",      # heatmap column name
    values_to = "value"         # cell value
  )
  |> mutate(value = factor(value, levels = statuses))
)

p <- (
  p
  + new_scale_fill()
  + geom_fruit(
    data    = dat2,
    geom    = geom_tile,
    mapping = aes(x = variable, y = taxa, fill = value),
    offset  = 0.05,   # distance from tree to heatmap
    pwidth  = 0.3    # relative width of the heatmap band
  )
    + scale_fill_manual(values=c("Final" = "green", "Undecided" = "red"), name="Status")
)

p <- (
    p + geom_tiplab(aes(color = OrganismColor), size=3, align=TRUE, offset=0.1)
    + scale_color_manual(values=pull(colors, HexCode2, OrganismColor), name="OrganismGroups")
)

ntaxa <- length(tree_seq_nwk$tip.label)
pdf(pdf_file, width=40, height=max(10, ntaxa/4))
print(p)
dev.off()

