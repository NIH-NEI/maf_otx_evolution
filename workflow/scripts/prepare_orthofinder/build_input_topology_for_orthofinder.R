#!/usr/bin/env Rscript
library(ape)
library(tidyverse)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop(paste("Need three input files. Example: Rscript this_script.R orglist.txt in_tree.nwk out_tree.nwk"))
}
org_file <- args[1]
infile <- args[2]
outfile <- args[3]

print(org_file)
print(infile)
print(outfile)


org_list <-(
    read.delim(org_file, header = FALSE)
    |> pull(V1)
)
print(org_list)

pdf("Rplots.pdf", height=50, width=20)

tree <- read.tree(infile)
plot(tree)
print(tree$tip.label)
tree$tip.label <- sub("\\|.*$", "", tree$tip.label)
print(tree$tip.label)
plot(tree)
tree <- drop.tip(tree, setdiff(tree$tip.label, org_list))
print(tree$tip.label)
plot(tree)
dev.off()
stopifnot(length(setdiff(org_list, tree$tip.label))==0)
stopifnot(length(setdiff(tree$tip.label, org_list))==0)
write.tree(tree,outfile)

tree<-read.tree(outfile)
stopifnot(length(setdiff(org_list, tree$tip.label))==0)
stopifnot(length(setdiff(tree$tip.label, org_list))==0)
