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

get_uid <- function(sp) {
  tryCatch(
    get_uuid(name = sp),
    error = function(msg) {
      message(paste("Error for list member:", sp, "there is no phylopic\n"))
      NA
    }
  )
}

df <- (
  read.delim("configs/species173_table.tsv")
  |> mutate(
    PhylopicUID = map(OrganismScientificName, get_uid)
  )
  |> select(OrganismID, PhylopicUID)
)

head(df)

df <- unnest(df, cols = c(PhylopicUID))

write.table(df, "configs/species173_phylopic.tsv", sep = "\t", row.names = FALSE)
