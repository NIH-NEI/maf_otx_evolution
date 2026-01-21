library(tidyverse)
library(tximport)
library(readr) # for read_tsv
library(edgeR)

under_snakemake <- TRUE

if (under_snakemake) {
samples_file <- snakemake@input[["samples"]]
species_file <- snakemake@input[["orgs"]] # "configs/species50_table.tsv"
ortholog_table_file <- snakemake@input[["ortho"]]
ortholog_expr_file <- snakemake@output[["norm"]]
dge_rds_file <- snakemake@output[["dge"]]
factors_file <- snakemake@output[["factor"]]

normalization <- snakemake@wildcards[["norm"]]

print(ortholog_table_file)
print(samples_file)
print(ortholog_expr_file)
print(dge_rds_file)
print(normalization)

samples <- (
    stack(snakemake@input)
    |> rename(SampleID := ind, quant_file := values)
    |> filter(SampleID != "") # remove non-named entries created by snakemake
    |> filter(!SampleID %in% c("ortho", "samples", "orgs"))
    |> left_join(read.delim(samples_file) |> select(SampleID, OrganismID))
    |> left_join(read.delim(species_file) |> select(OrganismShortName, OrganismID, SurrogateShortName))
)
head(samples)
write.table(samples, "sample_for_debug.tsv", sep = "\t")
} else {
    samples <- read.delim("sample_for_debug.tsv")

    ortholog_table_file <- "scratch/orthofinder_surrogate_bulkrna/ortholog_tables/orthotab_maf_focused.tsv"
    normalization <- "tmm_abundance_edger"
}

if (normalization == "tmm_tpm_edger") {
    expr_col = "TPM"
} else if (normalization == "tmm_abundance_edger") {
    expr_col = "abndance"
} else {
    stop(paste("Unknown normalization", normalization))
}

#samples <- samples[1:10,]

ortholog_table <- (
    read.delim(ortholog_table_file)
)
head(ortholog_table)
tail(ortholog_table)

quantify_sample <- function(species, surrogate, quant_file) {
    message(paste("Processing", species, surrogate, quant_file, "...."))
    expr <- (
        read.delim(quant_file)
        |> column_to_rownames("gene_id")
    )
    this_species = (
        ortholog_table
        |> rename(symbol := !!surrogate)
        |> select(OrthoSymbol, symbol)
        |> separate_rows(symbol, sep = ",")
    )
    this_species$expr = expr[this_species$symbol, expr_col]
    #print(this_species)
    ret_val <- (
        this_species
        |> group_by(OrthoSymbol)
        |> summarize(expr = median(expr))
        |> column_to_rownames("OrthoSymbol")
    )
    #print(head(ret_val))
    ret_val |> select(expr)
}

aligned_cbind <- function(dflist) {
  # union of all rownames
  all_rows <- Reduce(union, lapply(dflist, rownames))
  print(all_rows)

  # reorder each df to have the same rows (in same order)
  dflist_aligned <- lapply(dflist, function(df) {
    df[all_rows, , drop = FALSE]
  })

  # now cbind works properly
  do.call(cbind, dflist_aligned)
}

dfs <- lapply(1:nrow(samples), function(i) {
  quantify_sample(samples$OrganismShortName[i],
                  samples$SurrogateShortName[i],
                  samples$quant_file[i])
})

# Any NA rownames?
sapply(dfs, function(df) sum(is.na(rownames(df))))

# Any duplicated rownames?
sapply(dfs, function(df) any(duplicated(rownames(df))))

# Show the worst offenders
which.max(sapply(dfs, function(df) sum(is.na(rownames(df)))))


expr <- aligned_cbind(dfs)
#rownames(expr) <- row.names(ortholog_table)
colnames(expr) <- samples$SampleID


dim(expr)
dim(ortholog_table)
setdiff(row.names(expr), ortholog_table$OrthoSymbol)
setdiff(ortholog_table$OrthoSymbol,row.names(expr))


#### TMM Normalize the data
#Gene-level normalization
cat("Going to create DGEList ...\n")
na_positions <- is.na(expr)
na_positions[1:5,1:5]
expr[na_positions] <- 0
gene.dge <- DGEList(expr, samples = samples, genes = ortholog_table)
cat("Going to calculate norm factors ...\n")
gene.dge <- calcNormFactors(gene.dge)
gene.dge$cpm <- cpm(gene.dge)
gene.dge$lcpm <- log2(gene.dge$cpm + 1)
gene.dge$napos <- na_positions
gene.dge$cpm[na_positions] <- NA

head(gene.dge$lcpm)
head(gene.dge$cpm)
head(gene.dge[[2]])

factors <- (
    gene.dge[[2]]
    |> select(group, lib.size, norm.factors)
    |> rownames_to_column("SampleID")
)
head(factors)

saveRDS(gene.dge, dge_rds_file)
write.table(gene.dge$cpm, ortholog_expr_file, row.names = TRUE, sep = "\t")
write.table(factors, factors_file, row.names = FALSE, sep = "\t")
