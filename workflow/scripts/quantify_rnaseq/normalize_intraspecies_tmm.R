library(tidyverse)
library(tximport)
library(readr) # for read_tsv
library(edgeR)

#samples_file <- "experiment_samples.tsv"
#ortholog_table_file <- "ortholog_table.tsv"
#ortholog_expr_file <- "ortho_expression_quantized.csv"

samples_file <- snakemake@input[["samples"]]
ortholog_expr_file <- snakemake@output[["norm"]]
dge_rds_file <- snakemake@output[["dge"]]
factors_file <- snakemake@output[["factor"]]

normalization <- "tmm_tpm_edger"

if (normalization == "tmm_tpm_edger") {
    expr_col = "TPM"
} else if (normalization == "tmm_abundance_edger") {
    expr_col = "abndance"
} else {
    stop(paste("Unknown normalization", normalization))
}

print(samples_file)
print(ortholog_expr_file)
print(dge_rds_file)

samples <- (
    stack(snakemake@input)
    |> rename(SampleID := ind, quant_file := values)
    |> filter(SampleID != "") # remove non-named entries created by snakemake
    |> filter(!SampleID %in% c("ortho", "samples"))
    |> left_join(read.delim(samples_file))
)
samples

quantize_sample <- function(species, sample, quant_file) {
    expr <- (
        read.delim(quant_file)
        |> select(gene_id, !!expr_col)
    )
    colnames(expr) <- c("gene_id", sample)
    print(expr[1:5,])
    expr
}

expr <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE),
               lapply(1:nrow(samples), function(i) {
                   quantize_sample(samples$OrganismID[i], samples$SampleID[i], samples$quant_file[i])
               }))

expr[1:5,]

expr <- (
    expr
    |> column_to_rownames("gene_id")
    |> as.matrix()
)


#### TMM Normalize the data
#Gene-level normalization
cat("Going to create DGEList ...\n")
na_positions <- is.na(expr)
na_positions[1:5,]
expr[na_positions] <- 0
gene.dge <- DGEList(expr, samples = samples, genes = row.names(expr))
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
