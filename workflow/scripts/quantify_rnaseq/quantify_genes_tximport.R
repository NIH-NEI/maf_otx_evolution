library(tidyverse)
library(tximport)
library(readr) # for read_tsv
library(edgeR)

seq2gene_file <- snakemake@input[["seq2gene"]]
quants_file <- snakemake@input[["quants"]]
output_file <- snakemake@output[["gene"]]
print(output_file)

seq2gene <- (
    read.delim(seq2gene_file)
    |> select(cds_header, gene)
)

stopifnot(sum(is.na(seq2gene)) == 0)

#Import trans data
trans <- tximport(quants_file, txOut=T, type = "kallisto", importer = read_tsv, tx2gene = NULL)
print(names(trans))

#Summarize to gene
gene <- summarizeToGene(trans, seq2gene, countsFromAbundance = "lengthScaledTPM")
print(paste0("Abundance was estimated for ", dim(gene$counts)[1], " genes in this annotation."))

df = cbind(gene$abundance, gene$counts) #, gene$infReps[[1]], gene$length)
colnames(df) <- c("abndance", "TPM")
df <- (
    as.data.frame(df)
    |> rownames_to_column("gene_id")
)
head(df)
print(output_file)
write.table(df, output_file, row.names=FALSE, sep="\t")
