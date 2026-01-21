# Load libraries
library(GenomicFeatures)
library(tximport)


# Path to your NCBI RefSeq GFF file (e.g., GCF_XXXXX_genomic.gff or similar)
# Path to your Ensembl annotation file (usually .gtf.gz, e.g., Homo_sapiens.GRCh38.110.gtf.gz)
gff_file <- snakemake@input[[1]]

# output is tab separated file
out_file <- snakemake@output[[1]]

# Create a TxDb object from the GFF
txdb <- makeTxDbFromGFF(file = gff_file,
                        format = "auto")  # Detects GFF3 or GTF automatically

# Extract all transcript names
k <- keys(txdb, keytype = "TXNAME")

# Create tx2gene: column 1 = TXNAME (transcript ID), column 2 = GENEID (gene ID)
tx2gene <- select(txdb, keys = k, columns = "GENEID", keytype = "TXNAME")

# Reorder columns if needed (tximport expects transcript first, then gene)
tx2gene <- tx2gene[, c("TXNAME", "GENEID")]

## Remove version suffixes (e.g., .1 in NM_123.1 or ENST00000380152.7) if your Salmon/Kallisto quant.sf uses bare IDs
#tx2gene$TXNAME <- gsub("\\.[0-9]+$", "", tx2gene$TXNAME)  # Strip transcript version
#tx2gene$GENEID <- gsub("\\.[0-9]+$", "", tx2gene$GENEID)  # Strip gene version (ENSG...)

# Save for reuse
write.table(tx2gene, out_file, sep = "\t", row.names = FALSE, quote = FALSE)

