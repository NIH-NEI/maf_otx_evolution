library(Biostrings)

path2fasta <- "../eye_evolution/code/imports/ForGENESPACE/rawGenomes/human/GCF_000001405.40_GRCh38.p14_translated_cds.faa.gz"
path2fasta <- snakemake@input[["pep"]]
out_file <- snakemake@output[[1]]


troubleShoot <- TRUE
headerStripText <- "gene=|\\[|\\]"
headerEntryIndex <- 2
headerSep <- " "
convertSpecialCharacters <- "_"

print(path2fasta)

# -- read in the peptide annotations
fa <- Biostrings::readAAStringSet(path2fasta)
orig_names <- sapply(names(fa), function(y) {
    strsplit(y, headerSep)[[1]][1]
})

cat(head(orig_names), sep = "\n")

# -- parse the peptide headers
names(fa) <- gsub(headerStripText, "", sapply(names(fa), function(y) {
    strsplit(y, headerSep)[[1]][headerEntryIndex]
}))

# -- remove all special characters
names(fa) <- gsub("[^a-zA-Z0-9_.-]", convertSpecialCharacters, names(fa))

df <- data.frame(GeneSymbol = names(fa), Header = unname(orig_names), ProteinLength = width(fa))
head(df)

# -- ensure no duplicated headers
if(any(duplicated(names(fa)))){
    df <- df[order(-width(fa)),]
    fa <- fa[order(-width(fa)),]
    df <- df[!duplicated(names(fa)),]
    fa <- fa[!duplicated(names(fa)),]
}

cat(head(names(fa)), sep = "\n")
head(df)

write.table(df, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
