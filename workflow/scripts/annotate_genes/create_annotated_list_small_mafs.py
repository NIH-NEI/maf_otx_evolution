from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

orgs = (
    pd.read_csv("configs/species173_table.tsv", sep = "\t")
    [["OrganismShortName", "OrganismID", "OrganismColor", "IsVertebrate"]]
    .merge(
        (
            pd.read_csv("configs/species173_annotations.tsv", sep = "\t")
            [["OrganismShortName", "AnnotationSource"]]
        ),
        how  ="left",
    )
)
print(orgs)

def add_protein_seqs(group: pd.DataFrame) -> pd.DataFrame:
    # group.name is the value of the "pep" column for this group
    pep_file = group.name

    print(f"Processing {pep_file} ...")

    required = group.ProteinID.tolist()

    # Read all sequences from this fasta once
    seq_dict = {}
    with open(pep_file, "r") as spe:
        for record in SeqIO.parse(spe, "fasta"):
            if record.id in required:
                seq_dict[record.id] = str(record.seq)

    if len(required) != len(seq_dict):
        print("Error: Not all sequences are found!")
        print("Required: ", required)
        exit(-1)

    # Map ProteinID to sequences
    group["ProteinSeq"] = group["ProteinID"].map(seq_dict.get)

    return group


# Build the base table first (without ProteinSeq)
base = (
    pd.read_csv("exports/curated_genes/draft_annotated_MAFS_longformat.tsv", sep="\t")
    .merge(orgs, how="left")
    #.assign(ProteinID=lambda tdf: tdf.OrganismShortName + "__" + tdf.ProteinID)
    .assign(pep=lambda tdf: "scratch/canonical_peptides/" + tdf.OrganismShortName + ".faa")
)

# Group by pep and process each FASTA file once
annotated_mafs = (
    base
    .groupby("pep", group_keys=False)
    .apply(add_protein_seqs)
)

print(annotated_mafs)

def write_df_to_fasta(df, outfile):
    records = []
    for _, row in df.iterrows():
        if pd.isna(row["ProteinSeq"]):
            continue  # skip missing sequences

        rec = SeqRecord(
            Seq(row["ProteinSeq"]),
            id=row["ProteinID"],
            description=""   # no trailing description
        )
        records.append(rec)

    SeqIO.write(records, outfile, "fasta")
    print(f"Wrote {len(records)} sequences to {outfile}")


# Example usage:
write_df_to_fasta(
    annotated_mafs.query("IsVertebrate == 'No'"),
    "scratch/blast_query_proteins/query_invertebrate_MAFS.faa"
)

write_df_to_fasta(
    annotated_mafs.query("IsVertebrate == 'Yes'").query("GeneGroup == 'MAFF'"),
    "scratch/blast_query_proteins/query_vertebrate_MAFF.faa"
)

write_df_to_fasta(
    annotated_mafs.query("IsVertebrate == 'Yes'").query("GeneGroup == 'MAFG'"),
    "scratch/blast_query_proteins/query_vertebrate_MAFG.faa"
)

write_df_to_fasta(
    annotated_mafs.query("IsVertebrate == 'Yes'").query("GeneGroup == 'MAFK'"),
    "scratch/blast_query_proteins/query_vertebrate_MAFK.faa"
)

