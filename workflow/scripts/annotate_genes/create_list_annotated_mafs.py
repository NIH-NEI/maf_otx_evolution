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
    len_dict = {}
    with open(pep_file, "r") as spe:
        for record in SeqIO.parse(spe, "fasta"):
            if record.id in required:
                seq_dict[record.id] = str(record.seq)
                len_dict[record.id] = len(record.seq)

    if len(required) != len(seq_dict):
        print("Error: Not all sequences are found!")
        print("Required: ", required)
        exit(-1)

    # Map ProteinID to sequences
    group["ProteinSeq"] = group["ProteinID"].map(seq_dict.get)
    group["ProteinLen"] = group["ProteinID"].map(len_dict.get)

    return group


# Build the base table first (without ProteinSeq)
base = (
    pd.concat([
      (
          pd.read_csv("exports/curated_genes/draft_annotated_MAFL_longformat.tsv", sep="\t")
          .assign(Remark = "AnnotatedLargeMAF")
      ),
      (
          pd.read_csv("exports/curated_genes/draft_annotated_MAFS_longformat.tsv", sep="\t")
          .assign(Remark = "AnnotatedSmallMAF")
      ),
    ])
    .merge(orgs, how="left")
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

    tsvfile = outfile.replace(".faa", ".tsv")
    df[["ProteinID", "ProteinLen", "Remark", "ProteinSeq"]].to_csv(tsvfile, sep = "\t", index=False)

write_df_to_fasta(
    annotated_mafs.query("IsVertebrate == 'Yes'"),
    "scratch/annotated_proteins/query_vertebrate_ALLMAF.faa"
)


write_df_to_fasta(
    annotated_mafs.query("IsVertebrate == 'No'"),
    "scratch/annotated_proteins/query_invertebrate_ALLMAF.faa"
)

