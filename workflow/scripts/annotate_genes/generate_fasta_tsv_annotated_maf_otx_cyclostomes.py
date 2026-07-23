from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

orgs = (
    pd.read_csv("configs/species1XX_table.tsv", sep = "\t")
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

    from pathlib import Path

    file_path = Path(pep_file)

    if not file_path.exists():
        raise ValueError("File does not exist:", pep_file)

    print(f"Processing {pep_file} ...")

    required = group.ProteinID.tolist()

    # Read all sequences from this fasta once
    seq_dict = {}
    len_dict = {}
    with open(pep_file, "r") as spe:
        for record in SeqIO.parse(spe, "fasta"):
            record.id = record.id.replace("__gene:ENS", "__gene_ENS")
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
    pd.read_csv("exports/curated_genes/draft_annotated_maf_otx_longformat_cyclostomes.tsv", sep="\t")
    .assign(Remark = lambda tdf: "Annotated" + tdf.GeneGroup)
    .merge(orgs, how="left")
    .assign(pep=lambda tdf: "scratch/canonical_peptides_unique/" + tdf.OrganismShortName + ".faa")
)
print(base)

# Group by pep and process each FASTA file once
annotated_maf_otx = (
    base
    .groupby("pep", group_keys=False)
    .apply(add_protein_seqs)
)

print(annotated_maf_otx)
annotated_maf_otx.to_csv("scratch/grand_list/extra_cyclostomes.tsv", sep ="\t", index = False)
quit()

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
    annotated_maf_otx,
    "scratch/annotated_proteins/annotated_maf_otx.faa"
)

write_df_to_fasta(
    annotated_maf_otx.query("IsVertebrate == 'Yes'"),
    "scratch/annotated_proteins/annotated_vertebrate_maf_otx.faa"
)

write_df_to_fasta(
    annotated_maf_otx.query("IsVertebrate == 'No'"),
    "scratch/annotated_proteins/annotated_invertebrate_maf_otx.faa"
)
