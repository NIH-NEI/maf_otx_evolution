from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os
import argparse

organisms = (
    pd.read_csv("configs/species173_table.tsv", sep = "\t")
)
print(organisms)

vertebrates = organisms.query("IsVertebrate == 'Yes'").OrganismShortName.tolist()
invertebrates = organisms.query("IsVertebrate == 'No'").OrganismShortName.tolist()

print(vertebrates)
print(invertebrates)

def create_peptide_files(group, outdir):
    org = group.name
    outfile = f"{outdir}/{org}.faa"
    records = []
    for _, row in group.iterrows():
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

def tsv_to_multi_fasta(infile, outdir, species_group):
    df = (
        pd.read_csv(infile, sep = "\t")
        .assign(organism = lambda tdf: tdf.ProteinID.str.split("__").str[0])

    )
    if species_group == "invertebrate":
        df = df.query("organism in @invertebrates")
    elif species_group == "vertebrate":
        df = df.query("organism in @vertebrates")
    print(df)

    os.makedirs(outdir, exist_ok=True)

    df.groupby("organism").apply(lambda grp: create_peptide_files(grp, outdir))


def main():
    parser = argparse.ArgumentParser(description="Convert protein TSV to FASTA files.")
    parser.add_argument("protein_tsv", help="Input protein TSV file (with ProteinID and Sequence columns)")
    parser.add_argument("output_dir", help="Directory to save FASTA files")
    parser.add_argument("species_group", nargs = "?", default="ignore", help="Directory to save FASTA files")
    args = parser.parse_args()
    print(args)

    tsv_to_multi_fasta(args.protein_tsv, args.output_dir, args.species_group)

if __name__ == "__main__":
    main()
