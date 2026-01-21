from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os

import subprocess

def shell(cmd):
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)

infile = "scratch/for_mini_orthofinder/cow_combined_denovo_annotated.tsv"
top_dir = os.path.realpath("scratch/for_mini_orthofinder/cow")
outdir = f"{top_dir}/peptides"
print(outdir)

os.makedirs(outdir, exist_ok=True)

def create_peptide_files(group):
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

df = (
    pd.read_csv(infile, sep = "\t")
    .assign(organism = lambda tdf: tdf.ProteinID.str.split("__").str[0])
)
print(df)
print(df.organism.unique())

df.groupby("organism").apply(create_peptide_files)

shell(
    "module load OrthoFinder;"
    f"cd {top_dir};"
    f"orthofinder -f {outdir};"
)
