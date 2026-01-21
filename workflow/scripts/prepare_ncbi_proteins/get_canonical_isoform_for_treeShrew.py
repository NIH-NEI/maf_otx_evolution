#!/usr/bin/env python3
import argparse
import sys
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
import gzip

"""
Logic:
1) Read all GTF files to build:
   - gene_id -> set of transcript_ids (from CDS features)
   - gene_id -> gene_name (if present)
2) Read protein FASTA(s) to get protein sequences and lengths, indexed by transcript_id.
3) Per gene:
   - Pick the longest protein across transcripts of the gene.
4) Write selected protein FASTA entries, de-duplicated, with informative headers.

Notes:
- Custom TransDecoder GTF encodes:
  CDS lines with attributes: transcript_id "TSDBTID..."; gene_id "TSDBGID..."; gene_name "..."; biotype "protein coding"
- Each transcript corresponds to one protein, named by transcript_id in the FASTA header.
- No canonical tags; fallback to longest protein per gene.
- Only processes CDS features to identify coding transcripts.
"""

def parse_gtf_attributes(attr_str):
    """
    Parse GTF attributes string into a dict.
    Handles format: key "value"; key2 "value2";
    Strips quotes from values.
    """
    attrs = {}
    # Split on '; ' (semicolon followed by space)
    parts = [p.strip() for p in attr_str.split(';') if p.strip()]
    for part in parts:
        if ' ' in part:
            key, value = part.split(' ', 1)
            # Strip surrounding quotes if present
            if value.startswith('"') and value.endswith('"'):
                value = value[1:-1]
            attrs[key] = value
    return attrs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True, help="Species name to be added to protein header")
    ap.add_argument("--gtf", nargs="+", required=True, help="One or more GTF files (wildcards expanded by shell)")
    ap.add_argument("--protein-faa", nargs="+", required=True, help="One or more peptide FASTA files")
    ap.add_argument("--out", required=True, help="Output canonical peptide FASTA")
    ap.add_argument("--info", required=True, help="Output info about canonical peptide FASTA")
    args = ap.parse_args()

    # Gather GTF paths (filter non-existent)
    gtf_paths = [Path(p) for p in args.gtf if Path(p).exists()]
    if not gtf_paths:
        sys.exit("No GTF files found. Aborting.")

    # Maps
    gene_to_transcripts = defaultdict(set)  # gene_id -> set(transcript_ids)
    gene_symbols = {}                      # gene_id -> gene_name (if any)

    # Parse GTF(s) - only CDS features
    for gtf in gtf_paths:
        with gzip.open(gtf, "rt") as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    continue
                seqid, source, ftype, start, end, score, strand, phase, attrs_raw = cols
                if ftype != "CDS":
                    continue
                attrs = parse_gtf_attributes(attrs_raw)
                gene_id = attrs.get("gene_id")
                transcript_id = attrs.get("transcript_id")
                biotype = attrs.get("biotype")
                if gene_id and transcript_id and biotype == "protein coding":
                    gene_to_transcripts[gene_id].add(transcript_id)
                    # Store gene name if present
                    gene_name = attrs.get("gene_name", gene_id)
                    gene_symbols[gene_id] = gene_name

    print("Processed GTF with number of genes:", len(gene_to_transcripts))

    # Read peptide FASTAs (index by transcript_id)
    protein_records = {}
    for faa in args.protein_faa:
        p = Path(faa)
        if not p.exists():
            continue
        with gzip.open(faa, "rt") as handle:
            for rec in SeqIO.parse(handle, "fasta"):
                tid = rec.id
                # Keep the longest if duplicates (rare but safe)
                if tid in protein_records:
                    if len(rec.seq) > len(protein_records[tid].seq):
                        protein_records[tid] = rec
                else:
                    protein_records[tid] = rec

    # Helper to get length of a transcript_id (protein)
    def prot_len(tid):
        rec = protein_records.get(tid)
        return len(rec.seq) if rec else 0

    # Select per gene: longest protein
    selected_tids = []
    for gid, transcripts in gene_to_transcripts.items():
        all_tids = [t for t in transcripts if prot_len(t) > 0]
        if all_tids:
            best_tid = max(all_tids, key=prot_len)
            selected_tids.append((gid, best_tid))

    # Deduplicate by (gene, tid) - though unlikely
    unique = {}
    for gid, tid in selected_tids:
        unique[(gid, tid)] = 1
    if not unique:
        print("Error: No protein isoforms was selected!")
        sys.exit(1)
    else:
        print("Num of isoforms selected: ", len(unique.keys()))

    # Write FASTA
    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    info = []
    src = "TransDecoder"
    acc = "Custom"
    with open(outp, "w") as out_fh:
        for (gid, tid) in sorted(unique.keys()):
            rec = protein_records.get(tid)
            if not rec:
                continue
            sym = gene_symbols.get(gid, gid)
            ref = gid  # No separate Dbxref; use gene_id
            org = args.organism or ""
            # Header example:
            # >tid|org=species|gene=SYMBOL|gid=gene-123|dbxref=gene-123|assembly=GCF_...|source=TransDecoder
            #header = f"{tid}|org={org}|gene={sym}|gid={gid}|dbxref={ref}|assembly=Custom|source=TransDecoder"
            header = f"{org}__{sym}__{tid}"
            info.append([tid, org, sym, gid, ref, acc, src])
            out_fh.write(f">{header}\n")
            # Wrap sequence at 60 chars
            seq = str(rec.seq)
            for i in range(0, len(seq), 60):
                out_fh.write(seq[i:i+60] + "\n")
    infop = Path(args.info)
    infop.parent.mkdir(parents=True, exist_ok=True)
    import pandas as pd
    df = pd.DataFrame(info, columns = [
        "ProteinID",
        "OrganismShortName",
        "GeneSymbol",
        "GeneID",
        "dbxref",
        "assembly",
        "source",
    ])
    print(df)
    df.to_csv(infop, sep = "\t", index = False)

    # Report
    print(f"Wrote {outp}")

if __name__ == "__main__":
    main()

