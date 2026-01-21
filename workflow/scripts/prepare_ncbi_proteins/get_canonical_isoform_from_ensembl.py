#!/usr/bin/env python3
import argparse
import sys
from collections import defaultdict, namedtuple
from pathlib import Path
from Bio import SeqIO
import gzip
import re

"""
Logic:
1) Read all GFF3 files to build:
   - transcript -> tags (look for 'tag=Ensembl_canonical' or 'tag=MANE_Select')
   - transcript -> parent gene (gene ID)
   - transcript -> protein_id(s)
   - protein_id -> (transcript, gene)
2) Read protein FASTA(s) to get protein sequences and lengths.
3) Per gene:
   - If any transcript has a canonical tag, take its protein(s) (typically one).
   - Else pick the longest protein across transcripts of the gene.
4) Write selected protein FASTA entries, de-duplicated, with informative headers.

Notes:
- Ensembl GFF3 typically encodes:
  transcript lines with attributes: ID=ENST..., Parent=ENSG..., tag=Ensembl_canonical|MANE_Select (when present)
  CDS lines with attributes: Parent=<transcript-ID>, protein_id=ENSP_...
- Genes can be addressed by ID=ENSG...; symbols via Name=.
- Updated from NCBI RefSeq to Ensembl/GENCODE.
"""

def parse_attr(attr_field):
    attrs = {}
    for item in attr_field.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k] = v
        else:
            # Some GFFs may have standalone flags; ignore
            pass
    return attrs

# Accept Ensembl protein prefixes
_PID_PREFIX = r'(?:ENSP)'
_PID_CORE   = rf'{_PID_PREFIX}[0-9]+(?:\.[0-9]+)?'
RE_PROTEIN_ID_FIELD = re.compile(rf'protein_id=({_PID_CORE})')
RE_PID_TOKEN        = re.compile(rf'\b({_PID_CORE})\b')
# No need for _prot_ segment in Ensembl

def normalize_protein_id_from_seqrecord(rec):
    """
    Return the canonical protein_id (e.g., ENSP00000369572.7) from an Ensembl protein FASTA SeqRecord.
    Tries, in order:
      1) [protein_id=ENSP...] in description
      2) any ENSP token in description
      3) any such token in rec.id
    If none found, returns rec.id as a fallback.
    """
    desc = rec.description or ""
    rid  = rec.id or ""

    m = RE_PROTEIN_ID_FIELD.search(desc)
    if m:
        return m.group(1)

    m = RE_PID_TOKEN.search(desc)
    if m:
        return m.group(1)

    m = RE_PID_TOKEN.search(rid)
    if m:
        return m.group(1)

    return rid

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True, help="Species name to be added to protein header")
    ap.add_argument("--accession", required=True, help="Assembly accession for labeling")
    ap.add_argument("--gff", nargs="+", required=True, help="One or more GFF3 files (wildcards expanded by shell)")
    ap.add_argument("--protein-faa", nargs="+", required=True, help="One or more protein FASTA files")
    ap.add_argument("--out", required=True, help="Output canonical peptide FASTA")
    ap.add_argument("--info", required=True, help="Output info about canonical peptide FASTA")
    args = ap.parse_args()

    # Gather GFF paths (filter non-existent)
    gff_paths = [Path(p) for p in args.gff if Path(p).exists()]
    if not gff_paths:
        sys.exit("No GFF3 files found. Aborting.")

    # Maps
    transcript_tags = defaultdict(set)    # transcript_id -> set of 'tags'
    transcript_parent_gene = {}           # transcript_id -> gene_id
    gene_symbols = {}                     # gene_id -> gene symbol (if any)
    gene_dbxref = {}                      # gene_id -> dbxref (if any)
    transcript_to_protein = defaultdict(set)  # transcript_id -> set(protein_ids)
    protein_to_transcript = {}            # protein_id -> transcript_id

    # Parse GFF3(s)
    for gff in gff_paths:
        with gzip.open(gff, "rt") as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    continue
                seqid, source, ftype, start, end, score, strand, phase, attrs_raw = cols
                attrs = parse_attr(attrs_raw)

                if ftype == "gene":
                    gid = attrs.get("ID") or attrs.get("gene_id")
                    if gid:
                        # Best-effort store a symbol for nicer headers
                        sym = attrs.get("Name") or attrs.get("gene_name") or attrs.get("gene") or gid
                        gene_symbols[gid] = sym
                        dbxref = attrs.get("Dbxref") or gid
                        gene_dbxref[gid] = dbxref

                if ftype in ("mRNA", "transcript"):
                    tid = attrs.get("ID")
                    if not tid:
                        continue
                    parent_gene = attrs.get("Parent") or attrs.get("gene_id") or attrs.get("gene")
                    # Some files set Parent=ENSG-XYZ; prefer that exact value if present.
                    if parent_gene:
                        # If multiple parents are given, use the first
                        parent_gene = parent_gene.split(",")[0]
                    else:
                        # Last resort: use gene symbol as ID
                        parent_gene = attrs.get("gene") or attrs.get("gene_name") or attrs.get("Name") or tid
                    transcript_parent_gene[tid] = parent_gene

                    # parse tags (e.g., tag=Ensembl_canonical, tag=MANE_Select)
                    tag_field = attrs.get("tag", "")
                    if tag_field:
                        for t in tag_field.split(","):
                            transcript_tags[tid].add(t.strip())

                if ftype == "CDS":
                    parent = attrs.get("Parent")
                    pid = attrs.get("protein_id")
                    if parent and pid:
                        # Parent can be a list; take each transcript ID
                        for tid in parent.split(","):
                            tid = tid.strip()
                            transcript_to_protein[tid].add(pid)
                            protein_to_transcript[pid] = tid


    # Read protein FASTAs (index by normalized protein_id)
    protein_records = {}
    unmatched_examples = []

    for faa in args.protein_faa:
        p = Path(faa)
        if not p.exists():
            continue
        with gzip.open(faa, "rt") as handle:
            for rec in SeqIO.parse(handle, "fasta"):
                pid = normalize_protein_id_from_seqrecord(rec)
                # in ensemble gff3 the protein version is not kept
                # so we need to keep only the ensemble id, not version
                pid = pid.split(".")[0]
                if pid == rec.id:
                    # keep a couple examples to help debug if needed
                    if len(unmatched_examples) < 3:
                        unmatched_examples.append(rec.description)
                # Keep the longest if duplicates (rare but safe)
                if pid in protein_records:
                    if len(rec.seq) > len(protein_records[pid].seq):
                        protein_records[pid] = rec
                else:
                    protein_records[pid] = rec

    # Optional: quick debug note if nothing matched canonically-shaped IDs
    if not protein_records and unmatched_examples:
        sys.stderr.write("Warning: no normalized protein_ids found; examples:\n")
        for ex in unmatched_examples:
            sys.stderr.write(f"  {ex}\n")



    # Compute a mapping: gene -> candidate transcripts
    gene_to_transcripts = defaultdict(set)
    for tid, gid in transcript_parent_gene.items():
        gene_to_transcripts[gid].add(tid)

    # Helper to get length of a protein_id
    def prot_len(pid):
        rec = protein_records.get(pid)
        return len(rec.seq) if rec else 0

    # Select per gene
    selected_pids = []
    for gid, transcripts in gene_to_transcripts.items():
        # 1) Canonical preference
        canonical_transcripts = [t for t in transcripts if any(tag in ("Ensembl_canonical", "MANE_Select") for tag in transcript_tags.get(t, []))]
        candidate_pids = []
        if canonical_transcripts:
            for t in canonical_transcripts:
                candidate_pids.extend(list(transcript_to_protein.get(t, [])))
        else:
            # 2) Fallback: longest protein across all transcripts in this gene
            all_pids = []
            for t in transcripts:
                all_pids.extend(list(transcript_to_protein.get(t, [])))
            if all_pids:
                best = max(all_pids, key=prot_len)
                candidate_pids = [best]

        # Record
        for pid in candidate_pids:
            if pid in protein_records:
                selected_pids.append((gid, pid))


    # Deduplicate by (gene, pid)
    unique = {}
    for gid, pid in selected_pids:
        unique[(gid, pid)] = 1


    if not unique:
        print("Error: No protein isoforms was selected!")
        quit()
    else:
        print("Num of isoforms selected: ", len(unique.keys()))

    # Write FASTA
    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    info = []
    src = "Ensembl"
    with open(outp, "w") as out_fh:
        for (gid, pid) in sorted(unique.keys()):
            rec = protein_records.get(pid)
            if not rec:
                continue
            sym = gene_symbols.get(gid, gid)
            ref = gene_dbxref.get(gid, gid)
            acc = args.accession or ""
            org = args.organism or ""
            # Header example:
            # >pid|org=human|gene=SYMBOL|gid=ENSG...|dbxref=ENSG...|assembly=GRC...|source=Ensembl
            #header = f"{pid}|org={org}|gene={sym}|gid={gid}|dbxref={ref}|assembly={acc}|source=Ensembl"
            header = f"{org}__{sym}__{pid}"
            info.append([pid, org, sym, gid, ref, acc, src])
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
