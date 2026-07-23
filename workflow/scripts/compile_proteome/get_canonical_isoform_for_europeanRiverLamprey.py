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

def parse_attr(attr_field):
    attrs = {}
    print(attr_field)
    for item in attr_field.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k] = v.strip('"')
        else:
            # Some GFFs may have standalone flags; ignore
            pass
    return attrs

# Accept the common RefSeq protein prefixes
_PID_PREFIX = r'(?:NP|XP|YP|WP|AP)'
_PID_CORE   = rf'{_PID_PREFIX}_[0-9]+(?:\.[0-9]+)?'
RE_PROTEIN_ID_FIELD = re.compile(rf'protein_id=({_PID_CORE})')
RE_PID_TOKEN        = re.compile(rf'\b({_PID_CORE})\b')
RE_PROT_SEGMENT     = re.compile(rf'_prot_({_PID_CORE})')

def normalize_protein_id_from_seqrecord(rec):
    """
    Return the canonical protein_id (e.g., NP_123456.1) from an NCBI protein FASTA SeqRecord.
    Tries, in order:
      1) [protein_id=NP_...] in description
      2) any NP_/XP_/YP_/WP_/AP_ token in description
      3) any such token in rec.id
      4) _prot_NP_... pattern in rec.id
    If none found, returns rec.id as a fallback (but that likely won't match GFF protein_id).
    """
    desc = rec.description or ""
    # this is a hack as the protein file have only one protein per gene
    rid  = rec.id or ""

    print(rec)
    m = RE_PROTEIN_ID_FIELD.search(desc)
    if m:
        return m.group(1)

    m = RE_PID_TOKEN.search(desc)
    if m:
        return m.group(1)

    m = RE_PID_TOKEN.search(rid)
    if m:
        return m.group(1)

    m = RE_PROT_SEGMENT.search(rid)
    if m:
        return m.group(1)

    print(rid)
    return rid

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
    mrna_tags = defaultdict(set)          # mRNA_id -> set of 'tags'
    mrna_parent_gene = {}                 # mRNA_id -> gene_id
    gene_symbols = {}                     # gene_id -> gene symbol (if any)
    gene_dbxref = {}                     # gene_id -> gene symbol (if any)
    mrna_to_protein = defaultdict(set)    # mRNA_id -> set(protein_ids)
    protein_to_mrna = {}                  # protein_id -> mRNA_id

    # Parse GTF(s) - only CDS features
    for gtf in gtf_paths:
        #with gzip.open(gtf, "rt") as fh:
        with open(gtf, "r") as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                print(cols)
                if len(cols) < 9:
                    continue
                seqid, source, ftype, start, end, score, strand, phase, attrs_raw = cols
                attrs = parse_attr(attrs_raw)
                print(attrs)

                if ftype == "gene":
                    gid = attrs.get("ID")
                    if gid:
                        # Best-effort store a symbol for nicer headers
                        sym = attrs.get("Name") or gid
                        gene_symbols[gid] = sym
                        dbxref = attrs.get("Dbxref") or gid
                        gene_dbxref[gid] = dbxref

                if ftype in ("mRNA", "transcript"):
                    mid = attrs.get("ID")
                    print(ftype, mid)
                    if not mid:
                        continue
                    parent_gene = attrs.get("Parent") or attrs.get("gene_id") or attrs.get("locus_tag")
                    # Some files set Parent=gene-XYZ; prefer that exact value if present.
                    if parent_gene:
                        # If multiple parents are given, use the first
                        parent_gene = parent_gene.split(",")[0]
                    else:
                        # Last resort: use gene symbol as ID
                        parent_gene = attrs.get("gene") or attrs.get("locus_tag") or attrs.get("Name") or mid
                    mrna_parent_gene[mid] = parent_gene
                    print("myself:", mid, "parent:", parent_gene)

                    # parse tags (e.g., tag=MANE Select, tag=RefSeq Select)
                    tag_field = attrs.get("tag", "")
                    if tag_field:
                        for t in tag_field.split(","):
                            mrna_tags[mid].add(t.strip())

                if ftype == "CDS":
                    parent = attrs.get("Parent")
                    pid = attrs.get("protein_id")
                    if parent and pid:
                        # Parent can be a list; take each mRNA ID
                        for mid in parent.split(","):
                            mid = mid.strip()
                            mrna_to_protein[mid].add(pid)
                            protein_to_mrna[pid] = mid

    print("Processed GTF with number of genes:", len(mrna_parent_gene))

    # Read protein FASTAs (index by normalized protein_id)
    protein_records = {}
    unmatched_examples = []
    for faa in args.protein_faa:
        p = Path(faa)
        if not p.exists():
            continue
        with open(faa, "rt") as handle:
            for rec in SeqIO.parse(handle, "fasta"):
                pid = normalize_protein_id_from_seqrecord(rec)
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

    print("Read fasta with number of proteins:", len(protein_records))
    print(list(protein_records.keys())[:10])
    print(list(gene_symbols.items())[:10])

    # Optional: quick debug note if nothing matched canonically-shaped IDs
    if not protein_records and unmatched_examples:
        sys.stderr.write("Warning: no normalized protein_ids found; examples:\n")
        for ex in unmatched_examples:
            sys.stderr.write(f"  {ex}\n")


    # Compute a mapping: gene -> candidate proteins (by mRNA)
    gene_to_mrnas = defaultdict(set)
    for mid, gid in mrna_parent_gene.items():
        gene_to_mrnas[gid].add(mid)

    # Helper to get length of a protein_id
    def prot_len(pid):
        rec = protein_records.get(pid)
        return len(rec.seq) if rec else 0

    # Select per gene
    selected_pids = []
    for gid, mrnas in gene_to_mrnas.items():
        print("gid", gid, "mrna", mrnas)
        all_pids = [p for p in mrnas if prot_len(p) > 0]
        print("all_pids", all_pids)
        if all_pids:
            best = max(all_pids, key=prot_len)
            candidate_pids = [best]
        print("candidate_pids", candidate_pids)

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
        sys.exit(1)
    else:
        print("Num of isoforms selected: ", len(unique.keys()))

    # Write FASTA
    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    info = []
    src = "AH2P"
    acc = "Custom"
    with open(outp, "w") as out_fh:
        for (gid, pid) in sorted(unique.keys()):
            rec = protein_records.get(pid)
            if not rec:
                continue
            sym = gene_symbols.get(gid, gid)
            ref = gid  # No separate Dbxref; use gene_id
            org = args.organism or ""
            # Header example:
            # >pid|org=species|gene=SYMBOL|gid=gene-123|dbxref=gene-123|assembly=GCF_...|source=TransDecoder
            #header = f"{pid}|org={org}|gene={sym}|gid={gid}|dbxref={ref}|assembly=Custom|source=funanotate"
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

