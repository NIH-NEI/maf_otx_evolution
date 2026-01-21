#!/usr/bin/env python3
import argparse
import sys
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
import gzip
import json
import pandas as pd
from typing import Dict, List, Any

def parse_attributes(attr_str: str) -> Dict[str, str]:
    """
    Parse the attributes column from GFF/GTF.
    Handles both GTF (key "value";) and GFF3 (key=value;) formats.
    """
    if pd.isna(attr_str) or not attr_str:
        return {}
    
    attrs = {}
    parts = attr_str.split(';')
    for part in parts:
        part = part.strip()
        if not part:
            continue
        if '=' in part:
            if ' ' in part:
                # GTF style: key "value"
                key_value = part.split(' ', 1)
                key = key_value[0].strip()
                value = key_value[1].strip().strip('"')
            else:
                # GFF3 style: key=value
                key, value = part.split('=', 1)
                key = key.strip()
                value = value.strip().strip('"')
            attrs[key] = value
    return attrs

def parse_gff_gtf(file_path: str) -> tuple[Dict[str, List[str]], Dict[str, str], Dict[str, Dict[str, Any]], Dict[str, Dict[str, Any]]]:
    """
    Parse a GTF or GFF file to extract gene-transcript mappings and metadata.
    
    Returns:
    - gene_to_transcripts: Dict[gene_id, List[transcript_id]]
    - transcript_to_gene: Dict[transcript_id, gene_id]
    - gene_metadata: Dict[gene_id, Dict[str, Any]] (seqname, start, end, strand, gene_name, gene_type)
    - transcript_metadata: Dict[transcript_id, Dict[str, Any]] (seqname, start, end, strand, transcript_name, transcript_type, gene_id)
    """
    # Read the file, skipping header lines starting with '#'
    df = pd.read_csv(
        file_path,
        sep='\t',
        comment='#',
        header=None,
        names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes'],
        dtype={'start': 'int64', 'end': 'int64', 'score': 'float64', 'phase': str}
    )
    
    # Parse attributes into dicts
    df['attributes_dict'] = df['attributes'].apply(parse_attributes)
    
    # Filter for gene and transcript features
    # Handle 'transcript' or 'mRNA' for transcripts
    genes_df = df[df['feature'] == 'gene'].copy()
    transcripts_df = df[df['feature'].isin(['transcript', 'mRNA'])].copy()
    
    # Gene metadata
    gene_metadata = {}
    for _, row in genes_df.iterrows():
        attrs = row['attributes_dict']
        gene_id = attrs.get('gene_id') or attrs.get('ID')
        if not gene_id:
            continue
        gene_metadata[gene_id] = {
            'seqname': row['seqname'],
            'start': int(row['start']),
            'end': int(row['end']),
            'strand': row['strand'],
            'gene_name': attrs.get('gene_name') or attrs.get('Name') or attrs.get('gene_biotype'),
            'gene_type': attrs.get('gene_type') or attrs.get('biotype') or attrs.get('gene_biotype'),
            'source': row['source']
        }
    
    # Transcript metadata and mappings
    gene_to_transcripts = defaultdict(list)
    transcript_metadata = {}
    for _, row in transcripts_df.iterrows():
        attrs = row['attributes_dict']
        transcript_id = attrs.get('transcript_id') or attrs.get('ID')
        if not transcript_id:
            continue
        gene_id = attrs.get('gene_id') or attrs.get('Parent')
        if not gene_id:
            continue
        
        transcript_metadata[transcript_id] = {
            'seqname': row['seqname'],
            'start': int(row['start']),
            'end': int(row['end']),
            'strand': row['strand'],
            'transcript_name': attrs.get('transcript_name') or attrs.get('Name') or attrs.get('transcript_biotype'),
            'transcript_type': attrs.get('transcript_type') or attrs.get('biotype') or attrs.get('transcript_biotype'),
            'gene_id': gene_id,
            'source': row['source']
        }
        gene_to_transcripts[gene_id].append(transcript_id)
    
    # Bidirectional mapping
    transcript_to_gene = {tid: meta['gene_id'] for tid, meta in transcript_metadata.items()}
    
    return dict(gene_to_transcripts), transcript_to_gene, gene_metadata, transcript_metadata

def main():
    if len(sys.argv) != 2:
        print("Usage: python gff_parser.py <path_to_gtf_or_gff_file>", file=sys.stderr)
        sys.exit(1)
    
    file_path = sys.argv[1]
    
    try:
        gene_to_trans, trans_to_gene, gene_meta, trans_meta = parse_gff_gtf(file_path)
        
        # Output to JSON files
        base_name = file_path.rsplit('.', 1)[0]
        with open(f'{base_name}_gene_to_transcripts.json', 'w') as f:
            json.dump(gene_to_trans, f, indent=2)
        with open(f'{base_name}_transcript_to_gene.json', 'w') as f:
            json.dump(trans_to_gene, f, indent=2)
        with open(f'{base_name}_gene_metadata.json', 'w') as f:
            json.dump(gene_meta, f, indent=2)
        with open(f'{base_name}_transcript_metadata.json', 'w') as f:
            json.dump(trans_meta, f, indent=2)
        
        print(f"Successfully parsed {file_path}. Output files:")
        print(f"  - {base_name}_gene_to_transcripts.json")
        print(f"  - {base_name}_transcript_to_gene.json")
        print(f"  - {base_name}_gene_metadata.json")
        print(f"  - {base_name}_transcript_metadata.json")
        print(f"\nSummary: {len(gene_to_trans)} genes, {len(trans_meta)} transcripts.")
    
    except Exception as e:
        print(f"Error parsing {file_path}: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()




def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", required=True, help="Species name to be added to protein header")
    ap.add_argument("--gtf", nargs="+", required=True, help="One or more GTF files (wildcards expanded by shell)")
    ap.add_argument("--protein-faa", nargs="+", required=True, help="One or more peptide FASTA files")
    ap.add_argument("--out", required=True, help="Output canonical peptide FASTA")
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
        gene_to_trans, trans_to_gene, gene_meta, trans_meta = parse_gff_gtf(gtf)
        print(gene_to_trans)


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
            header = f"{tid}|org={org}|gene={sym}|gid={gid}|dbxref={ref}|assembly=Custom|source=TransDecoder"
            out_fh.write(f">{header}\n")
            # Wrap sequence at 60 chars
            seq = str(rec.seq)
            for i in range(0, len(seq), 60):
                out_fh.write(seq[i:i+60] + "\n")

    # Report
    print(f"Wrote {outp}")

if __name__ == "__main__":
    main()

