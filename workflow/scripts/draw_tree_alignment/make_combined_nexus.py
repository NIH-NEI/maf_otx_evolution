#!/usr/bin/env python3

"""
Combine multiple FASTA alignments and Newick trees into one NEXUS file.

Input index file format:

    name,type,path
    gene1,alignment,gene1.fasta
    gene2,alignment,gene2.fasta
    gene1_tree,tree,gene1.tree
    gene2_tree,tree,gene2.tree

The script:
  1. Reads all FASTA alignments.
  2. Verifies each alignment is valid.
  3. Builds a union of all taxa.
  4. Concatenates alignments in the order listed.
  5. Fills missing taxa with ? characters for that partition.
  6. Adds CHARSET partition definitions.
  7. Adds all trees into a TREES block.
  8. Writes one NEXUS file.

Requires:
    pip install biopython
"""

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import AlignIO


VALID_ALIGNMENT_TYPES = {"alignment", "align", "aln", "fasta", "fa"}
VALID_TREE_TYPES = {"tree", "newick", "nwk", "tre"}

top_dir = Path(".")


def quote_nexus_name(name: str) -> str:
    """
    Quote a NEXUS name if it contains spaces or special characters.
    Single quotes inside names are escaped by doubling them.
    """
    if re.fullmatch(r"[A-Za-z0-9_.:-]+", name):
        return name
    return "'" + name.replace("'", "''") + "'"


def clean_partition_name(name: str) -> str:
    """
    Make a safe NEXUS partition/tree identifier.
    """
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
    if not cleaned:
        cleaned = "partition"
    if re.match(r"^[0-9]", cleaned):
        cleaned = "p_" + cleaned
    return cleaned


def read_index(index_file: Path) -> List[Dict[str, str]]:
    """
    Read CSV or TSV index file with columns:
        name,type,path

    Delimiter is inferred from the first line.
    """
    text = index_file.read_text().splitlines()

    if not text:
        raise ValueError(f"Index file is empty: {index_file}")

    first_line = text[0]
    delimiter = "\t" if "\t" in first_line else ","

    with index_file.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)

        required = {"name", "type", "path"}
        if reader.fieldnames is None:
            raise ValueError("Index file has no header row.")

        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(
                f"Index file is missing required column(s): {', '.join(sorted(missing))}"
            )

        rows = []
        for i, row in enumerate(reader, start=2):
            name = row["name"].strip()
            file_type = row["type"].strip().lower()
            path = row["path"].strip()

            if not name or not file_type or not path:
                raise ValueError(
                    f"Index file row {i} has an empty name, type, or path."
                )

            rows.append(
                {
                    "name": name,
                    "type": file_type,
                    "path": path,
                    "row_number": str(i),
                }
            )

    return rows


def resolve_path(raw_path: str) -> Path:
    """
    Resolve file paths relative to the index file location.
    """
    path = Path(raw_path)
    if path.is_absolute():
        return path
    return top_dir / path


def read_fasta_alignment(path: Path) -> Dict[str, str]:
    """
    Read an aligned FASTA file and return {taxon: sequence}.
    """
    if not path.exists():
        raise FileNotFoundError(f"Alignment file not found: {path}")

    alignment = AlignIO.read(str(path), "fasta")

    if len(alignment) == 0:
        raise ValueError(f"No sequences found in alignment: {path}")

    aln_length = alignment.get_alignment_length()
    if aln_length == 0:
        raise ValueError(f"Alignment has zero length: {path}")

    sequences = {}

    for record in alignment:
        taxon = str(record.id).strip()
        seq = str(record.seq).upper()

        if not taxon:
            raise ValueError(f"Found empty taxon name in: {path}")

        if taxon in sequences:
            raise ValueError(f"Duplicate taxon '{taxon}' in alignment: {path}")

        if len(seq) != aln_length:
            raise ValueError(
                f"Sequence '{taxon}' in {path} has length {len(seq)}, "
                f"expected {aln_length}."
            )

        sequences[taxon] = seq

    return sequences


def read_newick_tree(path: Path) -> str:
    """
    Read a Newick tree from a file.
    """
    if not path.exists():
        raise FileNotFoundError(f"Tree file not found: {path}")

    text = path.read_text().strip()

    if not text:
        raise ValueError(f"Tree file is empty: {path}")

    # Remove comments or blank lines only if the entire file has simple line breaks.
    # This preserves normal Newick content.
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    text = " ".join(lines)

    if not text.endswith(";"):
        text += ";"

    return text


def build_concatenated_matrix(
    alignment_entries: List[Tuple[str, Path]]
) -> Tuple[List[str], Dict[str, str], List[Tuple[str, int, int]]]:
    """
    Concatenate alignments and return:
        taxa
        concatenated matrix {taxon: sequence}
        charsets [(name, start, end)]
    """
    all_taxa = []
    seen_taxa = set()

    parsed_alignments = []

    for name, path in alignment_entries:
        seqs = read_fasta_alignment(path)
        aln_len = len(next(iter(seqs.values())))

        for taxon in seqs:
            if taxon not in seen_taxa:
                seen_taxa.add(taxon)
                all_taxa.append(taxon)

        parsed_alignments.append((name, seqs, aln_len))

    if not parsed_alignments:
        raise ValueError("No alignment files were provided.")

    matrix = {taxon: "" for taxon in all_taxa}
    charsets = []

    start = 1

    for name, seqs, aln_len in parsed_alignments:
        end = start + aln_len - 1
        charsets.append((clean_partition_name(name), start, end))

        for taxon in all_taxa:
            if taxon in seqs:
                matrix[taxon] += seqs[taxon]
            else:
                matrix[taxon] += "?" * aln_len

        start = end + 1

    return all_taxa, matrix, charsets


def write_nexus(
    output_file: Path,
    taxa: List[str],
    matrix: Dict[str, str],
    charsets: List[Tuple[str, int, int]],
    tree_entries: List[Tuple[str, Path]],
    datatype: str = "DNA",
    interleave: bool = False,
    block_size: int = 80,
) -> None:
    """
    Write the final NEXUS file.
    """
    nchar = len(next(iter(matrix.values())))

    with output_file.open("w") as out:
        out.write("#NEXUS\n\n")

        out.write("BEGIN TAXA;\n")
        out.write(f"    DIMENSIONS NTAX={len(taxa)};\n")
        out.write("    TAXLABELS\n")
        for taxon in taxa:
            out.write(f"        {quote_nexus_name(taxon)}\n")
        out.write("    ;\n")
        out.write("END;\n\n")

        out.write("BEGIN CHARACTERS;\n")
        out.write(f"    DIMENSIONS NCHAR={nchar};\n")
        out.write(f"    FORMAT DATATYPE={datatype} MISSING=? GAP=-;\n")
        out.write("    MATRIX\n")

        if interleave:
            for start in range(0, nchar, block_size):
                end = min(start + block_size, nchar)
                for taxon in taxa:
                    out.write(
                        f"        {quote_nexus_name(taxon)}  "
                        f"{matrix[taxon][start:end]}\n"
                    )
                out.write("\n")
        else:
            for taxon in taxa:
                out.write(f"        {quote_nexus_name(taxon)}  {matrix[taxon]}\n")

        out.write("    ;\n")
        out.write("END;\n\n")

        out.write("BEGIN SETS;\n")
        for name, start, end in charsets:
            out.write(f"    CHARSET {name} = {start}-{end};\n")
        out.write("END;\n\n")

        if tree_entries:
            out.write("BEGIN TREES;\n")
            for name, path in tree_entries:
                tree_name = clean_partition_name(name)
                tree = read_newick_tree(path)
                out.write(f"    TREE {tree_name} = {tree}\n")
            out.write("END;\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Combine multiple FASTA alignments and Newick trees into one NEXUS file."
    )

    parser.add_argument(
        "-i",
        "--index",
        required=True,
        help="CSV or TSV file with columns: name,type,path",
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output NEXUS file.",
    )

    parser.add_argument(
        "--datatype",
        default="DNA",
        choices=["DNA", "RNA", "PROTEIN", "STANDARD"],
        help="NEXUS datatype. Default: DNA",
    )

    parser.add_argument(
        "--interleave",
        action="store_true",
        help="Write the alignment matrix in interleaved blocks.",
    )

    parser.add_argument(
        "--block-size",
        type=int,
        default=80,
        help="Block size for interleaved output. Default: 80",
    )

    args = parser.parse_args()

    index_file = Path(args.index).resolve()
    output_file = Path(args.output).resolve()

    try:
        rows = read_index(index_file)

        alignment_entries = []
        tree_entries = []

        for row in rows:
            name = row["name"]
            file_type = row["type"]
            path = resolve_path(row["path"])

            if file_type in VALID_ALIGNMENT_TYPES:
                alignment_entries.append((name, path))
            elif file_type in VALID_TREE_TYPES:
                tree_entries.append((name, path))
            else:
                raise ValueError(
                    f"Unknown type '{file_type}' on index row {row['row_number']}. "
                    f"Use 'alignment' or 'tree'."
                )

        taxa, matrix, charsets = build_concatenated_matrix(alignment_entries)

        write_nexus(
            output_file=output_file,
            taxa=taxa,
            matrix=matrix,
            charsets=charsets,
            tree_entries=tree_entries,
            datatype=args.datatype,
            interleave=args.interleave,
            block_size=args.block_size,
        )

        print(f"Wrote NEXUS file: {output_file}")
        print(f"Number of taxa: {len(taxa)}")
        print(f"Total alignment length: {len(next(iter(matrix.values())))}")
        print(f"Number of partitions: {len(charsets)}")
        print(f"Number of trees: {len(tree_entries)}")

    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
