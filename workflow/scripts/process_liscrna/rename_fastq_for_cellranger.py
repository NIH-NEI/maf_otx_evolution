#!/usr/bin/env python3
from __future__ import annotations
import argparse
import collections
import csv
import gzip
import os
import re
import sys
from pathlib import Path
from typing import Dict, Tuple, List, Optional

GSM_DIR_RE = re.compile(r"^(GSM\d+)$")
SRR_FASTQ_RE = re.compile(r"^(SRR\d+)_(1|2)\.fastq\.gz$")

# Illumina core: instrument:run:flowcell:lane:tile:x:y  (lane in group 2)
ILLUMINA_LANE_RE = re.compile(r"\b[A-Za-z0-9]+:\d+:[A-Za-z0-9]+:(\d+):\d+:\d+:\d+\b")

# BGI/DNBSEQ-ish token with lane like "...L1C001R001..."
BGI_LANE_RE = re.compile(r"L(\d+)C")


def sanitize_sample_name(name: str) -> str:
    # Cell Ranger is happiest with A-Z a-z 0-9 _ -
    name = name.strip()
    name = re.sub(r"\s+", "_", name)
    name = re.sub(r"[^A-Za-z0-9_\-]", "_", name)
    name = re.sub(r"_+", "_", name).strip("_")
    if not name:
        raise ValueError("Sample name became empty after sanitization.")
    return name

def read_gsm_map(tsv_path: Path) -> Dict[str, str]:
    import pandas as pd
    gsm_to_sample = (
        pd.read_csv(tsv_path, sep = "\t")
        .set_index("GSM")["SampleID"].to_dict()
    )
    if not gsm_to_sample:
        raise ValueError(f"No mappings loaded from {tsv_path}")
    return gsm_to_sample

def extract_lane_from_header(header_line: str) -> Optional[int]:
    """
    header examples:
      @SRR... A00265:1137:HGKLCDSX5:4:1101:1199:1000 length=150
      @SRR... V350140183L1C001R0010000011 length=20
    """
    m = ILLUMINA_LANE_RE.search(header_line)
    if m:
        return int(m.group(1))

    # Try BGI: lane appears as L1, L2, etc in the second token
    # Use a broad search, but prefer matches in the token after SRR.
    parts = header_line.strip().split()
    if len(parts) >= 2:
        #print("Trying BGI on part2 of header", parts[1])
        m2 = BGI_LANE_RE.search(parts[1])
        if m2:
            #print("Found match")
            return int(m2.group(1))
        else:
            print("No match")

    # Fallback: any L\d+ anywhere
    m3 = BGI_LANE_RE.search(header_line)
    if m3:
        return int(m3.group(1))

    return None

def infer_lane_from_r1(r1_path: Path, n_reads: int = 100) -> Optional[int]:
    lane_counts = collections.Counter()
    try:
        with gzip.open(r1_path, "rt", encoding="utf-8", errors="replace") as fh:
            for _ in range(n_reads):
                header = fh.readline()
                if not header:
                    break
                # skip remaining 3 lines of FASTQ record
                fh.readline(); fh.readline(); fh.readline()
                lane = extract_lane_from_header(header)
                if lane is not None:
                    lane_counts[lane] += 1
    except Exception as e:
        print(f"WARNING: could not infer lane from {r1_path}: {e}", file=sys.stderr)
        return None

    if not lane_counts:
        return None
    return lane_counts.most_common(1)[0][0]

def safe_symlink(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        # If already correct link, leave it
        try:
            if dest.is_symlink() and dest.resolve() == src.resolve():
                return
        except Exception:
            pass
        raise FileExistsError(f"Destination exists: {dest}")
    os.symlink(src, dest)

def safe_rename(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        raise FileExistsError(f"Destination exists: {dest}")
    src.rename(dest)

def main():
    ap = argparse.ArgumentParser(
        description="Create Cell Ranger-compatible FASTQ symlinks/renames from SRA FASTQs grouped by GSM."
    )
    ap.add_argument("--root", required=True, type=Path,
                    help="Root directory containing GSM*_SRR* subdirectories.")
    ap.add_argument("--map_tsv", required=True, type=Path,
                    help="TSV mapping GSM to sample name (2 columns or header).")
    ap.add_argument("--out", required=True, type=Path,
                    help="Output directory for Cell Ranger-compatible FASTQs.")
    ap.add_argument("--mode", choices=["symlink", "rename"], default="symlink",
                    help="Create symlinks (recommended) or rename/move files.")
    ap.add_argument("--reads_for_lane", type=int, default=100,
                    help="Number of reads (records) to sample from R1 for lane inference.")
    ap.add_argument("--default_lane", type=int, default=None,
                    help="If lane cannot be inferred, use this lane number (e.g. 1).")
    ap.add_argument("--dry_run", action="store_true",
                    help="Print actions without making changes.")
    args = ap.parse_args()

    gsm_to_sample = read_gsm_map(args.map_tsv)

    jobs: List[Tuple[str, str, Path]] = []

    gsm_dir = args.root

    m = GSM_DIR_RE.match(gsm_dir.name)
    if not m:
        print(f"ERROR: No GSM* directories found under {args.root}", file=sys.stderr)
        sys.exit(2)

    gsm = m.group(1)

    # find SRRs inside GSM directory
    srrs = set()
    for f in args.root.iterdir():
        if not f.is_file():
            continue
        m2 = SRR_FASTQ_RE.match(f.name)
        if m2:
            srrs.add(m2.group(1))

    for srr in sorted(srrs):
        jobs.append((gsm, srr, gsm_dir))

    if not jobs:
        print(f"ERROR: No GSM*_SRR* directories found under {args.root}", file=sys.stderr)
        sys.exit(2)

    # For Cell Ranger naming uniqueness: increment per (sample,lane)
    counter_per_sample_lane: Dict[Tuple[str, int], int] = collections.defaultdict(int)
    summary = []
    errors = 0

    for gsm, srr, d in jobs:
        sample = gsm_to_sample.get(gsm)
        if not sample:
            print(f"WARNING: GSM {gsm} not found in mapping; skipping {d.name}", file=sys.stderr)
            continue

        r1 = d / f"{srr}_1.fastq.gz"
        r2 = d / f"{srr}_2.fastq.gz"
        if not r1.exists() or not r2.exists():
            print(f"WARNING: Missing FASTQs for {d.name}; skipping", file=sys.stderr)
            continue

        lane = infer_lane_from_r1(r1, n_reads=args.reads_for_lane)
        if lane is None:
            if args.default_lane is None:
                print(f"WARNING: Could not infer lane for {r1}; no --default_lane. Skipping.", file=sys.stderr)
                continue
            lane = args.default_lane

        key = (sample, lane)
        counter_per_sample_lane[key] += 1
        idx = counter_per_sample_lane[key]  # 1-based
        idx_str = f"{idx:03d}"
        lane_str = f"{lane:03d}"

        out_r1 = args.out / f"{sample}_S1_L{lane_str}_R1_{idx_str}.fastq.gz"
        out_r2 = args.out / f"{sample}_S1_L{lane_str}_R2_{idx_str}.fastq.gz"

        action = safe_symlink if args.mode == "symlink" else safe_rename

        try:
            if args.dry_run:
                print(f"{args.mode.upper()}: {r1} -> {out_r1} (lane={lane})")
                print(f"{args.mode.upper()}: {r2} -> {out_r2} (lane={lane})")
            else:
                action(r1, out_r1)
                action(r2, out_r2)

            summary.append({
                "GSM": gsm,
                "SRR": srr,
                "sample_name": sample,
                "lane": lane,
                "src_r1": str(r1),
                "src_r2": str(r2),
                "out_r1": str(out_r1),
                "out_r2": str(out_r2),
                "mode": args.mode
            })
        except Exception as e:
            errors += 1
            print(f"ERROR: {d.name}: {e}", file=sys.stderr)

    summary_path = args.out / "cellranger_fastq_links.summary.tsv"
    if args.dry_run:
        print(f"\n(DRY RUN) Would write summary to: {summary_path}")
    else:
        if summary:
            args.out.mkdir(parents=True, exist_ok=True)
            with summary_path.open("w", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=list(summary[0].keys()), delimiter="\t")
                writer.writeheader()
                for row in summary:
                    writer.writerow(row)

    if errors:
        print(f"\nCompleted with {errors} errors. Summary: {summary_path}", file=sys.stderr)
        sys.exit(1)

    print(f"\nDone. Wrote {len(summary)} entries. Summary: {summary_path}")

if __name__ == "__main__":
    main()

