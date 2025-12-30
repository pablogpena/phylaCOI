# -*- coding: utf-8 -*-
# +
# -*- coding: utf-8 -*-
# -

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[2]
sys.path.append(str(PROJECT_ROOT))

from src.otu_utils.otu_processing import process_phylum

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate OTUs per folder by cleaning FASTA files and performing clustering with VSEARCH."
    )
    parser.add_argument(
        "-a", "--abundance",
        type=Path,
        required=True,
        help="Root directory containing per-phylum abundance outputs."
    )
    parser.add_argument(
        "-o", "--otus",
        type=Path,
        required=True,
        help="Root directory where per-phylum OTU outputs will be written."
    )
    parser.add_argument(
        "-i", "--identity",
        type=float,
        default=0.97,
        help="Sequence identity threshold for VSEARCH clustering (default: 0.97)."
    )
    parser.add_argument(
        "--no-vsearch",
        action="store_true",
        help="If set, skips the VSEARCH execution (only performs FASTA cleaning)."
    )

    args = parser.parse_args()
    abundance_root = args.abundance
    otus_root = args.otus
    run_vsearch_flag = not args.no_vsearch

    if not abundance_root.exists():
        raise SystemExit(f"[ERROR] Input folder not found: {abundance_root}")

    otus_root.mkdir(parents=True, exist_ok=True)

    for phylum_dir in abundance_root.iterdir():
        if not phylum_dir.is_dir():
            continue
        print(f"Processing: {phylum_dir.name}")
        output_dir = otus_root / phylum_dir.name
        process_phylum(phylum_dir, output_dir, args.identity, run_vsearch_flag)
