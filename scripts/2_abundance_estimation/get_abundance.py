# -*- coding: utf-8 -*-
# +
# -*- coding: utf-8 -*-
# -

"""
get_abundance.py

Iterates through subfolders (one per phylum) and executes the abundance and alignment
pipeline defined in `fasta_processing.py`. For each phylum, the script locates the
FASTA, mapping, and metadata files, generates a complete abundance table with
abundance values and coordinates, produces a unique FASTA file, and optionally
runs MAFFT for sequence alignment.

Usage
-----
python scripts/get_abundance/get_abundance.py \
    -f /path/to/root_folder \
    -n nombres_secuencias.txt \
    -m metadata.csv \
    [--skip-mafft]
"""

import os
import sys
from pathlib import Path
import argparse

sys.path.append('/workspace/PhylaCOI')
from src.fasta_processing.fasta_processing import run_pipeline

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run abundance and alignment pipeline for all phylum folders."
    )
    parser.add_argument(
        "-f", "--folder", required=True, type=Path,
        help="Root folder containing subfolders (one per phylum)."
    )
    parser.add_argument(
        "-n", "--names", required=True, type=Path,
        help="Full path to the sequence mapping file (e.g. '/data/raw/seq_headers.txt')."
    )
    parser.add_argument(
        "-m", "--metadata", required=True, type=Path,
        help="Full path to the metadata file (e.g. '/data/raw/KOI_metadata.csv')."
    )
    parser.add_argument(
        "-a", "--mafft", required=True, type=int, choices=[0, 1],
        help="Run MAFFT alignment: 1 = run MAFFT, 0 = skip."
    )

    args = parser.parse_args()
    base_dir = args.folder
    names_path = args.names
    meta_path = args.metadata
    run_mafft_flag = bool(args.mafft)

    if not base_dir.exists():
        print(f"[ERROR] Input folder not found: {base_dir}")
        sys.exit(1)
    if not names_path.exists():
        print(f"[ERROR] Sequence mapping file not found: {names_path}")
        sys.exit(1)
    if not meta_path.exists():
        print(f"[ERROR] Metadata file not found: {meta_path}")
        sys.exit(1)

    print(f"Starting pipeline in root folder: {base_dir}")

    total = 0
    processed = 0
    failed = 0

    for phylum_folder in os.listdir(base_dir):
        phylum_path = base_dir / phylum_folder

        if not phylum_path.is_dir():
            continue

        total += 1
        print(f"\nProcessing phylum: {phylum_folder}")

        fasta_files = [f for f in phylum_path.glob("*.fasta")]
        if not fasta_files:
            print(f"  No FASTA file found in {phylum_folder}")
            failed += 1
            continue

        fasta_path = fasta_files[0]
        output_dir = phylum_path / "output"
        output_dir.mkdir(exist_ok=True)

        print(f"  Running pipeline for {phylum_folder}")
        try:
            run_pipeline(
                fasta_path=fasta_path,
                names_path=names_path,
                meta_path=meta_path,
                out_dir=output_dir,
                run_mafft_flag=run_mafft_flag
            )
            processed += 1
        except Exception as e:
            print(f"  Error while processing {phylum_folder}: {e}")
            failed += 1

    print(f"\nCompleted processing.")
    print(f"Total folders: {total}")
    print(f"Successfully processed: {processed}")
    print(f"Failed: {failed}")
