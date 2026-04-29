#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Run abundance estimation for one phylum folder."""

import argparse
import sys
from pathlib import Path


def default_repo_dir() -> Path:
    """Return the repository root from this workflow helper location."""
    return Path(__file__).resolve().parents[3]


def main():
    parser = argparse.ArgumentParser(
        description="Generate abundance outputs for a single phylum."
    )
    parser.add_argument("--phylum-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--names", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--mafft", required=True, type=int, choices=[0, 1])
    parser.add_argument("--repo-dir", type=Path, default=default_repo_dir())
    args = parser.parse_args()

    repo_dir = args.repo_dir.resolve()
    sys.path.insert(0, str(repo_dir))

    from src.fasta_utils.fasta_processing import run_pipeline

    if not args.phylum_dir.exists():
        raise SystemExit(f"Phylum folder not found: {args.phylum_dir}")
    if not args.names.exists():
        raise SystemExit(f"Sequence names file not found: {args.names}")
    if not args.metadata.exists():
        raise SystemExit(f"Sample metadata file not found: {args.metadata}")

    fasta_files = sorted(args.phylum_dir.glob("*.fasta"))
    if not fasta_files:
        raise SystemExit(f"No FASTA file found in: {args.phylum_dir}")
    if len(fasta_files) > 1:
        print(f"Multiple FASTA files found in {args.phylum_dir}; using {fasta_files[0].name}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    run_pipeline(
        fasta_path=fasta_files[0],
        names_path=args.names,
        meta_path=args.metadata,
        out_dir=args.output_dir,
        run_mafft_flag=bool(args.mafft),
    )


if __name__ == "__main__":
    main()
