#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Run OTU generation for one phylum folder."""

import argparse
import sys
from pathlib import Path


def default_repo_dir() -> Path:
    """Return the repository root from this workflow helper location."""
    return Path(__file__).resolve().parents[3]


def main():
    parser = argparse.ArgumentParser(
        description="Generate OTU outputs for a single phylum."
    )
    parser.add_argument("--abundance-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--identity", type=float, default=0.97)
    parser.add_argument("--no-vsearch", action="store_true")
    parser.add_argument("--repo-dir", type=Path, default=default_repo_dir())
    args = parser.parse_args()

    repo_dir = args.repo_dir.resolve()
    sys.path.insert(0, str(repo_dir))

    from src.otu_utils.otu_processing import process_phylum

    if not args.abundance_dir.exists():
        raise SystemExit(f"Abundance folder not found: {args.abundance_dir}")

    process_phylum(
        input_dir=args.abundance_dir,
        output_dir=args.output_dir,
        identity=args.identity,
        run_vsearch_flag=not args.no_vsearch,
    )


if __name__ == "__main__":
    main()
