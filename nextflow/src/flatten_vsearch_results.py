#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Flatten VSEARCH taxonomy outputs into one directory per phylum."""

import argparse
import shutil
from pathlib import Path


def find_phylum_dirs(input_root: Path):
    """Return directories that contain phylum FASTA files."""
    candidates = []
    for folder in input_root.rglob("*"):
        if not folder.is_dir():
            continue
        fasta_files = list(folder.glob("*.fasta"))
        if fasta_files:
            candidates.append(folder)
    return sorted(candidates, key=lambda path: path.name)


def main():
    parser = argparse.ArgumentParser(
        description="Flatten nested VSEARCH results into one folder per phylum."
    )
    parser.add_argument("-i", "--input", required=True, type=Path)
    parser.add_argument("-o", "--output", required=True, type=Path)
    args = parser.parse_args()

    if not args.input.exists():
        raise SystemExit(f"Input folder not found: {args.input}")

    args.output.mkdir(parents=True, exist_ok=True)
    phylum_dirs = find_phylum_dirs(args.input)

    if not phylum_dirs:
        raise SystemExit(f"No phylum FASTA folders found under: {args.input}")

    for phylum_dir in phylum_dirs:
        destination = args.output / phylum_dir.name
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(phylum_dir, destination, symlinks=False)
        print(f"Copied phylum folder: {phylum_dir} -> {destination}")

    print(f"Flattened {len(phylum_dirs)} phylum folders into: {args.output}")


if __name__ == "__main__":
    main()
