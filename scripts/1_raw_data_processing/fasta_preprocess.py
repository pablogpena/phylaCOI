# -*- coding: utf-8 -*-
# +
# #!/usr/bin/env python3
# -

"""
Preprocess a FASTA file:
1. Extracts sequence headers into a text file.
2. Cleans FASTA identifiers by keeping only the sequence name.

Example:
    python fasta_preprocess.py -i input.fasta -t headers.txt -o clean.fasta
"""

import argparse
from pathlib import Path
from Bio import SeqIO

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Preprocess a FASTA file: extract headers and clean identifiers."
    )
    parser.add_argument(
        "-i", "--input", required=True, type=Path,
        help="Input FASTA file"
    )
    parser.add_argument(
        "-t", "--txt", required=True, type=Path,
        help="Output TXT file with sequence headers"
    )
    parser.add_argument(
        "-o", "--output", required=True, type=Path,
        help="Output FASTA file with cleaned identifiers"
    )

    args = parser.parse_args()

    if not args.input.exists():
        print(f"[ERROR] Input file not found: {args.input}")
        exit(1)

    # 1. Extract headers
    with args.input.open("r") as fasta_file, args.txt.open("w") as txt_file:
        for line in fasta_file:
            if line.startswith(">"):
                txt_file.write(line.strip()[1:] + "\n")
    print(f"[OK] Headers extracted to: {args.txt}")

    # 2. Clean FASTA identifiers
    records = []
    for record in SeqIO.parse(str(args.input), "fasta"):
        record.id = record.id.split(" ")[0]   
        record.description = ""              
        records.append(record)

    SeqIO.write(records, str(args.output), "fasta")
    print(f"[OK] Clean FASTA created at: {args.output}")
