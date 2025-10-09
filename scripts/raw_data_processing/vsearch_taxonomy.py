# -*- coding: utf-8 -*-
# +
# #!/usr/bin/env python3
# -

"""
Pipeline to classify FASTA sequences against a reference database using vsearch.

Steps:
1. Run vsearch (--usearch_global) for each FASTA file in the input folder (except the reference DB).
2. Save raw results in TXT and Excel formats.
3. Filter hits by identity threshold.
4. Group sequences by the 5th taxonomic rank and write group-specific FASTA files.
5. Save filtered results in Excel.

Usage:
    python vsearch_taxonomy.py -d reference.fasta -i 0.84 -f /path/to/fasta_folder -o /path/to/output_folder
"""

import subprocess
import argparse
import pandas as pd
from pathlib import Path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Classify FASTA sequences against a reference database using vsearch."
    )
    parser.add_argument(
        "-d", "--db", required=True, type=Path,
        help="Reference database FASTA file"
    )
    parser.add_argument(
        "-i", "--identity", required=True, type=float,
        help="Minimum identity threshold (e.g. 0.84)"
    )
    parser.add_argument(
        "-f", "--folder", required=True, type=Path,
        help="Folder containing FASTA files to process"
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=None,
        help="Optional output base folder for results (default: alongside each FASTA)"
    )
    args = parser.parse_args()

    if not args.db.exists():
        print(f"[ERROR] Reference database not found: {args.db}")
        exit(1)
    if not args.folder.exists():
        print(f"[ERROR] Input folder not found: {args.folder}")
        exit(1)

    # Get all FASTA files in the given folder except the reference DB
    fasta_files = [f for f in args.folder.glob("*.fasta") if f.name != args.db.name]

    if not fasta_files:
        print(f"[WARNING] No FASTA files found in {args.folder}")
        exit(0)

    for fasta_file in fasta_files:
        # Decide output folder: central (-o) or next to FASTA
        if args.output:
            folder = args.output / fasta_file.stem
        else:
            folder = fasta_file.parent / fasta_file.stem

        folder.mkdir(parents=True, exist_ok=True)

        # Output paths
        txt_out = folder / f"{fasta_file.stem}_results.txt"
        xlsx_out = folder / f"{fasta_file.stem}_results.xlsx"
        xlsx_filtered = folder / f"{fasta_file.stem}_results_filtered.xlsx"

        # Run vsearch
        cmd = [
            "vsearch",
            "--usearch_global", str(fasta_file),
            "--db", str(args.db),
            "--id", str(args.identity),
            "--blast6out", str(txt_out),
            "--quiet"
        ]
        subprocess.run(cmd, check=True)

        # Parse vsearch results
        results = []
        with open(txt_out, "r") as f:
            for line in f:
                fields = line.strip().split("\t")
                query, identity, subject = fields[0], float(fields[2]), fields[1]
                results.append((query, identity, subject))
        
        df = pd.DataFrame(results, columns=["Query", "Identity", "Subject"])
        df.to_excel(xlsx_out, index=False)

        # Filter results
        df_filtered = df[df["Identity"] >= args.identity].copy()
        df_filtered["Subject"] = df_filtered["Subject"].str.split(";")

        # Group by taxonomic rank 5 (index 4)
        for group_name, group_df in df_filtered.groupby(df_filtered["Subject"].str[4]):
            phylum_folder = folder / str(group_name)
            phylum_folder.mkdir(parents=True, exist_ok=True)
            
            # Write FASTA file inside the phylum folder
            fasta_out = phylum_folder / f"{group_name}_{fasta_file.stem}.fasta"
            seq_names = set(group_df["Query"].tolist())

            with open(fasta_file, "r") as original, open(fasta_out, "w") as out:
                for line in original:
                    if line.startswith(">"):
                        name = line.strip()[1:]
                        if name in seq_names:
                            out.write(line)
                            out.write(next(original))

        df_filtered.to_excel(xlsx_filtered, index=False)

        print(f"[OK] Processed {fasta_file.name}. Results in {folder}/")
