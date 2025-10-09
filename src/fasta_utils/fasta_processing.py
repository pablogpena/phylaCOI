# -*- coding: utf-8 -*-
# +
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import re
import subprocess
from pathlib import Path
from typing import List, Tuple


# -

def _extract_id_and_localities(line):
    """
    Extract sequence ID and localities from a text line.

    Returns
    -------
    tuple
        (seq_id, locs) where each element can be a string or None.
    """
    id_match = re.match(r"(\S+)", line)                     # First token = sequence ID
    loc_match = re.search(r"merged_sample=\{(.*)\}", line)  # Content inside merged_sample={}
    seq_id = id_match.group(1) if id_match else None
    locs = loc_match.group(1) if loc_match else None
    return seq_id, locs


def _run_mafft(in_fasta, out_fasta):
    """Run MAFFT on a FASTA file and save the aligned output. Raises RuntimeError if MAFFT fails."""
    cmd = [
        "mafft", "--auto", "--ep", "1.5", "--op", "3.0",
        "--maxiterate", "0", "--large", str(in_fasta)
    ]
    with out_fasta.open("w") as fout:
        proc = subprocess.run(cmd, stdout=fout, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"MAFFT exited with return code {proc.returncode}\n{proc.stderr}"
        )


def run_pipeline(fasta_path, names_path, meta_path, out_dir, run_mafft_flag=True):
    """Process FASTA, mapping, and metadata files to generate abundance tables (with abundances), unique sequences, and optional MAFFT alignment."""
    out_dir.mkdir(parents=True, exist_ok=True)

    # Read FASTA and collect IDs
    fasta_sequences = SeqIO.to_dict(SeqIO.parse(str(fasta_path), "fasta"))
    fasta_ids = set(fasta_sequences)
    print(f"Read {len(fasta_ids):,} sequences from {fasta_path.name}")

    # Read mapping names and filter by FASTA IDs
    with names_path.open() as fh:
        filtered_lines = [ln.strip() for ln in fh if ln.strip() and any(fid in ln for fid in fasta_ids)]
    print(f"Lines remaining after FASTA-ID filter: {len(filtered_lines):,}")

    # Build abundance table 
    abundance_rows = []
    for line in filtered_lines:
        seq_id, locs = _extract_id_and_localities(line)
        if seq_id and locs:
            # Extract pairs: 'LOC1':15, 'LOC2':8, ...
            for loc, reads in re.findall(r"'([^']+)'\s*:\s*(\d+)", locs):
                loc_name = loc.upper().strip()
                abundance_rows.append({
                    "ID": seq_id,
                    "Localities": loc_name,
                    "Abundance": int(reads)
                })

    abund_df = pd.DataFrame(abundance_rows)

    # Read metadata and extract coordinates
    meta_df = pd.read_csv(meta_path, sep=";", encoding="latin1")
    meta_df["id_sample"] = meta_df["id_sample"].str.upper()
    coords = (
        meta_df.set_index("id_sample")["coordenates2"]
        .str.replace(" ", "")
        .str.split(",", expand=True)
        .rename(columns={0: "lat", 1: "lon"})
        .apply(pd.to_numeric, errors="coerce")
    )

    # Merge abundance table with coordinates
    abund_df = abund_df.join(coords, on="Localities")

    # Drop exact duplicates (same ID + coordinates)
    abund_unique = abund_df.drop_duplicates(subset=["ID", "lat", "lon"]).copy()

    # Create UniqueID column
    abund_unique["UniqueID"] = abund_unique["ID"]
    dup_counts = abund_unique.groupby("ID").cumcount()
    abund_unique.loc[abund_unique["ID"].duplicated(keep=False), "UniqueID"] += "_dup" + (dup_counts + 1).astype(str)

    # Save a single, well-formed abundance table
    abund_csv = out_dir / "abundances.csv"
    abund_unique.to_csv(abund_csv, index=False)
    print(f"Wrote {abund_csv.name} (with Abundance and UniqueID columns)")

    # Write unique FASTA file
    unique_records = [
        SeqRecord(fasta_sequences[re.sub(r"_dup\d+$", "", uid)].seq, id=uid, description="")
        for uid in abund_unique["UniqueID"]
        if re.sub(r"_dup\d+$", "", uid) in fasta_sequences
    ]

    unique_fasta = out_dir / "unique_sequences.fasta"
    SeqIO.write(unique_records, unique_fasta, "fasta")
    print(f"Wrote FASTA with unique IDs → {unique_fasta.name}")

    # Run MAFFT alignment
    if run_mafft_flag:
        aligned_fasta = out_dir / "aligned_sequences_mafft.fasta"
        try:
            _run_mafft(unique_fasta, aligned_fasta)
            print(f"Alignment OK → {aligned_fasta.name}")
        except RuntimeError as err:
            print(f"[ERROR] MAFFT could not complete:\n{err}")
    else:
        print("Alignment step skipped (run_mafft_flag = False).")

