# -*- coding: utf-8 -*-
from __future__ import annotations

import re
import subprocess
from pathlib import Path
from typing import List, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def extract_id_and_localities(line: str) -> Tuple[str | None, str | None]:
    id_match  = re.match(r"(\S+)", line)
    loc_match = re.search(r"merged_sample=\{(.*)\}", line)
    seq_id = id_match.group(1) if id_match else None
    locs   = loc_match.group(1) if loc_match else None
    return seq_id, locs

def run_mafft(in_fasta: Path, out_fasta: Path) -> None:
    cmd: List[str] = [
        "mafft", "--auto", "--ep", "1.5", "--op", "3.0",
        "--maxiterate", "0", "--large", str(in_fasta)
    ]
    with out_fasta.open("w") as fout:
        proc = subprocess.run(cmd, stdout=fout, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"MAFFT exited with return code {proc.returncode}\n{proc.stderr}"
        )

def run_pipeline(fasta_path: Path, names_path: Path, meta_path: Path, out_dir: Path, run_mafft_flag: bool = True) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    fasta_sequences = SeqIO.to_dict(SeqIO.parse(str(fasta_path), "fasta"))
    fasta_ids = set(fasta_sequences.keys())
    print(f"Read {len(fasta_ids):,} sequences from {fasta_path.name}")

    with names_path.open() as fh:
        name_lines = [ln.strip() for ln in fh if ln.strip()]

    filtered_lines = [ln for ln in name_lines if any(fid in ln for fid in fasta_ids)]
    print(f"Lines remaining after FASTA‑ID filter: {len(filtered_lines):,}")

    abundance_rows: List[dict] = []
    for line in filtered_lines:
        seq_id, locs = extract_id_and_localities(line)
        if not seq_id or not locs:
            continue
        for loc in locs.split(","):
            loc_name = re.sub(r"'(.*?):.*", r"\1", loc).strip().upper()
            abundance_rows.append({"ID": seq_id, "Localities": loc_name})

    abund_df = pd.DataFrame(abundance_rows)
    abund_df["Localities"] = abund_df["Localities"].str.strip().str.replace("'", "")

    meta_df = pd.read_csv(meta_path, sep=";", encoding="latin1")
    meta_df["id_sample"] = meta_df["id_sample"].str.upper()

    coords = (
        meta_df.set_index("id_sample")["coordenates2"]
        .str.replace(" ", "")
        .str.split(",", expand=True)
        .rename(columns={0: "lat", 1: "lon"})
        .apply(pd.to_numeric, errors="coerce")
    )

    abund_df = abund_df.join(coords, on="Localities")

    abund_csv = out_dir / "abundances.csv"
    abund_df.to_csv(abund_csv, index=False)
    print(f"Wrote {abund_csv.name}")

    abund_unique = abund_df.drop_duplicates(subset=["ID", "lat", "lon"]).copy()

    dup_mask = abund_unique["ID"].duplicated(keep=False)
    abund_unique["UniqueID"] = abund_unique["ID"]
    abund_unique.loc[dup_mask, "UniqueID"] = (
        abund_unique.loc[dup_mask, "ID"]
        + "_dup"
        + (abund_unique.groupby("ID").cumcount() + 1).astype(str)
    )

    abund_unique_csv = out_dir / "abundances_unique.csv"
    abund_unique.to_csv(abund_unique_csv, index=False)
    print(f"Wrote {abund_unique_csv.name}")

    unique_records: List[SeqRecord] = []
    for _, row in abund_unique.iterrows():
        orig_id = re.sub(r"_dup\d+$", "", row["UniqueID"])
        if orig_id in fasta_sequences:
            rec = SeqRecord(
                fasta_sequences[orig_id].seq,
                id=row["UniqueID"],
                description=""
            )
            unique_records.append(rec)

    unique_fasta = out_dir / "unique_sequences.fasta"
    SeqIO.write(unique_records, unique_fasta, "fasta")
    print(f"Wrote FASTA with unique IDs → {unique_fasta.name}")

    if run_mafft_flag:
        aligned_fasta = out_dir / "aligned_sequences_mafft.fasta"
        try:
            run_mafft(unique_fasta, aligned_fasta)
            print(f"Alignment OK → {aligned_fasta.name}")
        except RuntimeError as err:
            print(f"[ERROR] MAFFT could not complete:\n{err}")
    else:
        print("Alignment step skipped (run_mafft_flag = False).")
