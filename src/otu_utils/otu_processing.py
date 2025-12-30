# -*- coding: utf-8 -*-
from pathlib import Path
import subprocess


def _clean_fasta(src: Path, dest: Path):
    """
    Remove gap characters from an aligned FASTA file.
    Writes the cleaned output to the destination path.
    """
    with src.open() as fin, dest.open("w") as fout:
        for line in fin:
            fout.write(line if line.startswith(">") else line.replace("-", ""))
    print(f"Cleaned FASTA written to: {dest}")


def _run_vsearch(fasta: Path, centroids: Path, uc_file: Path, identity: float):
    """
    Run VSEARCH to cluster sequences into OTUs.
    Writes centroid and UC outputs to the provided paths.
    """
    cmd = [
        "vsearch",
        "--cluster_fast", str(fasta),
        "--id", f"{identity:.2f}",
        "--centroids", str(centroids),
        "--uc", str(uc_file),
    ]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"VSEARCH error ({proc.returncode}):\n{proc.stderr}")
    print(f"VSEARCH completed successfully: {centroids}, {uc_file}")


def _build_mapping(uc_file: Path, dest: Path):
    """
    Build a two-column OTU mapping from a VSEARCH .uc file.
    Writes the mapping to the destination path.
    """
    with uc_file.open() as fin, dest.open("w") as fout:
        for line in fin:
            cols = line.strip().split("\t")
            if cols[0] in {"H", "C"}:
                fout.write(f"{cols[1]}\t{cols[8]}\n")
    print(f"Mapping file written to: {dest}")


def process_phylum(input_dir: Path, output_dir: Path, identity: float, run_vsearch_flag: bool):
    """
    Process one phylum folder and generate OTU outputs.
    Cleans the alignment, runs VSEARCH if enabled, and writes mapping files.
    """

    in_fasta = input_dir / "aligned_sequences_mafft.fasta"

    if not in_fasta.exists():
        print(f"No FASTA file found in {input_dir}. Skipping.")
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    cleaned_fasta = output_dir / "aligned_sequences_cleaned.fasta"
    otus_fasta = output_dir / "otus.fasta"
    otus_uc = output_dir / "otus.uc"
    mapping_txt = output_dir / "otus_mapping.txt"

    _clean_fasta(in_fasta, cleaned_fasta)

    if run_vsearch_flag:
        try:
            _run_vsearch(cleaned_fasta, otus_fasta, otus_uc, identity)
            _build_mapping(otus_uc, mapping_txt)
        except RuntimeError as err:
            print(f"VSEARCH failed in {output_dir}: {err}")
    else:
        print("RUN_VSEARCH = False, clustering skipped.")

    print(f"Processing completed for: {output_dir}\n")
