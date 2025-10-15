# -*- coding: utf-8 -*-
from pathlib import Path
import subprocess


def _clean_fasta(src: Path, dest: Path):
    """
    Remove gap characters ('-') from an aligned FASTA file and save the cleaned output.
    """
    with src.open() as fin, dest.open("w") as fout:
        for line in fin:
            fout.write(line if line.startswith(">") else line.replace("-", ""))
    print(f"Cleaned FASTA written to: {dest}")


def _run_vsearch(fasta: Path, centroids: Path, uc_file: Path, identity: float):
    """
    Execute VSEARCH to cluster sequences into OTUs.
    
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
    Generate a mapping file linking sequence IDs to their corresponding OTUs
    based on a VSEARCH .uc file.

    """
    with uc_file.open() as fin, dest.open("w") as fout:
        for line in fin:
            cols = line.strip().split("\t")
            if cols[0] in {"H", "C"}:
                fout.write(f"{cols[1]}\t{cols[8]}\n")
    print(f"Mapping file written to: {dest}")


def process_folder(folder: Path, identity: float, run_vsearch_flag: bool):
    """
    Process a folder containing aligned sequences and generate OTUs.
    
    Parameters
    ----------
    folder : Path
        Directory containing an 'output' subfolder with aligned sequences.
    identity : float
        Sequence identity threshold for clustering.
    run_vsearch_flag : bool
        Whether to execute VSEARCH (True) or skip clustering (False).
    """

    output_folder = folder / "output"
    in_fasta = output_folder / "aligned_sequences_mafft.fasta"

    if not in_fasta.exists():
        print(f"No FASTA file found in {output_folder}. Skipping.")
        return

    out_dir = output_folder / "otus"
    out_dir.mkdir(exist_ok=True)

    cleaned_fasta = out_dir / "aligned_sequences_mafft_cleaned.fasta"
    otus_fasta = out_dir / "otus.fasta"
    otus_uc = out_dir / "otus.uc"
    mapping_txt = out_dir / "otus_mapping.txt"

    _clean_fasta(in_fasta, cleaned_fasta)

    if run_vsearch_flag:
        try:
            _run_vsearch(cleaned_fasta, otus_fasta, otus_uc, identity)
            _build_mapping(otus_uc, mapping_txt)
        except RuntimeError as err:
            print(f"VSEARCH failed in {output_folder}: {err}")
    else:
        print("RUN_VSEARCH = False, clustering skipped.")

    print(f"Processing completed for: {output_folder}\n")
