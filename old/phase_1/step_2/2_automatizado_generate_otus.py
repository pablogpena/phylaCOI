from pathlib import Path
import subprocess

# Ruta raíz donde están todas las carpetas que contienen "output"
ROOT_DIR = Path("/mnt/e/Gonzalez-Miguens preprint 2025/000000_version2_maria/FASE1/3_magia_nombre_carpeta/eKOI_metabarcoding_database_FILOS")

# Configuración
ALIGNED_FASTA   = "aligned_sequences_mafft.fasta"
OUTPUT_DIR      = "otus"
CLEANED_FASTA   = "aligned_sequences_cleaned.fasta"
OTUS_FASTA      = "otus.fasta"
OTUS_UC         = "otus.uc"
OTU_MAPPING_TXT = "otus_mapping.txt"
IDENTITY        = 0.97
RUN_VSEARCH     = True

def clean_fasta(src: Path, dest: Path):
    with src.open() as fin, dest.open("w") as fout:
        for line in fin:
            fout.write(line if line.startswith(">") else line.replace("-", ""))
    print(f"✓ Clean FASTA → {dest}")

def run_vsearch(fasta: Path, centroids: Path, uc_file: Path, identity: float):
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
    print(f"✓ VSEARCH OK → {centroids}, {uc_file}")

def build_mapping(uc_file: Path, dest: Path):
    with uc_file.open() as fin, dest.open("w") as fout:
        for line in fin:
            cols = line.strip().split("\t")
            if cols[0] in {"H", "C"}:
                fout.write(f"{cols[1]}\t{cols[8]}\n")
    print(f"✓ Mapping → {dest}")

def process_folder(folder: Path):
    output_folder = folder / "output"
    in_fasta = output_folder / ALIGNED_FASTA

    if not in_fasta.exists():
        print(f"✗ No FASTA in {output_folder}, skipped.")
        return

    out_dir = output_folder / OUTPUT_DIR
    out_dir.mkdir(exist_ok=True)

    cleaned_fasta = out_dir / CLEANED_FASTA
    otus_fasta    = out_dir / OTUS_FASTA
    otus_uc       = out_dir / OTUS_UC
    mapping_txt   = out_dir / OTU_MAPPING_TXT

    clean_fasta(in_fasta, cleaned_fasta)

    if RUN_VSEARCH:
        try:
            run_vsearch(cleaned_fasta, otus_fasta, otus_uc, IDENTITY)
            build_mapping(otus_uc, mapping_txt)
        except RuntimeError as err:
            print(f"✗ VSEARCH failed in {output_folder}: {err}")
    else:
        print("• RUN_VSEARCH = False, clustering skipped")

    print(f"✓ Done → {output_folder}\n")

# Iterar sobre cada carpeta que contenga un subdirectorio "output"
for folder in ROOT_DIR.iterdir():
    if folder.is_dir() and (folder / "output").is_dir():
        print(f"▶ Procesando: {folder.name}")
        process_folder(folder)
