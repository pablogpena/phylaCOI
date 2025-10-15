<!-- #raw -->
# generate_otus.py

## Description
`generate_otus.py` iterates through subfolders within a root directory and generates Operational Taxonomic Units (OTUs) using [VSEARCH](https://github.com/torognes/vsearch).  
For each folder containing an `output` subdirectory, the script removes gaps from aligned FASTA files, clusters sequences into OTUs based on a user-defined identity threshold, and produces the corresponding centroid and mapping files.

---

## Input Requirements
The script assumes the following directory structure:

```
root_directory/
├── Phylum1/
│   └── output/
│       └── aligned_sequences_mafft.fasta
├── Phylum2/
│   └── output/
│       └── aligned_sequences_mafft.fasta
└── ...
```

Each `output` folder must contain one aligned FASTA file named **`aligned_sequences_mafft.fasta`**.  
The script will automatically create an additional subfolder `/otus/` inside each `output` directory to store the generated results.
<!-- #endraw -->

---

## Usage

### Basic Command

```bash
python scripts/generate_otus/generate_otus.py   -r /workspace/PhylaCOI/data/aligned_fastas/   -i 0.97
```

To skip VSEARCH execution and perform only FASTA cleaning:

```bash
python scripts/generate_otus/generate_otus.py   -r /workspace/PhylaCOI/data/aligned_fastas/   --no-vsearch
```

---

### Arguments

| Argument | Type | Required | Description |
|-----------|------|-----------|-------------|
| `-r`, `--root` | Path | Yes | Path to the root directory containing subfolders with aligned FASTA files. |
| `-i`, `--identity` | Float | No (default: 0.97) | Sequence identity threshold used for clustering with VSEARCH. |
| `--no-vsearch` | Flag | No | If set, skips the clustering step and only performs FASTA cleaning. |

---

## Output

For each processed phylum folder, an `otus/` subdirectory will be created within its `output/` directory containing:

| File | Description |
|-------|-------------|
| `aligned_sequences_cleaned.fasta` | Cleaned version of the original aligned FASTA file (gaps removed). |
| `otus.fasta` | OTU centroid sequences produced by VSEARCH. |
| `otus.uc` | VSEARCH clustering output (cluster membership information). |
| `otus_mapping.txt` | Tab-separated file mapping individual sequences to OTUs. |

### Example Output Structure
```
Phylum1/
└── output/
    ├── aligned_sequences_mafft.fasta
    └── otus/
        ├── aligned_sequences_cleaned.fasta
        ├── otus.fasta
        ├── otus.uc
        └── otus_mapping.txt
```

---

## Requirements

To run `generate_otus.py`, the following dependencies are required:

- **Python ≥ 3.8**
- **VSEARCH ≥ 2.21.0** (must be installed and accessible from the system PATH)

Python standard library modules used:
- `argparse`
- `pathlib`
- `subprocess`

---

## Notes
- The script automatically scans all subfolders under the specified root directory that contain an `output` subdirectory.  
- The expected input file name is fixed as `aligned_sequences_mafft.fasta`.  
- Ensure VSEARCH is correctly installed and callable via the command `vsearch`.

---

## License
This script is provided for research and academic use.  
Please cite appropriately if used in a publication.
