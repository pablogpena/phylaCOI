<!-- #raw -->
# generate_otus.py
## Description
`generate_otus.py` iterates through subfolders within a root directory and generates Operational Taxonomic Units (OTUs) using [VSEARCH](https://github.com/torognes/vsearch).  
For each phylum folder in the abundance output, the script removes gaps from aligned FASTA files, clusters sequences into OTUs based on a user-defined identity threshold, and produces the corresponding centroid and mapping files.

## Input Requirements
The script assumes the following directory structure:

```
abundance/
├── Phylum1/
│   └── aligned_sequences_mafft.fasta
├── Phylum2/
│   └── aligned_sequences_mafft.fasta
└── ...
```

Each phylum folder must contain one aligned FASTA file named **`aligned_sequences_mafft.fasta`**.  
The script will create a matching phylum folder under the OTU output root to store the results.


## Usage

### Basic Command

```bash
python scripts/3_otu_generation/generate_otus.py \
  -a data/abundance \
  -o data/otus \
  -i 0.97
```

To skip VSEARCH execution and perform only FASTA cleaning:

```bash
python scripts/3_otu_generation/generate_otus.py \
  -a data/abundance \
  -o data/otus \
  --no-vsearch
```


### Arguments

| Argument | Type | Required | Description |
|-----------|------|-----------|-------------|
| `-a`, `--abundance` | Path | Yes | Path to the root directory containing per-phylum abundance outputs. |
| `-o`, `--otus` | Path | Yes | Path to the root directory where per-phylum OTU outputs will be written. |
| `-i`, `--identity` | Float | No (default: 0.97) | Sequence identity threshold used for clustering with VSEARCH. |
| `--no-vsearch` | Flag | No | If set, skips the clustering step and only performs FASTA cleaning. |


## Output

For each processed phylum folder, a matching folder will be created within the OTU output root containing:

| File | Description |
|-------|-------------|
| `aligned_sequences_cleaned.fasta` | Cleaned version of the original aligned FASTA file (gaps removed). |
| `otus.fasta` | OTU centroid sequences produced by VSEARCH. |
| `otus.uc` | VSEARCH clustering output (cluster membership information). |
| `otus_mapping.txt` | Tab-separated file mapping individual sequences to OTUs. |

### Example Output Structure
```
otus/
└── Phylum1/
    ├── aligned_sequences_cleaned.fasta
    ├── otus.fasta
    ├── otus.uc
    └── otus_mapping.txt
```
This folder is used by `get_informative_otus.R` and step 4 (OTU metrics).


## Requirements

To run `generate_otus.py`, the following dependencies are required:

- **Python ≥ 3.8**
- **VSEARCH ≥ 2.21.0** (must be installed and accessible from the system PATH)

Python standard library modules used:
- `argparse`
- `pathlib`
- `subprocess`


## Notes
- The script automatically scans all phylum folders under the specified abundance root.  
- The expected input file name is fixed as `aligned_sequences_mafft.fasta`.  
- Ensure VSEARCH is correctly installed and callable via the command `vsearch`.



# get_informative_otus.R
## Description
`get_informative_otus.R` iterates through phylum subdirectories within a specified base directory, computes genetic and geographic distance matrices, and identifies *informative OTUs* based on distance thresholds.  
This script automates the post-processing step of the OTU pipeline, linking sequence data with spatial coordinates and filtering for meaningful genetic variation.

## Input Requirements
The script assumes the following directory structure for each phylum:

```
abundance/
└── Phylum1/
    ├── abundances.csv
    └── aligned_sequences_mafft.fasta

otus/
└── Phylum1/
    ├── otus_mapping.txt
    ├── otus.fasta
    ├── otus.uc
    └── aligned_sequences_cleaned.fasta
```

Required input files:
- **`abundances.csv`** – Table containing sample IDs, sequence abundance, latitude, and longitude.
- **`aligned_sequences_mafft.fasta`** – Aligned nucleotide sequences used to compute pairwise genetic distances.
- **`otus_mapping.txt`** – Two-column tab-separated file mapping unique sequence IDs to OTUs (output of `generate_otus.py`).

Each phylum directory must include the files above under the abundance and otus roots.

## Usage

### Basic Command

```bash
Rscript scripts/3_otu_generation/get_informative_otus.R \
  -a data/abundance \
  -o data/otus
```

The `-a` flag points to the abundance root and `-o` to the OTU root.

### Arguments

| Argument | Type | Required | Description |
|-----------|------|-----------|-------------|
| `-a`, `--abundance` | Path | Yes | Path to the base directory containing per-phylum abundance outputs. |
| `-o`, `--otus` | Path | Yes | Path to the base directory containing per-phylum OTU outputs. |

## Output

For each phylum, the script produces an output file inside its corresponding OTU folder:

| File | Description |
|-------|-------------|
| `informative_OTUs.txt` | List of OTUs with ≥ 0.01 genetic distance and ≥ 1 m geographic distance, considered informative for downstream analyses. |

### Example Output Structure
```
otus/
└── Phylum1/
    ├── otus_mapping.txt
    ├── informative_OTUs.txt
    ├── otus.fasta
    ├── otus.uc
    └── aligned_sequences_cleaned.fasta
```

## Requirements

To run `get_informative_otus.R`, the following dependencies are required:

- **R ≥ 4.0**
- R packages:
  - `Biostrings`
  - `ape`
  - `geosphere`
  - `dplyr`

## Notes
- The script automatically scans all phylum directories under the OTU root.  
- Only OTUs with at least four unique sequences are evaluated.  
- Informative OTUs are defined by both minimum genetic (≥ 0.01) and geographic (≥ 1 m) distance thresholds.  
- This script complements `generate_otus.py` by identifying the most relevant OTUs for diversity and connectivity analyses.
<!-- #endraw -->
