<!-- #raw -->
# generate_otus.py
## Description
`generate_otus.py` iterates through subfolders within a root directory and generates Operational Taxonomic Units (OTUs) using [VSEARCH](https://github.com/torognes/vsearch).  
For each folder containing an `output` subdirectory, the script removes gaps from aligned FASTA files, clusters sequences into OTUs based on a user-defined identity threshold, and produces the corresponding centroid and mapping files.

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


## Usage

### Basic Command

```bash
python scripts/generate_otus/generate_otus.py   -r /workspace/PhylaCOI/data/vsearch_results/   -i 0.97
```

To skip VSEARCH execution and perform only FASTA cleaning:

```bash
python scripts/generate_otus/generate_otus.py   -r /workspace/PhylaCOI/data/vsearch_results/   --no-vsearch
```


### Arguments

| Argument | Type | Required | Description |
|-----------|------|-----------|-------------|
| `-r`, `--root` | Path | Yes | Path to the root directory containing subfolders with aligned FASTA files. |
| `-i`, `--identity` | Float | No (default: 0.97) | Sequence identity threshold used for clustering with VSEARCH. |
| `--no-vsearch` | Flag | No | If set, skips the clustering step and only performs FASTA cleaning. |


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


## Requirements

To run `generate_otus.py`, the following dependencies are required:

- **Python ≥ 3.8**
- **VSEARCH ≥ 2.21.0** (must be installed and accessible from the system PATH)

Python standard library modules used:
- `argparse`
- `pathlib`
- `subprocess`


## Notes
- The script automatically scans all subfolders under the specified root directory that contain an `output` subdirectory.  
- The expected input file name is fixed as `aligned_sequences_mafft.fasta`.  
- Ensure VSEARCH is correctly installed and callable via the command `vsearch`.



# get_informative_otus.R
## Description
`get_informative_otus.R` iterates through phylum subdirectories within a specified base directory, computes genetic and geographic distance matrices, and identifies *informative OTUs* based on distance thresholds.  
This script automates the post-processing step of the OTU pipeline, linking sequence data with spatial coordinates and filtering for meaningful genetic variation.

## Input Requirements
The script assumes the following directory structure for each phylum:

```
Phylum1/
└── output/
    ├── abundances.csv
    ├── aligned_sequences_mafft.fasta
    └── otus/
        ├── otus_mapping.txt
        ├── otus.fasta
        ├── otus.uc
        └── aligned_sequences_cleaned.fasta
```

Required input files:
- **`abundances.csv`** – Table containing sample IDs, sequence abundance, latitude, and longitude.
- **`aligned_sequences_mafft.fasta`** – Aligned nucleotide sequences used to compute pairwise genetic distances.
- **`otus_mapping.txt`** – Two-column tab-separated file mapping unique sequence IDs to OTUs (output of `generate_otus.py`).

Each phylum directory must include an `output` subfolder with the files above.

## Usage

### Basic Command

```bash
Rscript get_informative_otus.R -i /workspace/PhylaCOI/data/vsearch_results/
```

The `-i` (or `--input`) flag specifies the base directory containing the phylum subfolders to process.

### Arguments

| Argument | Type | Required | Description |
|-----------|------|-----------|-------------|
| `-i`, `--input` | Path | Yes | Path to the base directory containing phylum subfolders with `output` data. |

## Output

For each phylum, the script produces an output file inside its corresponding OTU folder:

| File | Description |
|-------|-------------|
| `informative_OTUs.txt` | List of OTUs with ≥ 0.01 genetic distance and ≥ 1 m geographic distance, considered informative for downstream analyses. |

### Example Output Structure
```
Phylum1/
└── output/
    ├── abundances.csv
    ├── aligned_sequences_mafft.fasta
    └── otus/
        ├── otus_mapping.txt
        ├── informative_OTUs.txt
        ├── otus.fasta
        ├── otus.uc
        └── aligned_sequences_cleaned.fasta
```

## Requirements

To run `run_distance_pipeline.R`, the following dependencies are required:

- **R ≥ 4.0**
- R packages:
  - `Biostrings`
  - `ape`
  - `geosphere`
  - `dplyr`

## Notes
- The script automatically scans all phylum directories under the specified base path.  
- Only OTUs with at least four unique sequences are evaluated.  
- Informative OTUs are defined by both minimum genetic (≥ 0.01) and geographic (≥ 1 m) distance thresholds.  
- This script complements `generate_otus.py` by identifying the most relevant OTUs for diversity and connectivity analyses.
<!-- #endraw -->