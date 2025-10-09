---
# get_abundance.py

## Description
`get_abundance.py` iterates through subfolders representing different phyla and runs the abundance and alignment pipeline defined in `fasta_processing.py`.  
For each phylum, it identifies the FASTA file, mapping file, and metadata file, generates a complete abundance table with read counts and coordinates, creates a FASTA file with unique sequence identifiers, and optionally performs sequence alignment using MAFFT.


## Input Requirements
The script assumes the following directory structure:

```
vsearch_results/
├── *_results.txt # Raw VSEARCH output in BLAST6 format
├── *_results.xlsx # Full results in Excel format
├── *_results_filtered.xlsx # Filtered results by identity threshold
│
├── Phylum1/
│ └──  Phylum1_eKOI_metabarcoding_cleaned.fasta
│
├── Phylum2/
│ └── Phylum2_eKOI_metabarcoding_cleaned.fasta
│
└── ...
```
Each phylum folder must contain **one FASTA file**.  
The mapping and metadata files are provided once for all analyses.

---

## Usage

### Basic Command

```bash
python scripts/get_abundance/get_abundance.py \
  -f /workspace/PhylaCOI/data//vsearch_results/ \
  -n /workspace/PhylaCOI/data/raw/seq_headers.txt \
  -m /workspace/PhylaCOI/data/raw/KOI_metadata.csv \
  -a 1
```
### Arguments

  | Argument | Type | Required | Description |
|-----------|------|-----------|-------------|
| `-f`, `--folder` | Path | Yes | Path to the root folder containing the phylum subfolders (e.g. `/data/vsearch_results/`). |
| `-n`, `--names` | Path | Yes | Full path to the sequence mapping file (e.g. `/data/raw/seq_headers.txt`). |
| `-m`, `--metadata` | Path | Yes | Full path to the metadata file (e.g. `/data/raw/KOI_metadata.csv`). |
| `-a`, `--mafft` | Integer (0 or 1) | Yes | Whether to run MAFFT alignment (`1`) or skip it (`0`). |


## Output

For each phylum subfolder, a new folder `/output` will be created containing the following files:

| File | Description |
|-------|-------------|
| `abundances.csv` | Table with sequence IDs, localities, abundances, coordinates, and unique identifiers. |
| `unique_sequences.fasta` | FASTA file containing non-duplicated sequence IDs. |
| `aligned_sequences_mafft.fasta` | *(optional)* Alignment produced by MAFFT, only if `-a 1`. |

### Example Output Structure
```
Phylum1/
├── Phylum1.fasta
└── output/
    ├── abundances.csv
    ├── unique_sequences.fasta
    └── aligned_sequences_mafft.fasta
 ```
 ## Requirements

To run `get_abundance.py`, the following dependencies are required:

- **Python ≥ 3.8**
- **MAFFT** (must be installed and accessible in the system PATH)

```bash

```
