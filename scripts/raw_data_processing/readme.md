# Raw Data Processing

This folder contains two scripts:

- `fasta_preprocess.py`
- `vsearch_taxonomy.py`

They are intended to help with FASTA preprocessing and taxonomy assignment using **vsearch**.

---

## 1) `fasta_preprocess.py`

### Arguments
- `-i / --input` → Input FASTA file  
- `-t / --txt` → Output TXT file with sequence headers  
- `-o / --output` → Output FASTA file with cleaned identifiers  

### Example
```bash
python fasta_preprocess.py \
  -i /data/raw/eKOI_metabarcoding.fasta \
  -t /data/raw/seq_headers.txt \
  -o /data/procesed/eKOI_metabarcoding_cleaned.fasta
```  
###  Outputs

- `seq_headers.txt` → all sequence headers without >
- `eKOI_metabarcoding_cleaned.fasta` → FASTA file with simplified identifiers


## 2) `vsearch_taxonomy.py`

### Purpose
Classify FASTA sequences against a reference database using **vsearch**, filter results, and group sequences by taxonomy.

### Workflow
- Run `vsearch --usearch_global` for each FASTA file in a folder (except the reference DB).  
- Save raw results in TXT and Excel files.  
- Filter hits by identity threshold.  
- Split taxonomy strings (by `;`) and group sequences by the 5th rank.  
- Create group-specific FASTA files.  
- Save filtered results in Excel.  

### Usage
```bash
python vsearch_taxonomy.py -d reference.fasta -i 0.84 -f /path/to/fasta_folder -o /path/to/output_folder
```
### Arguments
- `-d / --db` → Reference database FASTA file  
- `-i / --identity` → Minimum identity threshold (e.g., 0.84)  
- `-f / --folder` → Folder containing FASTA files to process  
- `-o / --output` → Optional base folder for results (default: alongside each FASTA)  

### Examples

**Results created next to each FASTA file**
```bash
python vsearch_taxonomy.py -d /data/raw/eKOI_database.fasta -i 0.84 -f /data/procesed/eKOI_metabarcoding_cleaned.fasta
```
**Results centralized in one folder**
```bash
python vsearch_taxonomy.py -d /data/raw/eKOI_database.fasta -i 0.84 -f /data/procesed/eKOI_metabarcoding_cleaned.fasta -o /data/vsearch_results/
```

### Outputs per FASTA file
- `*_results.txt` → Raw vsearch output (BLAST6 format)  
- `*_results.xlsx` → Results in Excel format  
- `*_results_filtered.xlsx` → Filtered results (≥ identity threshold)  
- `{group_name}_*.fasta` → Group-specific FASTA files by taxonomic rank 5  

---

## Requirements

- Python 3.8+  
- [Biopython](https://biopython.org/)  
- [pandas](https://pandas.pydata.org/)  
- [vsearch](https://github.com/torognes/vsearch) installed and available in `$PATH`
