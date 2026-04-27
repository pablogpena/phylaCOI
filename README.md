# phylaCOI

`phylaCOI` is a reproducible workflow for processing COI metabarcoding sequences, assigning taxonomy, estimating sequence abundance, generating OTUs, computing OTU-level metrics, and running downstream metric and haplotype-clustering analyses.

Each pipeline step has its own README with detailed inputs, arguments, outputs, and requirements. This root README only summarizes the full workflow and the commands needed to run it from the repository root.

## Repository Structure

```
phylaCOI/
├── data/
│   ├── raw/
│   ├── procesed/
│   ├── vsearch_results/
│   ├── abundance/
│   ├── otus/
│   └── analysis/
├── scripts/
│   ├── 1_raw_data_processing/
│   ├── 2_abundance_estimation/
│   ├── 3_otu_generation/
│   ├── 4_otu_metrics/
│   ├── 5_metrics_analysis/
│   └── 6_haplotype_clustering/
└── src/
```

## Workflow Overview

```text
1_raw_data_processing
        ↓
2_abundance_estimation
        ↓
3_otu_generation
        ↓
4_otu_metrics
        ├── 5_metrics_analysis
        └── 6_haplotype_clustering
```

Step 4 is the branching point. It creates the combined OTU metrics table used by step 5 and the haplotype-network tables used by step 6.

## Full Pipeline

Run all commands from the repository root.

### 1. Raw Data Processing

Clean FASTA identifiers and save the original sequence headers:

```bash
python scripts/1_raw_data_processing/fasta_preprocess.py \
  -i data/raw/eKOI_metabarcoding.fasta \
  -t data/raw/seq_headers.txt \
  -o data/procesed/eKOI_metabarcoding_cleaned.fasta
```

Assign taxonomy with VSEARCH and split sequences into phylum folders:

```bash
python scripts/1_raw_data_processing/vsearch_taxonomy.py \
  -d data/raw/eKOI_database.fasta \
  -i 0.84 \
  -f data/procesed/eKOI_metabarcoding_cleaned.fasta \
  -o data/vsearch_results/
```

More details: [scripts/1_raw_data_processing/readme.md](scripts/1_raw_data_processing/readme.md)

### 2. Abundance Estimation

Create per-phylum abundance tables, unique-sequence FASTA files, and MAFFT alignments:

```bash
python scripts/2_abundance_estimation/get_abundance.py \
  -i data/vsearch_results \
  -o data/abundance \
  -n data/raw/seq_headers.txt \
  -m data/raw/KOI_metadata.csv \
  -a 1
```

More details: [scripts/2_abundance_estimation/readme.md](scripts/2_abundance_estimation/readme.md)

### 3. OTU Generation

Cluster aligned sequences into OTUs:

```bash
python scripts/3_otu_generation/generate_otus.py \
  -a data/abundance \
  -o data/otus \
  -i 0.97
```

Identify informative OTUs for downstream analyses:

```bash
Rscript scripts/3_otu_generation/get_informative_otus.R \
  -a data/abundance \
  -o data/otus
```

More details: [scripts/3_otu_generation/README.md](scripts/3_otu_generation/README.md)

### 4. OTU Metrics

Compute diversity, abundance, and haplotype-connection metrics for informative OTUs:

```bash
Rscript scripts/4_otu_metrics/get_div_abun_conn.R \
  -a data/abundance \
  -o data/otus
```

This step writes:

- per-phylum OTU metrics to `data/otus/<phylum>/informative_otus_metrics/`
- per-phylum haplotype-network inputs to `data/otus/<phylum>/haplotype_network/`
- the combined metrics table to `data/otus/div_abun_conn_master.csv`

More details: [scripts/4_otu_metrics/README.md](scripts/4_otu_metrics/README.md)

### 5. Metrics Analysis

Analyze the combined OTU metrics table across phyla:

```bash
Rscript scripts/5_metrics_analysis/analyze_otu_metrics.R \
  -i data/otus/div_abun_conn_master.csv \
  -o data/analysis/otu_metrics_summary \
  --min-otus 10
```

More details: [scripts/5_metrics_analysis/README.md](scripts/5_metrics_analysis/README.md)

### 6. Haplotype Clustering

Run all-haplotype clustering:

```bash
Rscript scripts/6_haplotype_clustering/run_all_haplotypes_clustering.R \
  -i data/otus \
  -o data/analysis/haplotype_clustering/all_haplotypes \
  --metadata data/raw/metadata_eKOI_ver2.csv
```

Run same-haplotype vs different-haplotype clustering:

```bash
Rscript scripts/6_haplotype_clustering/run_same_vs_diff_currents.R \
  -i data/otus \
  -o data/analysis/haplotype_clustering/same_vs_diff_currents \
  --all-output data/analysis/haplotype_clustering/all_haplotypes
```

More details: [scripts/6_haplotype_clustering/README.md](scripts/6_haplotype_clustering/README.md)

## Main Outputs

The main output roots are:

| Folder | Contents |
|--------|----------|
| `data/vsearch_results/` | VSEARCH taxonomy results and phylum-specific FASTA files. |
| `data/abundance/` | Per-phylum abundance tables, unique-sequence FASTA files, and MAFFT alignments. |
| `data/otus/` | OTU files, informative OTU lists, OTU metrics, haplotype-network tables, and `div_abun_conn_master.csv`. |
| `data/analysis/otu_metrics_summary/` | Global metric models, correlations, predictions, and plots. |
| `data/analysis/haplotype_clustering/` | Locality clustering, current-zone summaries, distance summaries, maps, and figures. |

## Requirements

The workflow uses Python, R, VSEARCH, and MAFFT. Each step README lists the exact packages and command-line tools required for that step.

Recommended reading before running a step:

- [Raw data processing](scripts/1_raw_data_processing/readme.md)
- [Abundance estimation](scripts/2_abundance_estimation/readme.md)
- [OTU generation](scripts/3_otu_generation/README.md)
- [OTU metrics](scripts/4_otu_metrics/README.md)
- [Metrics analysis](scripts/5_metrics_analysis/README.md)
- [Haplotype clustering](scripts/6_haplotype_clustering/README.md)
