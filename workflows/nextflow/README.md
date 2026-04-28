# Nextflow Workflow

This folder contains a phase 1 Nextflow wrapper for the existing `phylaCOI` pipeline.

The workflow does not replace or modify the scripts in `scripts/` or the functions in `src/`. It only orchestrates the current commands in the correct order.

## Files

| File | Description |
|------|-------------|
| `main.nf` | Main Nextflow workflow. |
| `nextflow.config` | Default parameters and local execution settings. |
| `params.yaml` | Example parameter file for custom runs. |

## Workflow

```text
RAW_PREPROCESS
        ↓
VSEARCH_TAXONOMY
        ↓
ABUNDANCE_ESTIMATION
        ↓
OTU_GENERATION
        ↓
INFORMATIVE_OTUS
        ↓
OTU_METRICS
        ├── METRICS_ANALYSIS
        └── HAPLOTYPE_CLUSTERING_ALL
                    ↓
             HAPLOTYPE_CLUSTERING_SAME_DIFF
```

This is a coarse-grained workflow: each process calls one of the current pipeline scripts. Several scripts still process all phyla internally, so this first version is mainly an orchestration and reproducibility layer.

## Input Files

By default, the workflow expects:

| Parameter | Default |
|-----------|---------|
| `raw_fasta` | `data/raw/eKOI_metabarcoding.fasta` |
| `reference_db` | `data/raw/eKOI_database.fasta` |
| `sample_metadata` | `data/raw/KOI_metadata.csv` |
| `ocean_metadata` | `data/raw/ocean_metadata.csv` |
| `output_root` | `data` |

## Basic Run

Run from the repository root:

```bash
nextflow run workflows/nextflow/main.nf \
  -c workflows/nextflow/nextflow.config
```

To resume after a failed or interrupted run:

```bash
nextflow run workflows/nextflow/main.nf \
  -c workflows/nextflow/nextflow.config \
  -resume
```

## Run With Custom Parameters

You can override parameters directly:

```bash
nextflow run workflows/nextflow/main.nf \
  -c workflows/nextflow/nextflow.config \
  --raw_fasta data/raw/my_metabarcoding.fasta \
  --reference_db data/raw/my_reference_db.fasta \
  --sample_metadata data/raw/my_sample_metadata.csv \
  --ocean_metadata data/raw/my_ocean_metadata.csv
```

Or use the example parameter file:

```bash
nextflow run workflows/nextflow/main.nf \
  -c workflows/nextflow/nextflow.config \
  -params-file workflows/nextflow/params.yaml
```

## Outputs

The workflow publishes outputs to the same project structure used by the standalone scripts:

| Folder | Contents |
|--------|----------|
| `data/raw/` | Generated `seq_headers.txt`. |
| `data/procesed/` | Cleaned FASTA file. |
| `data/vsearch_results/` | VSEARCH taxonomy outputs. |
| `data/abundance/` | Per-phylum abundance and alignment outputs. |
| `data/otus/` | OTU outputs, informative OTUs, metrics, haplotype networks, and `div_abun_conn_master.csv`. |
| `data/analysis/otu_metrics_summary/` | Global metric-analysis outputs. |
| `data/analysis/haplotype_clustering/` | Haplotype clustering outputs. |

## Requirements

Nextflow must be installed and available in `PATH`.

The runtime environment must also provide the same tools and packages required by the standalone pipeline:

| Tool | Used by |
|------|---------|
| Python | Steps 1, 2, and 3. |
| R / Rscript | Steps 3, 4, 5, and 6. |
| VSEARCH | Taxonomic assignment and OTU generation. |
| MAFFT | Sequence alignment during abundance estimation. |

See the README inside each `scripts/<step>/` folder for exact package requirements.

## Notes

- This phase 1 workflow intentionally keeps the existing code untouched.
- Processes are coarse-grained and do not parallelize by phylum yet.
- A future phase 2 Nextflow version could split the workflow into one task per phylum.
