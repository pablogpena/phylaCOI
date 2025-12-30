# OTU Metrics (diversity, abundance, connections)

## Description
`get_div_abun_conn.R` computes nucleotide diversity, abundance-distance, and haplotype-network
connectivity metrics for each OTU and geographic cluster. It expects abundance outputs and OTU
outputs from the previous steps.

## Input Structure
```
abundance/
└── Phylum1/
    ├── abundances.csv
    └── aligned_sequences_mafft.fasta

otus/
└── Phylum1/
    ├── otus_mapping.txt
    └── informative_OTUs.txt
```

## Usage
```bash
Rscript scripts/4_otu_metrics/get_div_abun_conn.R \
  -a data/abundance \
  -o data/otus
```

### Arguments

| Argument | Type | Required | Description |
|-----------|------|-----------|-------------|
| `-a`, `--abundance` | Path | Yes | Path to the root directory containing per-phylum abundance outputs. |
| `-o`, `--otus` | Path | Yes | Path to the root directory containing per-phylum OTU outputs. |
| `--no-maps` | Flag | No | Skip HTML map generation. |
| `--max-otu-seqs` | Integer | No (default: 500) | Maximum sequences per OTU used for haplotype assignment. |
| `--cluster-radius-km` | Float | No (default: 5) | Spatial clustering radius in km. |

## Outputs
For each phylum, results are written to `otus/<phylum>/informative_otus_metrics/`:
- `otu_metrics_diversity.csv`
- `diversity_vs_distance.csv`
- `otu_metrics_abundance.csv`
- `abundance_vs_distance.csv`
- `otu_metrics_connections.csv`
- `div_abun_conn_combined.csv`
- `otu_nucleotide_diversity_map.html`
- `otu_log_abundance_map.html`
- `otu_log_connections_map.html`

Additionally, a combined table across all phyla is written to:
- `otus/div_abun_conn_master.csv` (includes a `phylum` column)

### Example Output Structure
```
otus/
├── div_abun_conn_master.csv
└── Phylum1/
    └── informative_otus_metrics/
        ├── otu_metrics_diversity.csv
        ├── diversity_vs_distance.csv
        ├── otu_metrics_abundance.csv
        ├── abundance_vs_distance.csv
        ├── otu_metrics_connections.csv
        ├── div_abun_conn_combined.csv
        ├── otu_nucleotide_diversity_map.html
        ├── otu_log_abundance_map.html
        └── otu_log_connections_map.html
```

## Notes
- The script prints a simple progress bar while processing phyla.

## Requirements
- R >= 4.0
- R packages: Biostrings, ape, geosphere, pegas, mgcv, dplyr, tidyr, leaflet, scales, htmlwidgets, tibble
