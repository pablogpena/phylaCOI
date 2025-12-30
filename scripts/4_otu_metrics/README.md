otu_metrics_diversity.csv
diversity_vs_distance.csv
otu_metrics_abundance.csv
abundance_vs_distance.csv
otu_metrics_connections.csv
div_abun_conn_combined.csv
otu_nucleotide_diversity_map.html
otu_log_abundance_map.html
otu_log_connections_map.html# OTU Metrics (diversity, abundance, connections)

## Description
`get_div_abun_conn.R` computes nucleotide diversity, abundance-distance, and haplotype-network
connectivity metrics for each OTU and geographic cluster. It expects the same phylum directory
structure produced by the earlier pipeline steps.

## Input Structure
```
root/
├── Phylum1/
│   └── output/
│       ├── abundances.csv
│       ├── aligned_sequences_mafft.fasta
│       └── otus/
│           ├── otus_mapping.txt
│           └── informative_OTUs.txt
└── Phylum2/
    └── output/
        └── ...
```

## Usage
```bash
Rscript scripts/4_otu_metrics/get_div_abun_conn.R -i /workspace/PhylaCOI/data/vsearch_results
```

### Optional Flags
- `--no-maps` Skip HTML map generation.
- `--max-otu-seqs` Maximum sequences per OTU for haplotype assignment (default: 500).
- `--cluster-radius-km` Spatial clustering radius in km (default: 5).

## Outputs
For each phylum, results are written to `output/informative_otus_metrics/`:
- `otu_metrics_diversity.csv`
- `diversity_vs_distance.csv`
- `otu_metrics_abundance.csv`
- `abundance_vs_distance.csv`
- `otu_metrics_connections.csv`
- `div_abun_conn_combined.csv`
- `otu_nucleotide_diversity_map.html`
- `otu_log_abundance_map.html`
- `otu_log_connections_map.html`

## Requirements
- R >= 4.0
- R packages: Biostrings, ape, geosphere, pegas, mgcv, dplyr, tidyr, leaflet, scales, htmlwidgets, tibble
