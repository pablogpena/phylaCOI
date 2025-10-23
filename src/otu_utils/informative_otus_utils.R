# -*- coding: utf-8 -*-
suppressPackageStartupMessages({
  library(Biostrings)
  library(ape)
  library(geosphere)
  library(dplyr)
})

# Main function to process a single phylum
process_phylum <- function(phylum_path) {
  phylum_name <- basename(phylum_path)
  message("Processing: ", phylum_name)
  
  # Define input paths
  abundances_file    <- file.path(phylum_path, "output", "abundances.csv")
  aligned_fasta_file <- file.path(phylum_path, "output", "aligned_sequences_mafft.fasta")
  otus_file          <- file.path(phylum_path, "output", "otus", "otus_mapping.txt")
  
  # Check file existence
  if (!file.exists(abundances_file) || !file.exists(aligned_fasta_file) || !file.exists(otus_file)) {
    warning("Missing files in: ", phylum_name)
    return(NULL)
  }  

  # Load data
  abundances <- read.csv(abundances_file)
  otus_map   <- read.table(otus_file, col.names = c("OTU", "UniqueID"))
  abundances$OTU <- otus_map$OTU[match(abundances$UniqueID, otus_map$UniqueID)]
  
  alignment <- readDNAStringSet(aligned_fasta_file)
  alignment_dnabin <- as.DNAbin(alignment)
  genetic_dist_mat <- as.matrix(dist.dna(alignment_dnabin, model = "K80"))
  

  # Filter OTUs with ≥ 4 unique sequences
  abundances_no_dup <- abundances[!duplicated(abundances$ID), ]
  otu_counts <- table(abundances_no_dup$OTU)
  valid_otus <- names(otu_counts[otu_counts >= 4])
  abundances_filt <- abundances[abundances$OTU %in% valid_otus, ]
  
  if (length(valid_otus) == 0) {
    message("No valid OTUs found in: ", phylum_name)
    return(NULL)
  }
  
  # Build geographic distance matrix
  coords <- abundances[, c("UniqueID", "lat", "lon")]
  coords <- coords[match(rownames(genetic_dist_mat), coords$UniqueID), ]
  rownames(coords) <- coords$UniqueID
  coords$UniqueID <- NULL
  
  geo_dist_mat <- distm(coords[, c("lon", "lat")], fun = distHaversine)
  rownames(geo_dist_mat) <- colnames(geo_dist_mat) <- rownames(coords)

  # Combine distances per OTU
  combined_list <- list()
  for (otu in valid_otus) {
    ids <- abundances_filt$UniqueID[abundances_filt$OTU == otu]
    ids <- intersect(ids, rownames(genetic_dist_mat))
    if (length(ids) == 0) next
    
    gmat <- genetic_dist_mat[ids, ids, drop = FALSE]
    dmat <- geo_dist_mat[ids, ids, drop = FALSE]
    
    gdf <- as.data.frame(as.table(gmat), stringsAsFactors = FALSE)
    ddf <- as.data.frame(as.table(dmat), stringsAsFactors = FALSE)
    colnames(gdf) <- c("Item1", "Item2", "Genetic_Distance")
    colnames(ddf) <- c("Item1", "Item2", "Geographic_Distance")
    
    cdf <- merge(gdf, ddf, by = c("Item1", "Item2"))
    if (nrow(cdf) > 0) {
      cdf$OTU <- otu
      combined_list[[otu]] <- cdf
    }
  }
  
  combined_df <- do.call(rbind, combined_list)
  if (is.null(combined_df)) {
    message("No combined data for: ", phylum_name)
    return(NULL)
  }
  
  # Keep informative OTUs with >= 0.01 genetic and >= 1 m geographic distance
  otus_keep <- combined_df |>
    dplyr::group_by(OTU) |>
    dplyr::filter(any(Genetic_Distance >= 0.01) &
                    any(Geographic_Distance >= 1)) |>
    dplyr::pull(OTU) |> unique()
  
  final_df <- combined_df |> dplyr::filter(OTU %in% otus_keep)
  
  # -----------------------------------------------------------------------------#
  # Write results
  # -----------------------------------------------------------------------------#
  otus_keep_file <- file.path(phylum_path, "output", "otus", "informative_OTUs.txt")
  writeLines(otus_keep, otus_keep_file)
  
  message("Informative OTUs saved to: ", otus_keep_file)
