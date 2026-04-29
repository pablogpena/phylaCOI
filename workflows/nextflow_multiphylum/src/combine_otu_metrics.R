#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# Combine per-phylum OTU metrics into the global master table.

suppressPackageStartupMessages({
  library(dplyr)
})

get_flag_value <- function(args, flags, default = NULL, convert = identity) {
  idx <- which(args %in% flags)
  if (length(idx) == 0) {
    return(default)
  }
  if (idx[1] == length(args)) {
    stop("Missing value for ", flags[1])
  }
  convert(args[idx[1] + 1])
}

args <- commandArgs(trailingOnly = TRUE)
otus_root <- get_flag_value(args, c("-o", "--otus-root"))

if (is.null(otus_root)) {
  stop("Usage: Rscript combine_otu_metrics.R --otus-root DIR")
}

otus_root <- normalizePath(otus_root, mustWork = TRUE)
phylum_dirs <- list.dirs(otus_root, recursive = FALSE, full.names = TRUE)

combined_rows <- list()
for (phylum_path in phylum_dirs) {
  phylum_name <- basename(phylum_path)
  combined_path <- file.path(
    phylum_path,
    "informative_otus_metrics",
    "div_abun_conn_combined.csv"
  )

  if (file.exists(combined_path)) {
    df <- read.csv(combined_path, stringsAsFactors = FALSE)
    df$phylum <- phylum_name
    combined_rows[[length(combined_rows) + 1]] <- df
  } else {
    message("Combined metrics file not found for ", phylum_name, ": ", combined_path)
  }
}

if (length(combined_rows) == 0) {
  stop("No per-phylum combined metrics files were found under: ", otus_root)
}

combined_df <- dplyr::bind_rows(combined_rows)
output_file <- file.path(otus_root, "div_abun_conn_master.csv")
write.csv(combined_df, output_file, row.names = FALSE)

message("Master OTU metrics table written to: ", output_file)
