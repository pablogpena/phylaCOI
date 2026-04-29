#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# Run OTU metrics for one phylum.

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)])
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(file_arg[1])))
  }
  getwd()
}

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
phylum_name <- get_flag_value(args, c("-p", "--phylum"))
abundance_root <- get_flag_value(args, c("-a", "--abundance-root"))
otus_root <- get_flag_value(args, c("-o", "--otus-root"))
max_otu_sequences <- get_flag_value(args, "--max-otu-seqs", default = 500, convert = as.integer)
cluster_radius_km <- get_flag_value(args, "--cluster-radius-km", default = 5, convert = as.numeric)
write_maps <- !("--no-maps" %in% args)

if (is.null(phylum_name) || is.null(abundance_root) || is.null(otus_root)) {
  stop("Usage: Rscript run_otu_metrics_one_phylum.R --phylum NAME --abundance-root DIR --otus-root DIR")
}
if (is.na(max_otu_sequences) || max_otu_sequences < 1) {
  stop("--max-otu-seqs must be a positive integer.")
}
if (is.na(cluster_radius_km) || cluster_radius_km <= 0) {
  stop("--cluster-radius-km must be a positive number.")
}

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
source(file.path(project_root, "src", "otu_utils", "div_abun_conn_utils.R"))

abundance_root <- normalizePath(abundance_root, mustWork = TRUE)
otus_root <- normalizePath(otus_root, mustWork = TRUE)

process_phylum(
  phylum_name,
  abundance_root,
  otus_root,
  cluster_radius_m = cluster_radius_km * 1000,
  max_sequences_per_otu = max_otu_sequences,
  write_maps = write_maps
)

message("OTU metrics completed for: ", phylum_name)
