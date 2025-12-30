# -*- coding: utf-8 -*-
# Driver script for diversity, abundance, and connectivity metrics per OTU.

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)])
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(file_arg[1])))
  }
  return(getwd())
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

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."))
source(file.path(project_root, "src", "otu_utils", "div_abun_conn_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
input_dir <- get_flag_value(args, c("-i", "--input"))

if (is.null(input_dir)) {
  stop("Usage: Rscript get_div_abun_conn.R -i <phylum_root> [--no-maps] [--max-otu-seqs N] [--cluster-radius-km K]")
}

base_dir <- normalizePath(input_dir, mustWork = TRUE)
write_maps <- !("--no-maps" %in% args)
max_otu_sequences <- get_flag_value(args, "--max-otu-seqs", default = 500, convert = as.integer)
cluster_radius_km <- get_flag_value(args, "--cluster-radius-km", default = 5, convert = as.numeric)
cluster_radius_m <- cluster_radius_km * 1000

phyla_dirs <- list.dirs(base_dir, recursive = FALSE)
if (length(phyla_dirs) == 0) {
  stop("No phylum folders found in: ", base_dir)
}

for (phylum_path in phyla_dirs) {
  tryCatch({
    process_phylum(
      phylum_path,
      cluster_radius_m = cluster_radius_m,
      max_sequences_per_otu = max_otu_sequences,
      write_maps = write_maps
    )
  }, error = function(e) {
    message("Error processing ", basename(phylum_path), ": ", e$message)
  })
}

message("Pipeline completed successfully.")
