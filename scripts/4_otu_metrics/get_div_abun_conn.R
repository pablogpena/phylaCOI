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
abundance_root <- get_flag_value(args, c("-a", "--abundance"))
otus_root <- get_flag_value(args, c("-o", "--otus"))

if (is.null(abundance_root) || is.null(otus_root)) {
  stop("Usage: Rscript get_div_abun_conn.R -a <abundance_root> -o <otus_root> [--no-maps] [--max-otu-seqs N] [--cluster-radius-km K]")
}

abundance_root <- normalizePath(abundance_root, mustWork = TRUE)
otus_root <- normalizePath(otus_root, mustWork = TRUE)
write_maps <- !("--no-maps" %in% args)
max_otu_sequences <- get_flag_value(args, "--max-otu-seqs", default = 500, convert = as.integer)
cluster_radius_km <- get_flag_value(args, "--cluster-radius-km", default = 5, convert = as.numeric)
cluster_radius_m <- cluster_radius_km * 1000

phyla_dirs <- list.dirs(otus_root, recursive = FALSE)
if (length(phyla_dirs) == 0) {
  stop("No phylum folders found in: ", otus_root)
}

progress <- utils::txtProgressBar(min = 0, max = length(phyla_dirs), style = 3)
for (i in seq_along(phyla_dirs)) {
  phylum_path <- phyla_dirs[[i]]
  phylum_name <- basename(phylum_path)
  tryCatch({
    process_phylum(
      phylum_name,
      abundance_root,
      otus_root,
      cluster_radius_m = cluster_radius_m,
      max_sequences_per_otu = max_otu_sequences,
      write_maps = write_maps
    )
  }, error = function(e) {
    message("Error processing ", phylum_name, ": ", e$message)
  })
  utils::setTxtProgressBar(progress, i)
}
close(progress)

combined_rows <- list()
for (phylum_path in phyla_dirs) {
  phylum_name <- basename(phylum_path)
  combined_path <- file.path(
    otus_root,
    phylum_name,
    "informative_otus_metrics",
    "div_abun_conn_combined.csv"
  )

  if (file.exists(combined_path)) {
    df <- read.csv(combined_path, stringsAsFactors = FALSE)
    df$phylum <- phylum_name
    combined_rows[[length(combined_rows) + 1]] <- df
  } else {
    message("Combined file not found: ", combined_path)
  }
}

if (length(combined_rows) > 0) {
  combined_df <- dplyr::bind_rows(combined_rows)
  write.csv(
    combined_df,
    file = file.path(otus_root, "div_abun_conn_master.csv"),
    row.names = FALSE
  )
} else {
  message("No combined files were found; global CSV not written.")
}

message("Pipeline completed successfully.")
