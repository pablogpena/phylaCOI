#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

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

args <- commandArgs(trailingOnly = TRUE)
abundance_root <- get_flag_value(args, c("-a", "--abundance"))
otus_root <- get_flag_value(args, c("-o", "--otus"))

if (is.null(abundance_root) || is.null(otus_root)) {
  stop("Usage: Rscript get_informative_otus.R -a <abundance_root> -o <otus_root>")
}

abundance_root <- normalizePath(abundance_root, mustWork = TRUE)
otus_root <- normalizePath(otus_root, mustWork = TRUE)

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."))
source(file.path(project_root, "src", "otu_utils", "informative_otus_utils.R"))

phyla_dirs <- list.dirs(otus_root, recursive = FALSE)
if (length(phyla_dirs) == 0) {
  stop("No phylum folders were found in the specified OTU directory: ", otus_root)
}

for (phylum_path in phyla_dirs) {
  phylum_name <- basename(phylum_path)
  tryCatch({
    process_phylum(phylum_name, abundance_root, otus_root)
  }, error = function(e) {
    message("Error processing ", phylum_name, ": ", e$message)
  })
}

message("Pipeline completed successfully for all phyla.")
