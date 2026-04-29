#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# Run informative OTU selection for one phylum.

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

if (is.null(phylum_name) || is.null(abundance_root) || is.null(otus_root)) {
  stop("Usage: Rscript run_informative_otus_one_phylum.R --phylum NAME --abundance-root DIR --otus-root DIR")
}

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
source(file.path(project_root, "src", "otu_utils", "informative_otus_utils.R"))

abundance_root <- normalizePath(abundance_root, mustWork = TRUE)
otus_root <- normalizePath(otus_root, mustWork = TRUE)

process_phylum(phylum_name, abundance_root, otus_root)

message("Informative OTU selection completed for: ", phylum_name)
