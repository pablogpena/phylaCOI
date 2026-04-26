# -*- coding: utf-8 -*-

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

has_flag <- function(args, flags) {
  any(args %in% flags)
}

parse_numeric_list <- function(raw_value) {
  values <- strsplit(raw_value, ",", fixed = TRUE)[[1]]
  values <- trimws(values[nzchar(trimws(values))])
  parsed_values <- as.numeric(values)

  if (length(parsed_values) == 0 || any(is.na(parsed_values)) || any(parsed_values <= 0)) {
    stop("Expected a comma-separated list of positive numeric values.")
  }

  parsed_values
}

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."))
source(file.path(project_root, "src", "clustering_utils", "haplotype_clustering_utils.R"))

args <- commandArgs(trailingOnly = TRUE)

if (has_flag(args, c("-h", "--help"))) {
  cat(paste(
    "Usage: Rscript run_same_vs_diff_currents.R",
    "[-i <otus_root>] [-o <output_dir>] [--all-output <all_output_dir>]",
    "[--label LABEL] [--metadata <metadata_csv>] [--sigma-grid V1,V2,V3]",
    "[--k-moran N] [--seed N] [--no-current-figure]\n"
  ))
  quit(status = 0)
}

input_root <- get_flag_value(
  args,
  c("-i", "--input"),
  default = file.path(project_root, "data", "otus")
)
output_dir <- get_flag_value(
  args,
  c("-o", "--output"),
  default = file.path(project_root, "data", "analysis", "haplotype_clustering", "same_vs_diff_currents")
)
all_output_dir <- get_flag_value(
  args,
  c("--all-output"),
  default = file.path(dirname(output_dir), "all_haplotypes")
)
label <- get_flag_value(args, c("--label"), default = "metazoa")
metadata_file <- get_flag_value(args, c("--metadata"))
sigma_grid <- get_flag_value(
  args,
  c("--sigma-grid"),
  default = c(0.005, 0.010, 0.015),
  convert = parse_numeric_list
)
k_moran <- get_flag_value(args, c("--k-moran"), default = 5, convert = as.integer)
seed <- get_flag_value(args, c("--seed"), default = 42, convert = as.integer)
write_current_figure <- !has_flag(args, c("--no-current-figure"))

if (is.na(k_moran) || k_moran < 1) {
  stop("--k-moran must be a positive integer.")
}
if (is.na(seed)) {
  stop("--seed must be an integer.")
}

input_root <- normalizePath(input_root, mustWork = TRUE)
output_dir <- normalizePath(output_dir, mustWork = FALSE)
all_output_dir <- normalizePath(all_output_dir, mustWork = FALSE)
metadata_file <- if (!is.null(metadata_file)) normalizePath(metadata_file, mustWork = TRUE) else NULL

run_same_diff_current_analysis(
  otus_root = input_root,
  output_dir = output_dir,
  label = label,
  all_output_dir = all_output_dir,
  metadata_file = metadata_file,
  sigma_grid = sigma_grid,
  k_moran = k_moran,
  seed = seed,
  write_current_figure = write_current_figure
)

message("SAME vs DIFF current analysis completed successfully.")
