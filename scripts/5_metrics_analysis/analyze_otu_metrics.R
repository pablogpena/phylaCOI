# -*- coding: utf-8 -*-

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
source(file.path(project_root, "src", "analysis_utils", "otu_metrics_analysis_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
input_file <- get_flag_value(args, c("-i", "--input"))
output_dir <- get_flag_value(args, c("-o", "--output"))
comparisons_file <- get_flag_value(args, c("--comparisons"))
heatmap_model <- get_flag_value(args, c("--heatmap-model"), default = "nls")
heatmap_metrics_raw <- get_flag_value(args, c("--heatmap-metrics"))
min_otus <- get_flag_value(args, c("--min-otus"), default = 10, convert = as.integer)

if (is.null(input_file) || is.null(output_dir)) {
  stop(paste(
    "Usage: Rscript analyze_otu_metrics.R -i <metrics_csv> -o <output_dir>",
    "[--min-otus N] [--comparisons <comparison_csv>] [--heatmap-model MODEL] [--heatmap-metrics M1,M2]"
  ))
}

input_file <- normalizePath(input_file, mustWork = TRUE)
output_dir <- normalizePath(output_dir, mustWork = FALSE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

metrics_data <- read_metrics_table(input_file)
metrics_data <- add_scaled_metrics(metrics_data)

if (is.na(min_otus) || min_otus < 1) {
  stop("--min-otus must be a positive integer.")
}

# Keep only phyla with at least the requested number of unique OTUs.
otu_counts <- metrics_data %>%
  group_by(phylum) %>%
  summarise(otu_count = n_distinct(OTU), .groups = "drop")

eligible_phyla <- otu_counts %>%
  filter(otu_count >= min_otus) %>%
  pull(phylum)

filtered_data <- metrics_data %>% filter(phylum %in% eligible_phyla)

if (length(eligible_phyla) == 0) {
  stop("No phyla meet the minimum OTU threshold (", min_otus, ").")
}

if (length(eligible_phyla) < length(unique(metrics_data$phylum))) {
  removed <- setdiff(unique(metrics_data$phylum), eligible_phyla)
  message(
    "Excluded phyla with fewer than ",
    min_otus,
    " OTUs: ",
    paste(removed, collapse = ", ")
  )
}

# Model fitting per metric.
res_abundance <- fit_models_for_metric(filtered_data, "scaled_log_abundance", "log_abundance")
res_diversity <- fit_models_for_metric(filtered_data, "scaled_diversity", "diversity")
res_connections <- fit_models_for_metric(filtered_data, "scaled_log_connections", "log_connections")

model_summary <- bind_rows(
  res_abundance$results,
  res_diversity$results,
  res_connections$results
)
model_predictions <- bind_rows(
  res_abundance$predictions,
  res_diversity$predictions,
  res_connections$predictions
)

predictions_long <- build_predictions_long(model_predictions)

write.csv(
  model_summary,
  file = file.path(output_dir, "model_fit_summary.csv"),
  row.names = FALSE
)
write.csv(
  predictions_long,
  file = file.path(output_dir, "model_predictions_long.csv"),
  row.names = FALSE
)

plot_model_comparison(
  predictions_long,
  file.path(output_dir, "model_comparison_by_phylum.pdf")
)

# Spearman correlations and plots.
cor_outputs <- compute_spearman_correlations(filtered_data)
write.csv(
  cor_outputs$cor_total,
  file = file.path(output_dir, "spearman_correlations.csv"),
  row.names = FALSE
)
plot_spearman_bars(cor_outputs$cor_phylum_df, cor_outputs$comparisons, output_dir)
plot_spearman_heatmap(
  cor_outputs$cor_total,
  file.path(output_dir, "spearman_correlation_heatmap.pdf")
)
plot_spearman_network(
  cor_outputs$cor_global_df,
  file.path(output_dir, "spearman_correlation_network.pdf")
)

# Permutation tests for slope differences between metrics.
metric_pairs <- list(
  c("scaled_log_abundance", "scaled_diversity"),
  c("scaled_log_abundance", "scaled_log_connections"),
  c("scaled_diversity", "scaled_log_connections")
)
model_types <- c("lm", "power", "gam", "nls", "gompertz")

perm_results <- run_metric_permutation_tests(filtered_data, metric_pairs, model_types, B = 999)
write.csv(
  perm_results,
  file = file.path(output_dir, "metric_slope_permutation_tests.csv"),
  row.names = FALSE
)

slopes_all <- bind_rows(
  get_all_slopes(filtered_data, "lm"),
  get_all_slopes(filtered_data, "power"),
  get_all_slopes(filtered_data, "gam"),
  get_all_slopes(filtered_data, "nls"),
  get_all_slopes(filtered_data, "gompertz")
)
write.csv(
  slopes_all,
  file = file.path(output_dir, "metric_slopes_by_model.csv"),
  row.names = FALSE
)

# Optional phylum comparison heatmaps.
if (!is.null(comparisons_file)) {
  comparisons_file <- normalizePath(comparisons_file, mustWork = TRUE)
  comparisons <- read_phylum_comparisons(comparisons_file)

  if (!is.null(heatmap_metrics_raw)) {
    heatmap_metrics <- strsplit(heatmap_metrics_raw, ",")[[1]]
    heatmap_metrics <- trimws(heatmap_metrics[nzchar(heatmap_metrics)])
  } else {
    heatmap_metrics <- sort(unique(comparisons$metric))
  }

  for (metric_name in heatmap_metrics) {
    heatmap_path <- file.path(
      output_dir,
      paste0("phylum_comparisons_", heatmap_model, "_", metric_name, "_heatmap.pdf")
    )
    plot_phylum_heatmap(comparisons, metric_name, heatmap_model, heatmap_path)

    matrix_paths <- write_comparison_matrices(comparisons, metric_name, heatmap_model, output_dir)
    dendro_path <- file.path(
      output_dir,
      paste0("phylum_comparisons_", heatmap_model, "_", metric_name, "_heatmap_dendrogram.pdf")
    )
    plot_phylum_dendrogram(matrix_paths$diff_path, matrix_paths$pval_path, metric_name, heatmap_model, dendro_path)
  }
}

message("Analysis completed successfully.")
