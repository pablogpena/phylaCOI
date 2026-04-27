# -*- coding: utf-8 -*-
# Utilities for phase 2 haplotype clustering and current-zone analyses.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(igraph)
  library(ggplot2)
  library(geosphere)
  library(spdep)
  library(stringr)
  library(stringi)
})

# Check optional packages before running optional outputs.
# Returns TRUE when every package is available.
ensure_optional_packages <- function(packages, feature_name) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0) {
    warning(
      "Skipping ", feature_name, ". Missing R packages: ",
      paste(missing_packages, collapse = ", ")
    )
    return(FALSE)
  }

  TRUE
}

# Normalize locality names for robust metadata joins.
# Keeps comparisons stable across spaces, accents, and case.
normalize_locality_key <- function(values) {
  values %>%
    stringr::str_replace_all("\\s+", " ") %>%
    stringr::str_trim() %>%
    stringr::str_replace_all("\u00A0", " ") %>%
    stringi::stri_trans_general("NFD; [:Nonspacing Mark:] Remove; NFC") %>%
    toupper()
}

# Create a filesystem-safe output label.
# Used as suffix in phase 2 output filenames.
normalize_output_label <- function(label) {
  cleaned_label <- tolower(gsub("[^A-Za-z0-9]+", "_", label))
  cleaned_label <- gsub("^_+|_+$", "", cleaned_label)

  if (!nzchar(cleaned_label)) {
    return("haplotypes")
  }

  cleaned_label
}

# Load per-phylum haplotype network files from data/otus.
# Returns one edge table and one point table with a Filo column.
load_haplotype_network_data <- function(otus_root) {
  phylum_dirs <- list.dirs(otus_root, recursive = FALSE, full.names = TRUE)

  all_edges <- list()
  all_points <- list()

  for (phylum_dir in phylum_dirs) {
    phylum_name <- basename(phylum_dir)
    network_dir <- file.path(phylum_dir, "haplotype_network")
    edges_file <- file.path(network_dir, "edges_Mi_filo.csv")
    points_file <- file.path(network_dir, "points_Mi_filo.csv")

    if (!file.exists(edges_file) || !file.exists(points_file)) {
      message("Skipping ", phylum_name, ": missing haplotype network files.")
      next
    }

    message("Loading haplotype network data for: ", phylum_name)
    edges_df <- read.csv(edges_file, stringsAsFactors = FALSE)
    points_df <- read.csv(points_file, stringsAsFactors = FALSE)

    edges_df$Filo <- phylum_name
    points_df$Filo <- phylum_name

    all_edges[[phylum_name]] <- edges_df
    all_points[[phylum_name]] <- points_df
  }

  if (length(all_edges) == 0 || length(all_points) == 0) {
    stop("No haplotype network files found under: ", otus_root)
  }

  list(
    edges = dplyr::bind_rows(all_edges),
    points = dplyr::bind_rows(all_points)
  )
}

# Join edge endpoints with locality and OTU metadata.
# Drops intra-locality edges and cross-OTU edges.
build_edges_with_info <- function(edges_df_all, points_df_all) {
  required_edges <- c(
    "from", "to", "OTU_ID", "group", "distancia_genetica",
    "x", "y", "xend", "yend", "Filo"
  )
  required_points <- c("UniqueID", "Localities", "OTU_ID", "Filo")

  missing_edges <- setdiff(required_edges, names(edges_df_all))
  missing_points <- setdiff(required_points, names(points_df_all))

  if (length(missing_edges) > 0) {
    stop("Missing edge columns: ", paste(missing_edges, collapse = ", "))
  }
  if (length(missing_points) > 0) {
    stop("Missing point columns: ", paste(missing_points, collapse = ", "))
  }

  points_from <- points_df_all %>%
    dplyr::transmute(
      Filo = Filo,
      from = UniqueID,
      Locality_from = Localities,
      OTU_ID_from = as.character(OTU_ID)
    )

  points_to <- points_df_all %>%
    dplyr::transmute(
      Filo = Filo,
      to = UniqueID,
      Locality_to = Localities,
      OTU_ID_to = as.character(OTU_ID)
    )

  edges_df_all %>%
    dplyr::mutate(OTU_ID = as.character(OTU_ID)) %>%
    dplyr::left_join(points_from, by = c("Filo", "from")) %>%
    dplyr::left_join(points_to, by = c("Filo", "to")) %>%
    dplyr::filter(
      !is.na(Locality_from),
      !is.na(Locality_to),
      Locality_from != Locality_to,
      OTU_ID_from == OTU_ID_to
    ) %>%
    dplyr::mutate(
      OTU_ID = OTU_ID_from,
      Filo_from = Filo,
      Filo_to = Filo
    )
}

# Extract one coordinate per locality from both edge endpoints.
# Coordinates are averaged when a locality appears multiple times.
get_locality_coordinates <- function(edges_with_info) {
  dplyr::bind_rows(
    edges_with_info %>%
      dplyr::transmute(Locality = Locality_from, lon = x, lat = y),
    edges_with_info %>%
      dplyr::transmute(Locality = Locality_to, lon = xend, lat = yend)
  ) %>%
    dplyr::group_by(Locality) %>%
    dplyr::summarise(
      lon = mean(lon, na.rm = TRUE),
      lat = mean(lat, na.rm = TRUE),
      .groups = "drop"
    )
}

# Collapse sequence-level edges into locality-level links.
# Keeps mean genetic distance and the number of connections.
collapse_edges_to_localities <- function(edges_with_info) {
  edges_with_info %>%
    dplyr::filter(Locality_from != Locality_to) %>%
    dplyr::group_by(Locality_from, Locality_to) %>%
    dplyr::summarise(
      mean_dist = mean(distancia_genetica, na.rm = TRUE),
      n_connections = dplyr::n(),
      .groups = "drop"
    )
}

# Build an igraph locality network with the requested edge weights.
# Supported weights are inverse distance, counts, linear similarity, and RBF.
make_graph_with_weights <- function(
  edges_locality,
  mode = c("inverse", "count", "linear", "rbf"),
  sigma = NULL,
  eps = 1e-6
) {
  mode <- match.arg(mode)

  if (nrow(edges_locality) == 0) {
    return(igraph::make_empty_graph(directed = FALSE))
  }

  graph_obj <- igraph::graph_from_data_frame(
    edges_locality %>% dplyr::select(Locality_from, Locality_to),
    directed = FALSE
  )

  edge_distances <- edges_locality$mean_dist
  edge_distances[!is.finite(edge_distances)] <- 0

  if (mode == "inverse") {
    edge_weights <- 1 / (edge_distances + eps)
  } else if (mode == "count") {
    edge_weights <- edges_locality$n_connections
  } else if (mode == "linear") {
    max_distance <- max(edge_distances, na.rm = TRUE)
    if (!is.finite(max_distance) || max_distance == 0) {
      max_distance <- 1
    }
    edge_weights <- 1 - edge_distances / max_distance
  } else {
    if (is.null(sigma) || !is.finite(sigma) || sigma <= 0) {
      stop("A positive sigma is required for RBF weights.")
    }
    edge_weights <- exp(-edge_distances / sigma)
  }

  igraph::E(graph_obj)$weight <- edge_weights
  igraph::E(graph_obj)$dist <- edge_distances
  graph_obj
}

# Cluster a graph and attach locality coordinates.
# Returns the Louvain object and a Locality-lon-lat-cluster table.
cluster_and_coordinates <- function(
  graph_obj,
  edges_with_info,
  cluster_col = "Cluster",
  seed = 42
) {
  empty_coords <- tibble::tibble(
    Locality = character(),
    lon = numeric(),
    lat = numeric(),
    !!cluster_col := integer()
  )

  if (igraph::vcount(graph_obj) == 0 || igraph::ecount(graph_obj) == 0) {
    return(list(cluster = NULL, coords = empty_coords, modularity = NA_real_))
  }

  set.seed(seed)
  cluster_obj <- igraph::cluster_louvain(
    graph_obj,
    weights = igraph::E(graph_obj)$weight
  )

  membership_df <- tibble::tibble(
    Locality = names(igraph::membership(cluster_obj)),
    Cluster = as.integer(igraph::membership(cluster_obj))
  )

  coords <- get_locality_coordinates(edges_with_info)

  locality_coords <- membership_df %>%
    dplyr::left_join(coords, by = "Locality") %>%
    dplyr::relocate(lon, lat, .after = Locality)

  names(locality_coords)[names(locality_coords) == "Cluster"] <- cluster_col

  list(
    cluster = cluster_obj,
    coords = locality_coords,
    modularity = igraph::modularity(
      graph_obj,
      igraph::membership(cluster_obj),
      weights = igraph::E(graph_obj)$weight
    )
  )
}

# Compute adjusted Rand index without external dependencies.
# Returns NA when the two partitions cannot be compared.
adjusted_rand_index <- function(labels_a, labels_b) {
  valid_idx <- !is.na(labels_a) & !is.na(labels_b)
  labels_a <- labels_a[valid_idx]
  labels_b <- labels_b[valid_idx]

  if (length(labels_a) < 2 || length(labels_a) != length(labels_b)) {
    return(NA_real_)
  }

  contingency <- table(labels_a, labels_b)
  choose_two <- function(values) values * (values - 1) / 2

  total_pairs <- choose_two(sum(contingency))
  if (total_pairs == 0) {
    return(NA_real_)
  }

  sum_cells <- sum(choose_two(contingency))
  sum_rows <- sum(choose_two(rowSums(contingency)))
  sum_cols <- sum(choose_two(colSums(contingency)))
  expected_index <- sum_rows * sum_cols / total_pairs
  max_index <- 0.5 * (sum_rows + sum_cols)
  denominator <- max_index - expected_index

  if (denominator == 0) {
    return(ifelse(identical(labels_a, labels_b), 1, 0))
  }

  (sum_cells - expected_index) / denominator
}

# Compare two locality cluster tables with ARI.
# Localities are aligned before computing the score.
pairwise_ari <- function(coords_a, coords_b, col_a = "Cluster", col_b = "Cluster") {
  if (nrow(coords_a) == 0 || nrow(coords_b) == 0) {
    return(NA_real_)
  }

  comparison_df <- coords_a %>%
    dplyr::select(Locality, ClusterA = tidyselect::all_of(col_a)) %>%
    dplyr::inner_join(
      coords_b %>% dplyr::select(Locality, ClusterB = tidyselect::all_of(col_b)),
      by = "Locality"
    )

  adjusted_rand_index(comparison_df$ClusterA, comparison_df$ClusterB)
}

# Format small p-values consistently for output tables.
# Used by Moran's I summary tables.
format_small_p_values <- function(
  values,
  decimals_fixed = 8,
  tiny_thresh = 1e-12
) {
  ifelse(
    values < tiny_thresh,
    paste0("<", formatC(tiny_thresh, format = "e", digits = 0)),
    sprintf(paste0("%.", decimals_fixed, "f"), values)
  )
}

# Compute Moran's I for every cluster indicator.
# Returns a formatted table ready to write as CSV.
moran_table_by_cluster <- function(
  locality_coords_df,
  cluster_col = "Cluster",
  k = 5,
  jitter_amount = 1e-6,
  round_I = 3
) {
  empty_moran <- tibble::tibble(
    Cluster = character(),
    I = numeric(),
    p = character(),
    Freq = integer(),
    q_FDR = character()
  )

  required_cols <- c("Locality", "lon", "lat", cluster_col)
  if (!all(required_cols %in% names(locality_coords_df))) {
    return(empty_moran)
  }

  coords_all <- locality_coords_df %>%
    dplyr::select(Locality, lon, lat, Cluster = tidyselect::all_of(cluster_col)) %>%
    tidyr::drop_na(lon, lat, Cluster)

  if (nrow(coords_all) < 4 || length(unique(coords_all$Cluster)) < 2) {
    return(empty_moran)
  }

  coordinates_matrix <- as.matrix(coords_all[, c("lon", "lat")])
  if (jitter_amount > 0) {
    set.seed(1)
    coordinates_matrix <- jitter(coordinates_matrix, amount = jitter_amount)
  }

  neighbor_count <- max(1, min(k, nrow(coords_all) - 1))
  neighbors <- spdep::knn2nb(spdep::knearneigh(coordinates_matrix, k = neighbor_count))
  weights <- spdep::nb2listw(neighbors, style = "W", zero.policy = TRUE)

  cluster_labels <- factor(coords_all$Cluster)
  indicator_matrix <- model.matrix(~ cluster_labels - 1)
  colnames(indicator_matrix) <- levels(cluster_labels)

  moran_matrix <- t(vapply(
    seq_len(ncol(indicator_matrix)),
    function(cluster_idx) {
      test_result <- tryCatch(
        spdep::moran.test(
          as.numeric(indicator_matrix[, cluster_idx]),
          weights,
          zero.policy = TRUE
        ),
        error = function(error_condition) NULL
      )

      if (is.null(test_result)) {
        return(c(I = NA_real_, p = NA_real_))
      }

      c(
        I = unname(test_result$estimate["Moran I statistic"]),
        p = test_result$p.value
      )
    },
    numeric(2)
  ))

  moran_df <- as.data.frame(moran_matrix)
  moran_df$Cluster <- colnames(indicator_matrix)

  freq_df <- tibble::tibble(Cluster = as.character(cluster_labels)) %>%
    dplyr::count(Cluster, name = "Freq")

  moran_df %>%
    dplyr::select(Cluster, I, p) %>%
    dplyr::left_join(freq_df, by = "Cluster") %>%
    dplyr::mutate(q_FDR = p.adjust(p, method = "BH")) %>%
    dplyr::arrange(dplyr::desc(I)) %>%
    dplyr::mutate(
      I = round(I, round_I),
      p = format_small_p_values(p),
      q_FDR = format_small_p_values(q_FDR)
    )
}

# Run all-haplotype clustering with the four legacy weighting schemes.
# Returns cluster tables, modularity summaries, and selected RBF sigma.
run_all_haplotype_partitions <- function(
  edges_with_info,
  sigma_grid = c(0.005, 0.010, 0.015),
  eps = 1e-6,
  seed = 42
) {
  edges_locality <- collapse_edges_to_localities(edges_with_info)

  inverse_graph <- make_graph_with_weights(edges_locality, "inverse", eps = eps)
  count_graph <- make_graph_with_weights(edges_locality, "count", eps = eps)
  linear_graph <- make_graph_with_weights(edges_locality, "linear", eps = eps)

  inverse_result <- cluster_and_coordinates(
    inverse_graph,
    edges_with_info,
    "Cluster_inv",
    seed
  )
  count_result <- cluster_and_coordinates(
    count_graph,
    edges_with_info,
    "Cluster_cnt",
    seed
  )
  linear_result <- cluster_and_coordinates(
    linear_graph,
    edges_with_info,
    "Cluster_linear",
    seed
  )

  rbf_grid <- dplyr::bind_rows(lapply(sigma_grid, function(sigma_value) {
    rbf_graph_tmp <- make_graph_with_weights(
      edges_locality,
      "rbf",
      sigma = sigma_value,
      eps = eps
    )
    rbf_result_tmp <- cluster_and_coordinates(
      rbf_graph_tmp,
      edges_with_info,
      "Cluster_rbf",
      seed
    )
    tibble::tibble(
      sigma = sigma_value,
      modularidad = rbf_result_tmp$modularity,
      n_comun = ifelse(
        is.null(rbf_result_tmp$cluster),
        NA_integer_,
        length(unique(igraph::membership(rbf_result_tmp$cluster)))
      )
    )
  }))

  best_sigma <- rbf_grid$sigma[which.max(rbf_grid$modularidad)]
  if (length(best_sigma) == 0 || !is.finite(best_sigma)) {
    best_sigma <- sigma_grid[1]
  }

  rbf_graph <- make_graph_with_weights(
    edges_locality,
    "rbf",
    sigma = best_sigma,
    eps = eps
  )
  rbf_result <- cluster_and_coordinates(
    rbf_graph,
    edges_with_info,
    "Cluster_rbf",
    seed
  )

  comp_all <- inverse_result$coords %>%
    dplyr::select(Locality, Cluster_inv) %>%
    dplyr::inner_join(
      count_result$coords %>% dplyr::select(Locality, Cluster_cnt),
      by = "Locality"
    ) %>%
    dplyr::inner_join(
      linear_result$coords %>% dplyr::select(Locality, Cluster_linear),
      by = "Locality"
    ) %>%
    dplyr::inner_join(
      rbf_result$coords %>% dplyr::select(Locality, Cluster_rbf),
      by = "Locality"
    )

  ari_table <- tibble::tibble(
    Comparacion = c(
      "Inv vs Conex", "Inv vs Lineal", "Inv vs RBF",
      "Conex vs Lineal", "Conex vs RBF", "Lineal vs RBF"
    ),
    ARI = c(
      adjusted_rand_index(comp_all$Cluster_inv, comp_all$Cluster_cnt),
      adjusted_rand_index(comp_all$Cluster_inv, comp_all$Cluster_linear),
      adjusted_rand_index(comp_all$Cluster_inv, comp_all$Cluster_rbf),
      adjusted_rand_index(comp_all$Cluster_cnt, comp_all$Cluster_linear),
      adjusted_rand_index(comp_all$Cluster_cnt, comp_all$Cluster_rbf),
      adjusted_rand_index(comp_all$Cluster_linear, comp_all$Cluster_rbf)
    )
  )

  mod_summary <- dplyr::bind_rows(
    tibble::tibble(
      metodo = "Inverse_1_d",
      modularidad = inverse_result$modularity,
      n_comun = ifelse(
        is.null(inverse_result$cluster),
        NA_integer_,
        length(unique(igraph::membership(inverse_result$cluster)))
      )
    ),
    tibble::tibble(
      metodo = "Conexiones",
      modularidad = count_result$modularity,
      n_comun = ifelse(
        is.null(count_result$cluster),
        NA_integer_,
        length(unique(igraph::membership(count_result$cluster)))
      )
    ),
    tibble::tibble(
      metodo = "Lineal_1_minus_d_over_max",
      modularidad = linear_result$modularity,
      n_comun = ifelse(
        is.null(linear_result$cluster),
        NA_integer_,
        length(unique(igraph::membership(linear_result$cluster)))
      )
    ),
    tibble::tibble(
      metodo = paste0("RBF_sigma_", best_sigma),
      modularidad = rbf_result$modularity,
      n_comun = ifelse(
        is.null(rbf_result$cluster),
        NA_integer_,
        length(unique(igraph::membership(rbf_result$cluster)))
      )
    )
  )

  mod_results <- tibble::tibble(
    Metodo = c("Inverse_1_d", "Conexiones", "Lineal", "RBF"),
    Modularity = c(
      inverse_result$modularity,
      count_result$modularity,
      linear_result$modularity,
      rbf_result$modularity
    )
  )

  list(
    edges_locality = edges_locality,
    comp_all = comp_all,
    ari_table = ari_table,
    mod_summary = mod_summary,
    mod_results = mod_results,
    rbf_grid = rbf_grid,
    best_sigma = best_sigma,
    results = list(
      inverse = inverse_result,
      count = count_result,
      linear = linear_result,
      rbf = rbf_result
    )
  )
}

# Read optional ocean metadata and normalize locality keys.
# Expects id_sample and Oceans columns in a semicolon-separated CSV.
read_ocean_metadata <- function(metadata_file) {
  if (is.null(metadata_file)) {
    return(NULL)
  }

  ocean_df <- read.csv(metadata_file, sep = ";", stringsAsFactors = FALSE)
  required_cols <- c("id_sample", "Oceans")
  missing_cols <- setdiff(required_cols, names(ocean_df))

  if (length(missing_cols) > 0) {
    stop("Missing ocean metadata columns: ", paste(missing_cols, collapse = ", "))
  }

  ocean_df %>%
    dplyr::transmute(
      Locality = normalize_locality_key(id_sample),
      Oceans = stringr::str_trim(gsub("\\s+", " ", Oceans))
    ) %>%
    dplyr::distinct(Locality, .keep_all = TRUE)
}

# Compute ocean percentages per cluster.
# Returns NULL when metadata is unavailable.
ocean_share_by_cluster <- function(coords_df, cluster_col, ocean_df) {
  if (is.null(ocean_df)) {
    return(NULL)
  }

  required_cols <- c("Locality", cluster_col)
  if (!all(required_cols %in% names(coords_df))) {
    return(NULL)
  }

  coords_norm <- coords_df %>%
    dplyr::select(Locality, Cluster = tidyselect::all_of(cluster_col)) %>%
    dplyr::mutate(Locality = normalize_locality_key(Locality))

  missing_in_ocean <- coords_norm %>%
    dplyr::anti_join(ocean_df, by = "Locality") %>%
    dplyr::distinct(Locality)

  if (nrow(missing_in_ocean) > 0) {
    message(
      "[", cluster_col, "] Localities missing in ocean metadata: ",
      nrow(missing_in_ocean)
    )
  }

  coords_norm %>%
    dplyr::inner_join(ocean_df, by = "Locality") %>%
    dplyr::mutate(Cluster = as.integer(Cluster)) %>%
    dplyr::group_by(Cluster, Oceans) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
    dplyr::mutate(pct = 100 * n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::select(Cluster, Oceans, pct) %>%
    tidyr::pivot_wider(
      names_from = Oceans,
      values_from = pct,
      values_fill = 0
    ) %>%
    dplyr::mutate(dplyr::across(-Cluster, ~ round(.x, 1))) %>%
    dplyr::arrange(Cluster)
}

# Assign localities to broad current/front zones.
# Uses the same longitude/latitude thresholds as the legacy scripts.
build_current_zones <- function(locality_coords_df) {
  locality_coords_df %>%
    dplyr::select(Locality, lon, lat) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      Zone = dplyr::case_when(
        lon < 0 & lat >= 37 ~ "ATL_N_UPW",
        lon < -2 & lat < 37 ~ "ATL_S_UPW_ALM",
        lon >= -2 ~ "MED",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(Locality, Zone)
}

# Compute current-zone percentages per cluster.
# Ensures every expected zone column exists in the output.
zone_share_by_cluster <- function(coords_df, cluster_col, zones_df) {
  zone_levels <- c("ATL_N_UPW", "ATL_S_UPW_ALM", "MED")

  if (!all(c("Locality", cluster_col) %in% names(coords_df))) {
    return(tibble::tibble(Cluster = integer()))
  }

  zone_share <- coords_df %>%
    dplyr::select(Locality, Cluster = tidyselect::all_of(cluster_col)) %>%
    dplyr::inner_join(zones_df, by = "Locality") %>%
    dplyr::filter(!is.na(Zone)) %>%
    dplyr::count(Cluster, Zone) %>%
    dplyr::group_by(Cluster) %>%
    dplyr::mutate(pct = 100 * n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-n) %>%
    tidyr::pivot_wider(
      names_from = Zone,
      values_from = pct,
      values_fill = 0
    )

  for (zone_name in zone_levels) {
    if (!zone_name %in% names(zone_share)) {
      zone_share[[zone_name]] <- 0
    }
  }

  zone_share %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 1))) %>%
    dplyr::select(Cluster, tidyselect::all_of(zone_levels)) %>%
    dplyr::arrange(Cluster)
}

# Add coordinates to the combined cluster table.
# Used by map plotting and compact downstream exports.
attach_coordinates_to_clusters <- function(comp_all, edges_with_info) {
  coords_ref <- get_locality_coordinates(edges_with_info)

  comp_all %>%
    dplyr::left_join(coords_ref, by = "Locality")
}

# Build a stable color palette across all clustering methods.
# Starts with the legacy 16 colors and expands if needed.
make_cluster_palette <- function(comp_all, hcl_palette = "Dark 3") {
  base_palette <- c(
    "1" = "#8DD3C7", "2" = "#FDB462", "3" = "#B3DE69", "4" = "#FCCDE5",
    "5" = "#D9D9D9", "6" = "#BC80BD", "7" = "#CCEBC5", "8" = "#FFED6F",
    "9" = "#7F7F7F", "10" = "#FFFFB3", "11" = "#BEBADA", "12" = "#fb8072",
    "13" = "#80B1D3", "14" = "#C49C94", "15" = "#AEC7E8", "16" = "#304847"
  )

  cluster_cols <- c("Cluster_inv", "Cluster_cnt", "Cluster_linear", "Cluster_rbf")
  cluster_levels <- sort(unique(unlist(lapply(
    comp_all[, cluster_cols, drop = FALSE],
    function(values) as.character(unique(values))
  ))))
  cluster_levels <- cluster_levels[!is.na(cluster_levels)]

  missing_levels <- setdiff(cluster_levels, names(base_palette))
  if (length(missing_levels) > 0) {
    extra_colors <- grDevices::hcl.colors(length(missing_levels), palette = hcl_palette)
    names(extra_colors) <- missing_levels
    base_palette <- c(base_palette, extra_colors)
  }

  base_palette[cluster_levels]
}

# Prepare Iberian map layers for static clustering maps.
# Uses Natural Earth polygons and clips to the legacy bounding box.
prepare_base_layers <- function(
  xmin = -10,
  ymin = 35.5,
  xmax = 4,
  ymax = 44.5,
  crs_planar = 3035
) {
  if (!ensure_optional_packages(c("sf", "rnaturalearth"), "static maps")) {
    return(NULL)
  }

  old_s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)

  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  bbox_wgs <- sf::st_as_sfc(
    sf::st_bbox(c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), crs = 4326)
  )

  iberia_wgs <- world %>%
    dplyr::filter(admin %in% c("Spain", "Portugal")) %>%
    sf::st_union() %>%
    sf::st_as_sf()
  france_wgs <- world %>%
    dplyr::filter(admin == "France") %>%
    sf::st_union() %>%
    sf::st_as_sf()
  andorra_wgs <- world %>%
    dplyr::filter(admin == "Andorra" | sovereignt == "Andorra") %>%
    sf::st_make_valid()

  iberia_wgs_clip <- suppressWarnings(sf::st_intersection(iberia_wgs, bbox_wgs))
  france_wgs_clip <- suppressWarnings(sf::st_intersection(france_wgs, bbox_wgs))
  andorra_wgs_clip <- suppressWarnings(sf::st_intersection(andorra_wgs, bbox_wgs))
  land_wgs_clip <- suppressWarnings(
    sf::st_union(iberia_wgs_clip, france_wgs_clip, andorra_wgs_clip)
  ) %>%
    sf::st_as_sf()

  list(
    iberia_wgs_clip = iberia_wgs_clip,
    france_wgs_clip = france_wgs_clip,
    andorra_wgs_clip = andorra_wgs_clip,
    land_wgs_clip = land_wgs_clip,
    bbox_wgs = bbox_wgs,
    crs_planar = crs_planar
  )
}

# Build one static cluster map.
# Combines marine Voronoi shading and hex-cell cluster pies.
make_map_for_case <- function(
  comp_all,
  cluster_col,
  title_text,
  base_layers = NULL,
  cellsize = 60000,
  pal = NULL,
  voronoi_alpha = 0.25
) {
  if (!ensure_optional_packages(c("sf", "scatterpie"), "static maps")) {
    return(NULL)
  }
  if (!all(c("Locality", "lon", "lat", cluster_col) %in% names(comp_all))) {
    return(NULL)
  }

  df_case <- comp_all %>%
    dplyr::select(Locality, lon, lat, Cluster = tidyselect::all_of(cluster_col)) %>%
    tidyr::drop_na(lon, lat, Cluster)

  if (nrow(df_case) < 2) {
    warning("Skipping map for ", cluster_col, ": fewer than two localities.")
    return(NULL)
  }

  if (is.null(pal)) {
    cluster_levels <- sort(unique(as.character(df_case$Cluster)))
    pal <- grDevices::hcl.colors(length(cluster_levels), palette = "Dark 3")
    names(pal) <- cluster_levels
  }

  df_case$Cluster <- factor(as.character(df_case$Cluster), levels = names(pal))

  if (is.null(base_layers)) {
    base_layers <- prepare_base_layers()
  }
  if (is.null(base_layers)) {
    return(NULL)
  }

  points_wgs <- sf::st_as_sf(df_case, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

  iberia <- sf::st_transform(base_layers$iberia_wgs_clip, base_layers$crs_planar)
  france <- sf::st_transform(base_layers$france_wgs_clip, base_layers$crs_planar)
  andorra <- sf::st_transform(base_layers$andorra_wgs_clip, base_layers$crs_planar)
  land <- sf::st_transform(base_layers$land_wgs_clip, base_layers$crs_planar)
  bbox_pr <- sf::st_transform(base_layers$bbox_wgs, base_layers$crs_planar)
  points_projected <- sf::st_transform(points_wgs, base_layers$crs_planar)

  andorra <- suppressWarnings(sf::st_make_valid(andorra))
  need_andorra_fallback <- nrow(andorra) == 0 ||
    any(sf::st_is_empty(andorra)) ||
    as.numeric(sf::st_area(sf::st_union(andorra))) < 5e7

  if (need_andorra_fallback) {
    andorra_point <- sf::st_as_sf(
      data.frame(id = 1, lon = 1.6, lat = 42.55),
      coords = c("lon", "lat"),
      crs = 4326
    )
    andorra <- andorra_point %>%
      sf::st_transform(base_layers$crs_planar) %>%
      sf::st_buffer(dist = 12000) %>%
      sf::st_as_sf()
  }

  voronoi_polygons <- sf::st_voronoi(sf::st_union(points_projected)) %>%
    sf::st_collection_extract("POLYGON") %>%
    sf::st_as_sf()
  sf::st_crs(voronoi_polygons) <- sf::st_crs(points_projected)

  nearest_points <- sf::st_nearest_feature(voronoi_polygons, points_projected)
  voronoi_polygons <- voronoi_polygons %>%
    dplyr::mutate(
      Cluster = as.character(points_projected$Cluster[nearest_points]),
      Cluster = factor(Cluster, levels = names(pal))
    )

  sea_mask <- sf::st_difference(bbox_pr, sf::st_buffer(land, 0))
  voronoi_sea <- suppressWarnings(sf::st_intersection(voronoi_polygons, sea_mask))
  if (inherits(voronoi_sea, "sfc")) {
    voronoi_sea <- sf::st_as_sf(voronoi_sea)
  }
  if (!inherits(voronoi_sea, "sf") || nrow(voronoi_sea) == 0) {
    voronoi_sea <- NULL
  }

  hex_grid <- sf::st_make_grid(bbox_pr, cellsize = cellsize, what = "polygons", square = FALSE)
  hex_grid <- sf::st_as_sf(hex_grid)
  sf::st_crs(hex_grid) <- sf::st_crs(bbox_pr)
  hex_grid <- sf::st_intersection(hex_grid, bbox_pr) %>%
    dplyr::mutate(id = dplyr::row_number())
  hex_centroids <- sf::st_centroid(hex_grid)

  joined_points <- sf::st_join(points_projected, hex_grid, join = sf::st_within)
  cluster_counts <- joined_points %>%
    sf::st_drop_geometry() %>%
    dplyr::filter(!is.na(id)) %>%
    dplyr::count(id, Cluster = as.character(Cluster), name = "n")

  totals <- cluster_counts %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(total = sum(n), .groups = "drop")

  cluster_wide <- cluster_counts %>%
    dplyr::group_by(id, Cluster) %>%
    dplyr::summarise(n = sum(n), .groups = "drop_last") %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(
      id_cols = id,
      names_from = Cluster,
      values_from = prop,
      values_fill = 0
    )

  pie_cols <- names(pal)[names(pal) %in% names(cluster_wide)]

  hex_data <- hex_centroids %>%
    dplyr::left_join(cluster_wide, by = "id") %>%
    dplyr::left_join(totals, by = "id")

  hex_coords <- sf::st_coordinates(sf::st_geometry(hex_data))
  hex_df <- hex_data %>%
    sf::st_drop_geometry() %>%
    dplyr::mutate(x = hex_coords[, 1], y = hex_coords[, 2]) %>%
    dplyr::mutate(dplyr::across(tidyselect::all_of(pie_cols), ~ tidyr::replace_na(as.numeric(.), 0))) %>%
    dplyr::mutate(
      total = tidyr::replace_na(total, 0),
      radius = 15000 + sqrt(total) * 1500
    ) %>%
    dplyr::filter(total > 0)

  plot_limits <- sf::st_bbox(bbox_pr)

  plot_obj <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = france, fill = "white", color = "grey75", size = 0.25)

  if (!is.null(voronoi_sea)) {
    plot_obj <- plot_obj +
      ggplot2::geom_sf(
        data = voronoi_sea,
        ggplot2::aes(fill = Cluster),
        alpha = voronoi_alpha,
        color = NA
      )
  }

  plot_obj +
    ggplot2::geom_sf(data = iberia, fill = "white", color = "grey30", size = 0.4) +
    ggplot2::geom_sf(data = andorra, fill = "white", color = "grey30", size = 0.35) +
    scatterpie::geom_scatterpie(
      data = hex_df,
      ggplot2::aes(x = x, y = y, r = radius),
      cols = pie_cols,
      color = NA,
      alpha = 1
    ) +
    ggplot2::scale_fill_manual(values = pal, name = "Cluster", drop = FALSE) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 1))) +
    ggplot2::coord_sf(
      xlim = c(plot_limits["xmin"], plot_limits["xmax"]),
      ylim = c(plot_limits["ymin"], plot_limits["ymax"]),
      expand = FALSE
    ) +
    ggplot2::labs(
      title = title_text,
      subtitle = "Voronoi shading + cluster proportions by hex cell",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey85", size = 0.2)
    )
}

# Write the four all-haplotype static maps.
# Skips cleanly when optional map packages are unavailable.
write_all_haplotype_maps <- function(comp_all, output_dir, label) {
  if (!ensure_optional_packages(c("sf", "scatterpie", "rnaturalearth"), "static maps")) {
    return(invisible(NULL))
  }

  figures_dir <- file.path(output_dir, "figures")
  dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

  palette_all <- make_cluster_palette(comp_all)
  base_layers <- prepare_base_layers()
  if (is.null(base_layers)) {
    return(invisible(NULL))
  }

  map_specs <- list(
    inverse = list(
      col = "Cluster_inv",
      title = "Clusters by genetic distance (Inverse 1/d)",
      file = paste0("map_inverse_1d_", label, ".pdf")
    ),
    count = list(
      col = "Cluster_cnt",
      title = "Clusters by number of connections",
      file = paste0("map_connections_", label, ".pdf")
    ),
    linear = list(
      col = "Cluster_linear",
      title = "Clusters by linear similarity (1 - d/max)",
      file = paste0("map_linear_", label, ".pdf")
    ),
    rbf = list(
      col = "Cluster_rbf",
      title = "Clusters by RBF similarity (exp(-d/sigma))",
      file = paste0("map_rbf_", label, ".pdf")
    )
  )

  for (map_spec in map_specs) {
    map_plot <- tryCatch(
      make_map_for_case(
        comp_all,
        map_spec$col,
        map_spec$title,
        base_layers = base_layers,
        pal = palette_all,
        voronoi_alpha = 0.25
      ),
      error = function(error_condition) {
        warning("Could not build map ", map_spec$file, ": ", conditionMessage(error_condition))
        NULL
      }
    )

    if (!is.null(map_plot)) {
      ggplot2::ggsave(
        filename = file.path(figures_dir, map_spec$file),
        plot = map_plot,
        width = 8,
        height = 7
      )
    }
  }

  invisible(NULL)
}

# Run one same/diff case for a selected edge subset.
# Returns clustering, modularity, selected sigma, and Moran's I.
analyze_haplotype_case <- function(
  edges_sub,
  mode = c("inv", "count", "linear", "rbf"),
  sigma_grid = c(0.005, 0.010, 0.015),
  eps = 1e-6,
  k_moran = 5,
  seed = 42
) {
  mode <- match.arg(mode)
  edges_locality <- collapse_edges_to_localities(edges_sub)

  if (nrow(edges_locality) == 0) {
    empty_coords <- tibble::tibble(
      Locality = character(),
      lon = numeric(),
      lat = numeric(),
      Cluster = integer()
    )
    return(list(
      mode = mode,
      sigma = NA_real_,
      graph = igraph::make_empty_graph(directed = FALSE),
      cluster = NULL,
      modularity = NA_real_,
      locality_coords = empty_coords,
      moran = moran_table_by_cluster(empty_coords, "Cluster", k = k_moran)
    ))
  }

  if (mode == "rbf") {
    rbf_grid <- dplyr::bind_rows(lapply(sigma_grid, function(sigma_value) {
      graph_tmp <- make_graph_with_weights(
        edges_locality,
        "rbf",
        sigma = sigma_value,
        eps = eps
      )
      result_tmp <- cluster_and_coordinates(graph_tmp, edges_sub, "Cluster", seed)
      tibble::tibble(
        sigma = sigma_value,
        modularidad = result_tmp$modularity,
        n_comun = ifelse(
          is.null(result_tmp$cluster),
          NA_integer_,
          length(unique(igraph::membership(result_tmp$cluster)))
        )
      )
    }))
    best_sigma <- rbf_grid$sigma[which.max(rbf_grid$modularidad)]
    if (length(best_sigma) == 0 || !is.finite(best_sigma)) {
      best_sigma <- sigma_grid[1]
    }
    graph_obj <- make_graph_with_weights(
      edges_locality,
      "rbf",
      sigma = best_sigma,
      eps = eps
    )
  } else {
    graph_mode <- dplyr::case_when(
      mode == "inv" ~ "inverse",
      mode == "count" ~ "count",
      mode == "linear" ~ "linear",
      TRUE ~ "inverse"
    )
    best_sigma <- NA_real_
    graph_obj <- make_graph_with_weights(edges_locality, graph_mode, eps = eps)
  }

  result <- cluster_and_coordinates(graph_obj, edges_sub, "Cluster", seed)
  result$moran <- moran_table_by_cluster(result$coords, "Cluster", k = k_moran)

  list(
    mode = mode,
    sigma = best_sigma,
    graph = graph_obj,
    cluster = result$cluster,
    modularity = result$modularity,
    locality_coords = result$coords,
    moran = result$moran
  )
}

# Count communities in a clustering result.
# Returns NA when no clustering was produced.
count_communities <- function(result) {
  if (is.null(result$cluster)) {
    return(NA_integer_)
  }

  length(unique(igraph::membership(result$cluster)))
}

# Run all same-vs-diff clustering comparisons.
# Produces results for SAME and DIFF across all four methods.
run_same_diff_partitions <- function(
  edges_with_info,
  sigma_grid = c(0.005, 0.010, 0.015),
  eps = 1e-6,
  k_moran = 5,
  seed = 42
) {
  edges_same <- edges_with_info %>% dplyr::filter(group == "same")
  edges_diff <- edges_with_info %>% dplyr::filter(group == "diff")

  same_results <- list(
    inverse = analyze_haplotype_case(edges_same, "inv", sigma_grid, eps, k_moran, seed),
    count = analyze_haplotype_case(edges_same, "count", sigma_grid, eps, k_moran, seed),
    linear = analyze_haplotype_case(edges_same, "linear", sigma_grid, eps, k_moran, seed),
    rbf = analyze_haplotype_case(edges_same, "rbf", sigma_grid, eps, k_moran, seed)
  )

  diff_results <- list(
    inverse = analyze_haplotype_case(edges_diff, "inv", sigma_grid, eps, k_moran, seed),
    count = analyze_haplotype_case(edges_diff, "count", sigma_grid, eps, k_moran, seed),
    linear = analyze_haplotype_case(edges_diff, "linear", sigma_grid, eps, k_moran, seed),
    rbf = analyze_haplotype_case(edges_diff, "rbf", sigma_grid, eps, k_moran, seed)
  )

  mod_summary <- dplyr::bind_rows(
    tibble::tibble(
      grupo = "SAME",
      metodo = "Inverse_1/d",
      modularidad = same_results$inverse$modularity,
      n_comun = count_communities(same_results$inverse)
    ),
    tibble::tibble(
      grupo = "SAME",
      metodo = "Count",
      modularidad = same_results$count$modularity,
      n_comun = count_communities(same_results$count)
    ),
    tibble::tibble(
      grupo = "SAME",
      metodo = "Linear",
      modularidad = same_results$linear$modularity,
      n_comun = count_communities(same_results$linear)
    ),
    tibble::tibble(
      grupo = "SAME",
      metodo = paste0("RBF (sigma=", sprintf("%.3f", same_results$rbf$sigma), ")"),
      modularidad = same_results$rbf$modularity,
      n_comun = count_communities(same_results$rbf)
    ),
    tibble::tibble(
      grupo = "DIFF",
      metodo = "Inverse_1/d",
      modularidad = diff_results$inverse$modularity,
      n_comun = count_communities(diff_results$inverse)
    ),
    tibble::tibble(
      grupo = "DIFF",
      metodo = "Count",
      modularidad = diff_results$count$modularity,
      n_comun = count_communities(diff_results$count)
    ),
    tibble::tibble(
      grupo = "DIFF",
      metodo = "Linear",
      modularidad = diff_results$linear$modularity,
      n_comun = count_communities(diff_results$linear)
    ),
    tibble::tibble(
      grupo = "DIFF",
      metodo = paste0("RBF (sigma=", sprintf("%.3f", diff_results$rbf$sigma), ")"),
      modularidad = diff_results$rbf$modularity,
      n_comun = count_communities(diff_results$rbf)
    )
  )

  ari_same <- tibble::tibble(
    Comparacion = c(
      "Inv vs Count", "Inv vs Linear", "Inv vs RBF",
      "Count vs Linear", "Count vs RBF", "Linear vs RBF"
    ),
    ARI = c(
      pairwise_ari(same_results$inverse$locality_coords, same_results$count$locality_coords),
      pairwise_ari(same_results$inverse$locality_coords, same_results$linear$locality_coords),
      pairwise_ari(same_results$inverse$locality_coords, same_results$rbf$locality_coords),
      pairwise_ari(same_results$count$locality_coords, same_results$linear$locality_coords),
      pairwise_ari(same_results$count$locality_coords, same_results$rbf$locality_coords),
      pairwise_ari(same_results$linear$locality_coords, same_results$rbf$locality_coords)
    )
  )

  ari_diff <- tibble::tibble(
    Comparacion = c(
      "Inv vs Count", "Inv vs Linear", "Inv vs RBF",
      "Count vs Linear", "Count vs RBF", "Linear vs RBF"
    ),
    ARI = c(
      pairwise_ari(diff_results$inverse$locality_coords, diff_results$count$locality_coords),
      pairwise_ari(diff_results$inverse$locality_coords, diff_results$linear$locality_coords),
      pairwise_ari(diff_results$inverse$locality_coords, diff_results$rbf$locality_coords),
      pairwise_ari(diff_results$count$locality_coords, diff_results$linear$locality_coords),
      pairwise_ari(diff_results$count$locality_coords, diff_results$rbf$locality_coords),
      pairwise_ari(diff_results$linear$locality_coords, diff_results$rbf$locality_coords)
    )
  )

  ari_same_vs_diff <- tibble::tibble(
    Metodo = c("Inverse_1/d", "Count", "Linear", "RBF"),
    ARI = c(
      pairwise_ari(same_results$inverse$locality_coords, diff_results$inverse$locality_coords),
      pairwise_ari(same_results$count$locality_coords, diff_results$count$locality_coords),
      pairwise_ari(same_results$linear$locality_coords, diff_results$linear$locality_coords),
      pairwise_ari(same_results$rbf$locality_coords, diff_results$rbf$locality_coords)
    )
  )

  list(
    same = same_results,
    diff = diff_results,
    mod_summary = mod_summary,
    ari_same = ari_same,
    ari_diff = ari_diff,
    ari_same_vs_diff = ari_same_vs_diff
  )
}

# Compute geographical distance summaries for all/same/diff edges.
# Returns per-OTU medians, summary stats, and Welch-test inputs.
compute_geographical_distance_outputs <- function(edges_with_info) {
  empty_result <- list(
    med_geo_by_otu = tibble::tibble(),
    resumen_geo_otu = tibble::tibble(),
    dist_long = tibble::tibble(),
    sum_stats = tibble::tibble(),
    stats_summary = tibble::tibble(),
    t_test = NULL
  )

  if (nrow(edges_with_info) == 0) {
    return(empty_result)
  }

  edges_base <- edges_with_info %>%
    dplyr::mutate(
      geo_dist_m = geosphere::distHaversine(cbind(x, y), cbind(xend, yend)),
      geo_dist_km = geo_dist_m / 1000
    ) %>%
    dplyr::filter(Locality_from != Locality_to, geo_dist_km > 0)

  if (nrow(edges_base) == 0) {
    return(empty_result)
  }

  med_geo_by_otu <- edges_base %>%
    dplyr::group_by(OTU_ID) %>%
    dplyr::summarise(
      n_edges = dplyr::n(),
      n_diff = sum(group == "diff"),
      n_same = sum(group == "same"),
      med_km_all = median(geo_dist_km, na.rm = TRUE),
      med_km_diff = ifelse(
        n_diff > 0,
        median(geo_dist_km[group == "diff"], na.rm = TRUE),
        NA_real_
      ),
      med_km_same = ifelse(
        n_same > 0,
        median(geo_dist_km[group == "same"], na.rm = TRUE),
        NA_real_
      ),
      .groups = "drop"
    )

  resumen_geo_otu <- med_geo_by_otu %>%
    dplyr::summarise(
      n_OTUs = dplyr::n(),
      n_OTUs_con_diff = sum(!is.na(med_km_diff)),
      n_OTUs_con_same = sum(!is.na(med_km_same)),
      med_of_meds_all_km = median(med_km_all, na.rm = TRUE),
      med_of_meds_diff_km = median(med_km_diff, na.rm = TRUE),
      med_of_meds_same_km = median(med_km_same, na.rm = TRUE)
    )

  edges_geo <- edges_base %>% dplyr::select(group, geo_dist_km)
  dist_long <- dplyr::bind_rows(
    edges_geo %>% dplyr::mutate(tipo = "all"),
    edges_geo %>% dplyr::filter(group == "same") %>% dplyr::mutate(tipo = "same"),
    edges_geo %>% dplyr::filter(group == "diff") %>% dplyr::mutate(tipo = "diff")
  ) %>%
    dplyr::mutate(tipo = factor(tipo, levels = c("all", "same", "diff")))

  sum_stats <- dist_long %>%
    dplyr::group_by(tipo) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median_km = median(geo_dist_km, na.rm = TRUE),
      q1 = stats::quantile(geo_dist_km, 0.25, na.rm = TRUE),
      q3 = stats::quantile(geo_dist_km, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  dist_test <- dist_long %>%
    dplyr::filter(tipo %in% c("same", "diff")) %>%
    droplevels()

  stats_summary <- dist_test %>%
    dplyr::group_by(tipo) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median_km = median(geo_dist_km, na.rm = TRUE),
      IQR_km = IQR(geo_dist_km, na.rm = TRUE),
      mean_km = mean(geo_dist_km, na.rm = TRUE),
      .groups = "drop"
    )

  group_counts <- table(dist_test$tipo)
  t_test <- NULL
  if (length(group_counts) == 2 && all(group_counts > 1)) {
    t_test <- stats::t.test(geo_dist_km ~ tipo, data = dist_test, var.equal = FALSE)
  }

  list(
    med_geo_by_otu = med_geo_by_otu,
    resumen_geo_otu = resumen_geo_otu,
    dist_long = dist_long,
    sum_stats = sum_stats,
    stats_summary = stats_summary,
    t_test = t_test
  )
}

# Plot geographical distances for all/same/diff connections.
# Uses medians, IQR bars, and a Welch-test annotation when possible.
plot_geographical_distances <- function(geo_outputs, output_file) {
  if (nrow(geo_outputs$dist_long) == 0 || nrow(geo_outputs$sum_stats) == 0) {
    warning("Skipping geographical distance plot: no distance data.")
    return(invisible(NULL))
  }

  palette_values <- c(
    all = "#bdbdbd",
    same = "#756bb1",
    diff = "#de2d26"
  )

  p_value <- if (!is.null(geo_outputs$t_test)) geo_outputs$t_test$p.value else NA_real_
  significance_label <- dplyr::case_when(
    is.na(p_value) ~ "NA",
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "ns"
  )

  y_top <- max(geo_outputs$dist_long$geo_dist_km, na.rm = TRUE)
  y_label <- y_top * 1.05
  p_label <- ifelse(
    is.na(p_value),
    "Welch t-test: not available",
    paste0("Welch t-test: p = ", signif(p_value, 3), " (", significance_label, ")")
  )

  distance_plot <- ggplot2::ggplot() +
    ggplot2::geom_jitter(
      data = geo_outputs$dist_long,
      ggplot2::aes(x = tipo, y = geo_dist_km, color = tipo),
      width = 0.15,
      alpha = 0.10,
      size = 1,
      show.legend = FALSE
    ) +
    ggplot2::geom_col(
      data = geo_outputs$sum_stats,
      ggplot2::aes(x = tipo, y = median_km, fill = tipo),
      width = 0.6,
      alpha = 0.8,
      color = NA
    ) +
    ggplot2::geom_errorbar(
      data = geo_outputs$sum_stats,
      ggplot2::aes(x = tipo, ymin = q1, ymax = q3),
      width = 0.16,
      linewidth = 0.7
    ) +
    ggplot2::annotate("text", x = 2.5, y = y_label, label = p_label, size = 4.2) +
    ggplot2::scale_fill_manual(values = palette_values, guide = "none") +
    ggplot2::scale_color_manual(values = palette_values, guide = "none") +
    ggplot2::labs(
      title = "Geographical distances of connections",
      subtitle = "Bars = median; lines = IQR (Q1-Q3); points = all connections",
      x = "Connection set",
      y = "Geographical distance (km)"
    ) +
    ggplot2::expand_limits(y = y_label * 1.05) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())

  ggplot2::ggsave(
    filename = output_file,
    plot = distance_plot,
    width = 8,
    height = 6
  )

  invisible(distance_plot)
}

# Parse zone-share filenames from all/same/diff outputs.
# Supports both legacy all filenames and same/diff filenames.
parse_zone_share_filename <- function(file_path) {
  filename <- basename(file_path) %>%
    stringr::str_remove("^zone_share_") %>%
    stringr::str_remove("\\.csv$")
  filename_parts <- strsplit(filename, "_", fixed = TRUE)[[1]]

  if (length(filename_parts) < 2) {
    return(list(type = NA_character_, method = NA_character_, group = NA_character_))
  }

  if (filename_parts[1] %in% c("all", "same", "diff")) {
    type_value <- filename_parts[1]
    method_value <- filename_parts[2]
    group_value <- paste(filename_parts[-c(1, 2)], collapse = "_")
  } else {
    type_value <- "all"
    method_value <- filename_parts[1]
    group_value <- paste(filename_parts[-1], collapse = "_")
  }

  list(type = type_value, method = method_value, group = group_value)
}

# Plot current-zone percentages across methods and groups.
# Reads zone_share_*.csv files from the provided directories.
plot_current_zone_summary <- function(zone_dirs, output_file) {
  existing_dirs <- zone_dirs[dir.exists(zone_dirs)]
  zone_files <- unique(unlist(lapply(
    existing_dirs,
    list.files,
    pattern = "^zone_share_.*\\.csv$",
    full.names = TRUE
  )))

  if (length(zone_files) == 0) {
    warning("Skipping current-zone summary figure: no zone_share_*.csv files found.")
    return(invisible(NULL))
  }

  zone_levels <- c("ATL_N_UPW", "ATL_S_UPW_ALM", "MED")

  all_data <- dplyr::bind_rows(lapply(zone_files, function(file_path) {
    zone_df <- read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)
    file_info <- parse_zone_share_filename(file_path)

    for (zone_name in zone_levels) {
      if (!zone_name %in% names(zone_df)) {
        zone_df[[zone_name]] <- 0
      }
    }

    zone_df %>%
      dplyr::mutate(
        type = file_info$type,
        method = file_info$method,
        group = file_info$group
      )
  }))

  needed_cols <- c("Cluster", zone_levels, "type", "method", "group")
  missing_cols <- setdiff(needed_cols, names(all_data))
  if (length(missing_cols) > 0) {
    stop("Missing columns in zone-share files: ", paste(missing_cols, collapse = ", "))
  }

  long_data <- all_data %>%
    tidyr::pivot_longer(
      cols = tidyselect::all_of(zone_levels),
      names_to = "Zone",
      values_to = "pct"
    ) %>%
    dplyr::mutate(
      type = factor(type, levels = c("all", "same", "diff")),
      method = factor(method, levels = c("inverse", "count", "linear", "rbf")),
      group = factor(
        group,
        levels = c(
          "metazoa", "arqueaplastidia", "protists",
          "multicellular", "unicellular",
          setdiff(sort(unique(group)), c(
            "metazoa", "arqueaplastidia", "protists",
            "multicellular", "unicellular"
          ))
        )
      ),
      Zone = factor(Zone, levels = zone_levels),
      Cluster = factor(Cluster)
    )

  current_plot <- ggplot2::ggplot(
    long_data,
    ggplot2::aes(x = Cluster, y = pct, fill = Zone)
  ) +
    ggplot2::geom_col(width = 0.9) +
    ggplot2::facet_grid(
      group ~ method + type,
      scales = "free_x",
      space = "free_x"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        ATL_N_UPW = "#4575b4",
        ATL_S_UPW_ALM = "#91bfdb",
        MED = "#8DD3C7"
      ),
      name = "Oceanographic zone"
    ) +
    ggplot2::labs(x = "Genetic cluster", y = "% of localities") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.spacing = grid::unit(0.4, "lines"),
      legend.position = "right"
    )

  ggplot2::ggsave(
    filename = output_file,
    plot = current_plot,
    width = 18,
    height = 11,
    dpi = 300
  )

  invisible(current_plot)
}

# Run the complete all-haplotype clustering workflow.
# Writes all legacy summary tables plus optional maps.
run_all_haplotype_clustering <- function(
  otus_root,
  output_dir,
  label = "metazoa",
  metadata_file = NULL,
  sigma_grid = c(0.005, 0.010, 0.015),
  k_moran = 5,
  seed = 42,
  write_maps = TRUE
) {
  label <- normalize_output_label(label)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  network_data <- load_haplotype_network_data(otus_root)
  edges_with_info <- build_edges_with_info(network_data$edges, network_data$points)
  if (nrow(edges_with_info) == 0) {
    stop("No valid inter-locality, within-OTU edges were found.")
  }

  clustering_outputs <- run_all_haplotype_partitions(
    edges_with_info,
    sigma_grid = sigma_grid,
    seed = seed
  )
  comp_all_coords <- attach_coordinates_to_clusters(
    clustering_outputs$comp_all,
    edges_with_info
  )

  write.csv(edges_with_info, file.path(output_dir, "edges_with_info.csv"), row.names = FALSE)
  write.csv(network_data$points, file.path(output_dir, "points_all.csv"), row.names = FALSE)
  write.csv(
    clustering_outputs$comp_all,
    file.path(output_dir, paste0("clusters_", label, ".csv")),
    row.names = FALSE
  )
  write.csv(
    comp_all_coords,
    file.path(output_dir, paste0("clusters_", label, "_with_coordinates.csv")),
    row.names = FALSE
  )
  write.csv(
    clustering_outputs$ari_table,
    file.path(output_dir, paste0("ari_", label, ".csv")),
    row.names = FALSE
  )
  write.csv(
    clustering_outputs$mod_summary,
    file.path(output_dir, paste0("modularity_summary_", label, ".csv")),
    row.names = FALSE
  )
  write.csv(
    clustering_outputs$mod_results,
    file.path(output_dir, paste0("modularity_results_", label, ".csv")),
    row.names = FALSE
  )
  write.csv(
    clustering_outputs$rbf_grid,
    file.path(output_dir, paste0("rbf_sigma_grid_", label, ".csv")),
    row.names = FALSE
  )

  moran_outputs <- list(
    inverse = moran_table_by_cluster(clustering_outputs$results$inverse$coords, "Cluster_inv", k = k_moran),
    count = moran_table_by_cluster(clustering_outputs$results$count$coords, "Cluster_cnt", k = k_moran),
    linear = moran_table_by_cluster(clustering_outputs$results$linear$coords, "Cluster_linear", k = k_moran),
    rbf = moran_table_by_cluster(clustering_outputs$results$rbf$coords, "Cluster_rbf", k = k_moran)
  )
  write.csv(moran_outputs$inverse, file.path(output_dir, paste0("moran_inverse_", label, ".csv")), row.names = FALSE)
  write.csv(moran_outputs$count, file.path(output_dir, paste0("moran_count_", label, ".csv")), row.names = FALSE)
  write.csv(moran_outputs$linear, file.path(output_dir, paste0("moran_linear_", label, ".csv")), row.names = FALSE)
  write.csv(moran_outputs$rbf, file.path(output_dir, paste0("moran_rbf_", label, ".csv")), row.names = FALSE)

  ocean_df <- read_ocean_metadata(metadata_file)
  if (!is.null(ocean_df)) {
    write.csv(
      ocean_share_by_cluster(clustering_outputs$results$inverse$coords, "Cluster_inv", ocean_df),
      file.path(output_dir, paste0("tabla_coords_inverse_", label, ".csv")),
      row.names = FALSE
    )
    write.csv(
      ocean_share_by_cluster(clustering_outputs$results$count$coords, "Cluster_cnt", ocean_df),
      file.path(output_dir, paste0("tabla_coords_count_", label, ".csv")),
      row.names = FALSE
    )
    write.csv(
      ocean_share_by_cluster(clustering_outputs$results$linear$coords, "Cluster_linear", ocean_df),
      file.path(output_dir, paste0("tabla_coords_linear_", label, ".csv")),
      row.names = FALSE
    )
    write.csv(
      ocean_share_by_cluster(clustering_outputs$results$rbf$coords, "Cluster_rbf", ocean_df),
      file.path(output_dir, paste0("tabla_coords_rbf_", label, ".csv")),
      row.names = FALSE
    )
  }

  zones_df <- build_current_zones(clustering_outputs$results$inverse$coords)
  message("Current-zone assignment summary:")
  print(table(is.na(zones_df$Zone)))
  write.csv(
    zones_df,
    file.path(output_dir, "locality_zones_currents_all.csv"),
    row.names = FALSE
  )

  write.csv(
    zone_share_by_cluster(clustering_outputs$results$inverse$coords, "Cluster_inv", zones_df),
    file.path(output_dir, paste0("zone_share_inverse_", label, ".csv")),
    row.names = FALSE
  )
  write.csv(
    zone_share_by_cluster(clustering_outputs$results$count$coords, "Cluster_cnt", zones_df),
    file.path(output_dir, paste0("zone_share_count_", label, ".csv")),
    row.names = FALSE
  )
  write.csv(
    zone_share_by_cluster(clustering_outputs$results$linear$coords, "Cluster_linear", zones_df),
    file.path(output_dir, paste0("zone_share_linear_", label, ".csv")),
    row.names = FALSE
  )
  write.csv(
    zone_share_by_cluster(clustering_outputs$results$rbf$coords, "Cluster_rbf", zones_df),
    file.path(output_dir, paste0("zone_share_rbf_", label, ".csv")),
    row.names = FALSE
  )

  if (write_maps) {
    write_all_haplotype_maps(comp_all_coords, output_dir, label)
  }

  message("All-haplotype clustering outputs written to: ", output_dir)
  invisible(list(edges_with_info = edges_with_info, clustering_outputs = clustering_outputs))
}

# Run the complete SAME vs DIFF current analysis workflow.
# Writes clustering comparisons, Moran tables, zones, and distance summaries.
run_same_diff_current_analysis <- function(
  otus_root,
  output_dir,
  label = "metazoa",
  all_output_dir = NULL,
  metadata_file = NULL,
  sigma_grid = c(0.005, 0.010, 0.015),
  k_moran = 5,
  seed = 42,
  write_current_figure = TRUE
) {
  label <- normalize_output_label(label)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  network_data <- load_haplotype_network_data(otus_root)
  edges_with_info <- build_edges_with_info(network_data$edges, network_data$points)
  if (nrow(edges_with_info) == 0) {
    stop("No valid inter-locality, within-OTU edges were found.")
  }

  write.csv(edges_with_info, file.path(output_dir, "edges_with_info.csv"), row.names = FALSE)

  locality_coords_df <- get_locality_coordinates(edges_with_info)
  zones_df <- build_current_zones(locality_coords_df)
  message("Current-zone assignment summary:")
  print(table(is.na(zones_df$Zone)))
  write.csv(zones_df, file.path(output_dir, "locality_zones_currents.csv"), row.names = FALSE)

  same_diff_outputs <- run_same_diff_partitions(
    edges_with_info,
    sigma_grid = sigma_grid,
    k_moran = k_moran,
    seed = seed
  )

  write.csv(
    same_diff_outputs$mod_summary,
    file.path(output_dir, paste0("modularity_same_diff_summary_", label, ".csv")),
    row.names = FALSE
  )
  write.csv(
    same_diff_outputs$ari_same,
    file.path(output_dir, paste0("ari_same_within_methods_", label, ".csv")),
    row.names = FALSE
  )
  write.csv(
    same_diff_outputs$ari_diff,
    file.path(output_dir, paste0("ari_diff_within_methods_", label, ".csv")),
    row.names = FALSE
  )
  write.csv(
    same_diff_outputs$ari_same_vs_diff,
    file.path(output_dir, paste0("ari_same_vs_diff_by_method_", label, ".csv")),
    row.names = FALSE
  )

  for (group_name in c("same", "diff")) {
    for (method_name in c("inverse", "count", "linear", "rbf")) {
      result_obj <- same_diff_outputs[[group_name]][[method_name]]
      write.csv(
        result_obj$moran,
        file.path(output_dir, paste0("moran_", group_name, "_", method_name, "_", label, ".csv")),
        row.names = FALSE
      )
      write.csv(
        result_obj$locality_coords,
        file.path(output_dir, paste0("locality_coords_", group_name, "_", method_name, ".csv")),
        row.names = FALSE
      )
    }
  }

  ocean_df <- read_ocean_metadata(metadata_file)
  if (!is.null(ocean_df)) {
    for (group_name in c("same", "diff")) {
      for (method_name in c("inverse", "count", "linear", "rbf")) {
        ocean_share <- ocean_share_by_cluster(
          same_diff_outputs[[group_name]][[method_name]]$locality_coords,
          "Cluster",
          ocean_df
        )
        write.csv(
          ocean_share,
          file.path(output_dir, paste0("ocean_share_", group_name, "_", method_name, "_", label, ".csv")),
          row.names = FALSE
        )
      }
    }
  }

  for (group_name in c("same", "diff")) {
    for (method_name in c("inverse", "count", "linear", "rbf")) {
      zone_share <- zone_share_by_cluster(
        same_diff_outputs[[group_name]][[method_name]]$locality_coords,
        "Cluster",
        zones_df
      )
      write.csv(
        zone_share,
        file.path(output_dir, paste0("zone_share_", group_name, "_", method_name, "_", label, ".csv")),
        row.names = FALSE
      )
    }
  }

  geo_outputs <- compute_geographical_distance_outputs(edges_with_info)
  write.csv(
    geo_outputs$med_geo_by_otu,
    file.path(output_dir, paste0("medianas_geo_por_OTU_", label, ".csv")),
    row.names = FALSE
  )
  write.csv(
    geo_outputs$resumen_geo_otu,
    file.path(output_dir, paste0("resumen_medianas_geo_por_OTU_", label, ".csv")),
    row.names = FALSE
  )
  write.csv(
    geo_outputs$stats_summary,
    file.path(output_dir, "estadisticos_dist_geograficas.csv"),
    row.names = FALSE
  )
  plot_geographical_distances(
    geo_outputs,
    file.path(output_dir, paste0("p_dist_", label, ".pdf"))
  )

  if (write_current_figure) {
    zone_dirs <- output_dir
    if (!is.null(all_output_dir)) {
      zone_dirs <- c(zone_dirs, all_output_dir)
    }
    plot_current_zone_summary(
      zone_dirs,
      file.path(output_dir, "Figure_currents_all_methods_ALL_same_diff_ALLgroups.png")
    )
  }

  message("SAME vs DIFF current-analysis outputs written to: ", output_dir)
  invisible(list(edges_with_info = edges_with_info, same_diff_outputs = same_diff_outputs))
}
