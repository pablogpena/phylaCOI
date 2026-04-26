# -*- coding: utf-8 -*-
# Utilities for computing nucleotide diversity, abundance, and haplotype network connections.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(geosphere)
  library(pegas)
  library(ape)
  library(leaflet)
  library(scales)
  library(htmlwidgets)
  library(mgcv)
  library(Biostrings)
  library(tibble)
})

# Locate required inputs for a phylum and validate they exist.
# Returns a named list of paths or NULL if missing.
resolve_input_files <- function(phylum_name, abundance_root, otus_root) {
  abundance_dir <- file.path(abundance_root, phylum_name)
  otus_dir <- file.path(otus_root, phylum_name)

  abundances_file <- file.path(abundance_dir, "abundances.csv")
  informative_otus_file <- file.path(otus_dir, "informative_OTUs.txt")
  alignment_file <- file.path(abundance_dir, "aligned_sequences_mafft.fasta")
  otu_mapping_file <- file.path(otus_dir, "otus_mapping.txt")

  if (!file.exists(abundances_file)) {
    warning("Missing abundance file in: ", phylum_name)
    return(NULL)
  }
  if (!file.exists(informative_otus_file)) {
    warning("Missing informative OTU list in: ", phylum_name)
    return(NULL)
  }
  if (!file.exists(alignment_file) || !file.exists(otu_mapping_file)) {
    warning("Missing alignment or OTU mapping in: ", phylum_name)
    return(NULL)
  }

  list(
    abundances_file = abundances_file,
    informative_otus_file = informative_otus_file,
    alignment_file = alignment_file,
    otu_mapping_file = otu_mapping_file
  )
}

# Ensure the Abundance column exists and is numeric.
# Stops early if the column is missing.
normalize_abundance_column <- function(abundances) {
  # Ensure the refactored abundance column exists.
  if (!"Abundance" %in% names(abundances)) {
    stop("Missing abundance column: Abundance.")
  }

  abundances$Abundance <- as.numeric(abundances$Abundance)
  abundances
}

# Assign haplotypes within each OTU using the alignment.
# Skips large OTUs and returns filtered abundances and tracking lists.
assign_haplotypes <- function(abundances, alignment_dnabin, max_sequences_per_otu = 500) {
  # Assign haplotypes within each OTU, skipping OTUs that are too large.
  abundances$Haplotype <- NA_integer_

  otus_with_network <- character()
  otus_skipped <- character()

  if (nrow(abundances) == 0) {
    return(list(
      abundances = abundances,
      otus_with_network = otus_with_network,
      otus_skipped = otus_skipped
    ))
  }

  for (otu in unique(abundances$OTU)) {
    otu_ids <- abundances$UniqueID[abundances$OTU == otu]

    if (length(otu_ids) > max_sequences_per_otu) {
      otus_skipped <- c(otus_skipped, otu)
      next
    }

    otu_alignment <- alignment_dnabin[names(alignment_dnabin) %in% otu_ids]
    if (length(otu_alignment) > 1) {
      haps <- pegas::haplotype(otu_alignment)
      idx <- attr(haps, "index")
      seq_names <- names(otu_alignment)

      for (h in seq_along(idx)) {
        abundances$Haplotype[abundances$UniqueID %in% seq_names[idx[[h]]]] <- h
      }

      otus_with_network <- c(otus_with_network, otu)
    }
  }

  abundances <- dplyr::filter(abundances, OTU %in% otus_with_network)

  list(
    abundances = abundances,
    otus_with_network = otus_with_network,
    otus_skipped = otus_skipped
  )
}

# Cluster points by geographic radius using haversine distance.
# Adds a cluster label per OTU and returns the updated data.
cluster_points_by_distance <- function(df, threshold_m) {
  # Group sequences into spatial clusters using a fixed radius (meters).
  if (nrow(df) == 0) {
    return(df)
  }

  coords <- df %>% dplyr::select(lon, lat) %>% as.matrix()
  dist_matrix <- geosphere::distm(coords, fun = distHaversine)

  clusters <- rep(NA_integer_, nrow(coords))
  cluster_id <- 1L

  for (i in seq_len(nrow(dist_matrix))) {
    if (is.na(clusters[i])) {
      nearby <- which(dist_matrix[i, ] <= threshold_m)
      clusters[nearby] <- cluster_id
      cluster_id <- cluster_id + 1L
    }
  }

  df$cluster <- paste0(df$OTU[1], "_", clusters)
  df
}

# Compute abundance summaries, distances, and model fits per OTU.
# Returns per-cluster data plus metric summaries.
compute_abundance_metrics <- function(abundances_grouped) {
  empty_results <- tibble::tibble(
    OTU = character(),
    cluster = character(),
    total_abundance = numeric(),
    lat = numeric(),
    lon = numeric(),
    log_abundance = numeric()
  )
  empty_distances <- tibble::tibble(
    OTU = character(),
    cluster = character(),
    total_abundance = numeric(),
    lat = numeric(),
    lon = numeric(),
    log_abundance = numeric(),
    max_log_abundance = numeric(),
    center_lat = numeric(),
    center_lon = numeric(),
    distance_km = numeric()
  )
  empty_metrics <- tibble::tibble(
    OTU = character(),
    max_log_abundance = numeric(),
    mean_log_abundance = numeric(),
    max_distance = numeric(),
    lm_slope = numeric(),
    spearman_rho = numeric(),
    nls_a = numeric(),
    nls_b = numeric(),
    gam_r2 = numeric()
  )

  if (nrow(abundances_grouped) == 0) {
    return(list(
      abundance_results = empty_results,
      abundance_distances = empty_distances,
      abundance_metrics = empty_metrics
    ))
  }

  abundance_results <- abundances_grouped %>%
    dplyr::group_by(OTU, cluster) %>%
    dplyr::summarise(
      total_abundance = sum(Abundance, na.rm = TRUE),
      lat = mean(lat),
      lon = mean(lon),
      .groups = "drop"
    ) %>%
    dplyr::mutate(log_abundance = log1p(total_abundance))

  abundance_distances <- abundance_results %>%
    dplyr::group_by(OTU) %>%
    dplyr::mutate(
      max_log_abundance = max(log_abundance, na.rm = TRUE),
      center_lat = mean(lat[log_abundance == max_log_abundance], na.rm = TRUE),
      center_lon = mean(lon[log_abundance == max_log_abundance], na.rm = TRUE),
      distance_km = geosphere::distHaversine(
        matrix(c(lon, lat), ncol = 2),
        matrix(c(center_lon, center_lat), ncol = 2)
      ) / 1000
    ) %>%
    dplyr::ungroup()

  abundance_model_data <- abundance_distances %>%
    dplyr::group_by(OTU) %>%
    dplyr::filter(dplyr::n() >= 3) %>%
    dplyr::ungroup()

  abundance_metrics <- empty_metrics
  if (nrow(abundance_model_data) > 0) {
    abundance_metrics <- abundance_model_data %>%
      dplyr::group_by(OTU) %>%
      dplyr::group_modify(~{
        df <- .x
        lm_model <- lm(log_abundance ~ distance_km, data = df)
        cor_spear <- suppressWarnings(cor(df$log_abundance, df$distance_km, method = "spearman"))
        nls_model <- tryCatch({
          nls(log_abundance ~ a * exp(-b * distance_km),
              data = df, start = list(a = max(df$log_abundance), b = 0.01))
        }, error = function(e) NULL)
        gam_model <- tryCatch({
          mgcv::gam(log_abundance ~ s(distance_km), data = df)
        }, error = function(e) NULL)

        tibble::tibble(
          max_log_abundance = max(df$log_abundance, na.rm = TRUE),
          mean_log_abundance = mean(df$log_abundance, na.rm = TRUE),
          max_distance = max(df$distance_km, na.rm = TRUE),
          lm_slope = coef(lm_model)[2],
          spearman_rho = cor_spear,
          nls_a = if (!is.null(nls_model)) coef(nls_model)["a"] else NA,
          nls_b = if (!is.null(nls_model)) coef(nls_model)["b"] else NA,
          gam_r2 = if (!is.null(gam_model)) summary(gam_model)$r.sq else NA
        )
      }) %>%
      dplyr::ungroup()
  }

  list(
    abundance_results = abundance_results,
    abundance_distances = abundance_distances,
    abundance_metrics = abundance_metrics
  )
}

# Build haplotype network connections and per-cluster connectivity metrics.
# Returns connection distances plus edge/point tables.
compute_connections_metrics <- function(abundances_grouped, alignment_dnabin, otus_with_network) {
  empty_connections <- tibble::tibble(
    OTU = character(),
    cluster = character(),
    lat = numeric(),
    lon = numeric(),
    n_connections = integer(),
    log_connections = numeric(),
    max_conn = numeric(),
    center_lat = numeric(),
    center_lon = numeric(),
    distance_km_connections = numeric()
  )
  empty_edges <- tibble::tibble(
    from = character(),
    to = character(),
    OTU_ID = character(),
    group = character(),
    distancia_genetica = numeric(),
    x = numeric(),
    y = numeric(),
    xend = numeric(),
    yend = numeric()
  )
  empty_points <- tibble::tibble(
    UniqueID = character(),
    ID = character(),
    Localities = character(),
    OTU_ID = character(),
    Haplotype = integer(),
    cluster = character(),
    lat = numeric(),
    lon = numeric()
  )

  if (length(otus_with_network) == 0 || nrow(abundances_grouped) == 0) {
    return(list(
      connections_distances = empty_connections,
      network_edges = empty_edges,
      network_points = empty_points
    ))
  }

  connections_by_cluster <- list()
  network_edges_by_otu <- list()
  network_points <- abundances_grouped %>%
    dplyr::filter(OTU %in% otus_with_network) %>%
    dplyr::transmute(
      UniqueID = UniqueID,
      ID = ID,
      Localities = Localities,
      OTU_ID = OTU,
      Haplotype = Haplotype,
      cluster = cluster,
      lat = lat,
      lon = lon
    ) %>%
    dplyr::distinct()

  for (otu in otus_with_network) {
    message("Calculating connections for OTU: ", otu)

    otu_data <- dplyr::filter(abundances_grouped, OTU == otu)
    ids <- otu_data$UniqueID
    aln <- alignment_dnabin[names(alignment_dnabin) %in% ids]

    if (length(aln) < 2) {
      message("  -> Fewer than 2 sequences. Skipping.")
      next
    }

    haps <- pegas::haplotype(aln)
    if (nrow(haps) < 2) {
      message("  -> Single haplotype only. Skipping.")
      next
    }

    net <- pegas::haploNet(haps)
    labels <- attr(net, "labels")
    haplo_to_seq <- lapply(attr(haps, "index"), function(idx) names(aln)[idx])
    names(haplo_to_seq) <- labels

    edges <- apply(as.matrix(net), 1, function(edge) {
      from <- haplo_to_seq[[labels[edge[1]]]]
      to <- haplo_to_seq[[labels[edge[2]]]]
      base::expand.grid(from = from, to = to, stringsAsFactors = FALSE)
    })
    seq_edges_df <- do.call(rbind, edges)
    if (is.null(seq_edges_df)) {
      seq_edges_df <- tibble::tibble(from = character(), to = character())
    }
    seq_edges_df <- seq_edges_df %>%
      dplyr::mutate(group = "diff")

    haplo_assignments <- tibble::tibble(
      UniqueID = unlist(haplo_to_seq),
      haplotype = rep(names(haplo_to_seq), lengths(haplo_to_seq))
    )

    haplo_coords <- dplyr::left_join(haplo_assignments, otu_data, by = "UniqueID")

    haplo_edges_df <- haplo_coords %>%
      dplyr::group_by(haplotype) %>%
      dplyr::filter(dplyr::n_distinct(cluster) > 1) %>%
      tidyr::expand(from = UniqueID, to = UniqueID) %>%
      dplyr::filter(from != to) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(group = "same")

    all_edges_df <- dplyr::bind_rows(seq_edges_df, haplo_edges_df) %>%
      dplyr::distinct()

    coord_map <- otu_data %>%
      dplyr::select(ID = UniqueID, cluster, lat, lon)
    # Helper to fetch a coordinate attribute for a given sequence ID.
    get_attr <- function(id, attr) coord_map[[attr]][coord_map$ID == id][1]

    genetic_dist_mat <- as.matrix(ape::dist.dna(aln, model = "K80"))
    # Helper to fetch genetic distance for a sequence pair.
    get_genetic_distance <- function(from_id, to_id) {
      if (!from_id %in% rownames(genetic_dist_mat) || !to_id %in% colnames(genetic_dist_mat)) {
        return(NA_real_)
      }
      genetic_dist_mat[from_id, to_id]
    }

    all_edges_df <- all_edges_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        from_cluster = get_attr(from, "cluster"),
        to_cluster = get_attr(to, "cluster"),
        from_lat = get_attr(from, "lat"),
        from_lon = get_attr(from, "lon"),
        to_lat = get_attr(to, "lat"),
        to_lon = get_attr(to, "lon"),
        distancia_genetica = get_genetic_distance(from, to)
      ) %>%
      dplyr::ungroup()

    all_edges_df <- dplyr::filter(all_edges_df, from_cluster != to_cluster)

    network_edges_by_otu[[as.character(otu)]] <- all_edges_df %>%
      dplyr::transmute(
        from = from,
        to = to,
        OTU_ID = otu,
        group = group,
        distancia_genetica = distancia_genetica,
        x = from_lon,
        y = from_lat,
        xend = to_lon,
        yend = to_lat
      )

    all_clusters <- otu_data %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(
        lat = mean(lat, na.rm = TRUE),
        lon = mean(lon, na.rm = TRUE),
        .groups = "drop"
      )

    if (nrow(all_edges_df) > 0) {
      edge_counts <- all_edges_df %>%
        dplyr::group_by(from_cluster) %>%
        dplyr::summarise(
          n_connections = dplyr::n(),
          lat = mean(from_lat, na.rm = TRUE),
          lon = mean(from_lon, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::rename(cluster = from_cluster)
    } else {
      edge_counts <- tibble::tibble(
        cluster = character(),
        n_connections = integer(),
        lat = numeric(),
        lon = numeric()
      )
    }

    edge_counts_full <- dplyr::full_join(all_clusters, edge_counts,
                                         by = c("cluster", "lat", "lon")) %>%
      dplyr::mutate(
        n_connections = dplyr::coalesce(n_connections, 0),
        OTU = otu,
        log_connections = log1p(n_connections)
      )

    connections_by_cluster[[as.character(otu)]] <- edge_counts_full
  }

  connections_raw <- dplyr::bind_rows(connections_by_cluster)
  network_edges <- dplyr::bind_rows(network_edges_by_otu)
  if (nrow(network_edges) == 0) {
    network_edges <- empty_edges
  }
  if (nrow(network_points) == 0) {
    network_points <- empty_points
  }

  if (nrow(connections_raw) == 0) {
    return(list(
      connections_distances = empty_connections,
      network_edges = network_edges,
      network_points = network_points
    ))
  }

  connections_summary <- connections_raw %>%
    dplyr::group_by(OTU, cluster) %>%
    dplyr::summarise(
      lat = mean(lat, na.rm = TRUE),
      lon = mean(lon, na.rm = TRUE),
      n_connections = sum(n_connections, na.rm = TRUE),
      log_connections = log1p(n_connections),
      .groups = "drop"
    )

  connections_distances <- connections_summary %>%
    dplyr::group_by(OTU) %>%
    dplyr::mutate(
      max_conn = max(n_connections, na.rm = TRUE),
      center_lat = lat[n_connections == max_conn][1],
      center_lon = lon[n_connections == max_conn][1],
      distance_km_connections = geosphere::distHaversine(
        matrix(c(lon, lat), ncol = 2),
        matrix(c(center_lon, center_lat), ncol = 2)
      ) / 1000
    ) %>%
    dplyr::ungroup()

  list(
    connections_distances = connections_distances,
    network_edges = network_edges,
    network_points = network_points
  )
}

# Pick a center cluster using diversity, connections, and spatial distance.
# Returns the index of the selected cluster.
choose_center_cluster <- function(diversity, n_connections, lon, lat) {
  if (all(is.na(diversity))) {
    return(1)
  }

  idx_max <- which(diversity == max(diversity, na.rm = TRUE))
  if (length(idx_max) == 1) {
    return(idx_max)
  }

  n_connections[is.na(n_connections)] <- 0
  conn_values <- n_connections[idx_max]
  idx_best <- idx_max[which.max(conn_values)]

  if (sum(conn_values == max(conn_values)) > 1) {
    idx_candidates <- idx_max[conn_values == max(conn_values)]
    mean_distances <- sapply(idx_candidates, function(i) {
      dist_vals <- geosphere::distHaversine(
        matrix(c(lon, lat), ncol = 2),
        matrix(c(lon[i], lat[i]), ncol = 2)
      )
      mean(dist_vals, na.rm = TRUE)
    })
    idx_best <- idx_candidates[which.min(mean_distances)]
  }

  idx_best
}

# Compute nucleotide diversity per cluster and distance-based summaries.
# Returns per-cluster data and model metrics.
compute_diversity_metrics <- function(abundances_grouped, alignment_dnabin, connections_distances) {
  empty_results <- tibble::tibble(
    OTU = character(),
    cluster = character(),
    seqs = list(),
    lat = numeric(),
    lon = numeric(),
    diversity = numeric(),
    n_connections = integer()
  )
  empty_distances <- tibble::tibble(
    OTU = character(),
    cluster = character(),
    seqs = list(),
    lat = numeric(),
    lon = numeric(),
    diversity = numeric(),
    n_connections = integer(),
    max_div = numeric(),
    center_cluster = integer(),
    center_lat = numeric(),
    center_lon = numeric(),
    distance_km = numeric()
  )
  empty_metrics <- tibble::tibble(
    OTU = character(),
    max_diversity = numeric(),
    mean_diversity = numeric(),
    max_distance = numeric(),
    lm_slope = numeric(),
    spearman_rho = numeric(),
    nls_a = numeric(),
    nls_b = numeric(),
    gam_r2 = numeric()
  )

  if (nrow(abundances_grouped) == 0) {
    return(list(
      diversity_results = empty_results,
      diversity_distances = empty_distances,
      diversity_metrics = empty_metrics
    ))
  }

  diversity_results <- abundances_grouped %>%
    dplyr::group_by(OTU, cluster) %>%
    dplyr::summarise(
      seqs = list(UniqueID),
      lat = mean(lat),
      lon = mean(lon),
      .groups = "drop"
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      diversity = {
        selected <- names(alignment_dnabin) %in% seqs
        aln <- alignment_dnabin[selected]
        if (length(aln) > 1) pegas::nuc.div(aln) else 0
      }
    ) %>%
    dplyr::ungroup()

  diversity_results <- diversity_results %>%
    dplyr::left_join(
      connections_distances %>% dplyr::select(OTU, cluster, n_connections),
      by = c("OTU", "cluster")
    )

  diversity_distances <- diversity_results %>%
    dplyr::group_by(OTU) %>%
    dplyr::mutate(
      n_connections = dplyr::coalesce(n_connections, 0),
      max_div = max(diversity, na.rm = TRUE),
      center_cluster = choose_center_cluster(diversity, n_connections, lon, lat),
      center_lat = lat[center_cluster],
      center_lon = lon[center_cluster],
      distance_km = geosphere::distHaversine(
        matrix(c(lon, lat), ncol = 2),
        matrix(c(center_lon, center_lat), ncol = 2)
      ) / 1000
    ) %>%
    dplyr::ungroup()

  diversity_model_data <- diversity_distances %>%
    dplyr::group_by(OTU) %>%
    dplyr::filter(dplyr::n() >= 3) %>%
    dplyr::ungroup()

  diversity_metrics <- empty_metrics
  if (nrow(diversity_model_data) > 0) {
    diversity_metrics <- diversity_model_data %>%
      dplyr::group_by(OTU) %>%
      dplyr::group_modify(~{
        df <- .x
        lm_model <- lm(diversity ~ distance_km, data = df)
        cor_spear <- suppressWarnings(cor(df$diversity, df$distance_km, method = "spearman"))
        nls_model <- tryCatch({
          nls(diversity ~ a * exp(-b * distance_km),
              data = df, start = list(a = max(df$diversity), b = 0.01))
        }, error = function(e) NULL)
        gam_model <- tryCatch({
          mgcv::gam(diversity ~ s(distance_km), data = df)
        }, error = function(e) NULL)

        tibble::tibble(
          max_diversity = max(df$diversity, na.rm = TRUE),
          mean_diversity = mean(df$diversity, na.rm = TRUE),
          max_distance = max(df$distance_km, na.rm = TRUE),
          lm_slope = coef(lm_model)[2],
          spearman_rho = cor_spear,
          nls_a = if (!is.null(nls_model)) coef(nls_model)["a"] else NA,
          nls_b = if (!is.null(nls_model)) coef(nls_model)["b"] else NA,
          gam_r2 = if (!is.null(gam_model)) summary(gam_model)$r.sq else NA
        )
      }) %>%
      dplyr::ungroup()
  }

  list(
    diversity_results = diversity_results,
    diversity_distances = diversity_distances,
    diversity_metrics = diversity_metrics
  )
}

# Create leaflet maps for diversity, abundance, and connections.
# Writes HTML files to the output directory.
build_maps <- function(output_dir, diversity_results, abundance_distances, connections_distances) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (nrow(diversity_results) > 0) {
    pal_div <- leaflet::colorNumeric("viridis", domain = diversity_results$diversity)
    map_div <- leaflet::leaflet() %>%
      leaflet::addProviderTiles(leaflet::providers$CartoDB.Positron)

    for (otu in unique(diversity_results$OTU)) {
      map_div <- map_div %>%
        leaflet::addCircleMarkers(
          data = dplyr::filter(diversity_results, OTU == otu),
          lng = ~lon, lat = ~lat,
          radius = ~scales::rescale(diversity, to = c(4, 12)),
          color = ~pal_div(diversity),
          stroke = TRUE, fillOpacity = 0.8,
          label = ~paste0(
            "OTU: ", OTU,
            "<br>Cluster: ", cluster,
            "<br>Diversity: ", round(diversity, 5)
          ),
          group = paste0("OTU ", otu)
        )
    }

    map_div <- map_div %>%
      leaflet::addLegend(
        "bottomright", pal = pal_div,
        values = diversity_results$diversity,
        title = "Nucleotide diversity"
      ) %>%
      leaflet::addLayersControl(
        overlayGroups = paste0("OTU ", unique(diversity_results$OTU)),
        options = leaflet::layersControlOptions(collapsed = FALSE)
      )

    htmlwidgets::saveWidget(
      map_div,
      file.path(output_dir, "otu_nucleotide_diversity_map.html"),
      selfcontained = TRUE
    )
  } else {
    message("No diversity data available for maps.")
  }

  if (nrow(abundance_distances) > 0) {
    pal_abun <- leaflet::colorNumeric("viridis", domain = abundance_distances$log_abundance)
    map_abun <- leaflet::leaflet() %>%
      leaflet::addProviderTiles(leaflet::providers$CartoDB.Positron)

    for (otu in unique(abundance_distances$OTU)) {
      map_abun <- map_abun %>%
        leaflet::addCircleMarkers(
          data = dplyr::filter(abundance_distances, OTU == otu),
          lng = ~lon, lat = ~lat,
          radius = ~scales::rescale(log_abundance, to = c(4, 12)),
          color = ~pal_abun(log_abundance),
          stroke = TRUE, fillOpacity = 0.8,
          label = ~paste0(
            "OTU: ", OTU,
            "<br>Cluster: ", cluster,
            "<br>Log(Abundance): ", round(log_abundance, 3)
          ),
          group = paste0("OTU ", otu)
        )
    }

    map_abun <- map_abun %>%
      leaflet::addLegend(
        "bottomright", pal = pal_abun,
        values = abundance_distances$log_abundance,
        title = "Log(Abundance)"
      ) %>%
      leaflet::addLayersControl(
        overlayGroups = paste0("OTU ", unique(abundance_distances$OTU)),
        options = leaflet::layersControlOptions(collapsed = FALSE)
      )

    htmlwidgets::saveWidget(
      map_abun,
      file.path(output_dir, "otu_log_abundance_map.html"),
      selfcontained = TRUE
    )
  } else {
    message("No abundance data available for maps.")
  }

  if (nrow(connections_distances) > 0) {
    pal_conn <- leaflet::colorNumeric("viridis", domain = connections_distances$log_connections)
    map_conn <- leaflet::leaflet() %>%
      leaflet::addProviderTiles(leaflet::providers$CartoDB.Positron)

    for (otu in unique(connections_distances$OTU)) {
      map_conn <- map_conn %>%
        leaflet::addCircleMarkers(
          data = dplyr::filter(connections_distances, OTU == otu),
          lng = ~lon, lat = ~lat,
          radius = ~scales::rescale(log_connections, to = c(4, 12)),
          color = ~pal_conn(log_connections),
          stroke = TRUE, fillOpacity = 0.8,
          label = ~paste0(
            "OTU: ", OTU,
            "<br>Cluster: ", cluster,
            "<br>Log(Connections): ", round(log_connections, 3)
          ),
          group = paste0("OTU ", otu)
        )
    }

    map_conn <- map_conn %>%
      leaflet::addLegend(
        "bottomright", pal = pal_conn,
        values = connections_distances$log_connections,
        title = "Log(Connections)"
      ) %>%
      leaflet::addLayersControl(
        overlayGroups = paste0("OTU ", unique(connections_distances$OTU)),
        options = leaflet::layersControlOptions(collapsed = FALSE)
      )

    htmlwidgets::saveWidget(
      map_conn,
      file.path(output_dir, "otu_log_connections_map.html"),
      selfcontained = TRUE
    )
  } else {
    message("No connection data available for maps.")
  }
}

# Run the full OTU metrics pipeline for one phylum.
# Reads inputs, computes metrics, writes CSVs, and optional maps.
process_phylum <- function(
  phylum_name,
  abundance_root,
  otus_root,
  cluster_radius_m = 5000,
  max_sequences_per_otu = 500,
  write_maps = TRUE
) {
  message("\n========================")
  message("Processing phylum: ", phylum_name)
  message("========================")

  input_files <- resolve_input_files(phylum_name, abundance_root, otus_root)
  if (is.null(input_files)) {
    return(NULL)
  }

  abundances <- read.csv(input_files$abundances_file, stringsAsFactors = FALSE)
  abundances <- normalize_abundance_column(abundances)

  if (!"UniqueID" %in% names(abundances)) {
    warning("Missing UniqueID column in: ", phylum_name)
    return(NULL)
  }
  if (!all(c("lat", "lon") %in% names(abundances))) {
    warning("Missing lat/lon columns in: ", phylum_name)
    return(NULL)
  }

  otu_map <- read.table(
    input_files$otu_mapping_file,
    col.names = c("OTU", "UniqueID"),
    stringsAsFactors = FALSE
  )
  informative_otus <- scan(input_files$informative_otus_file, what = character(), quiet = TRUE)

  if (length(informative_otus) == 0) {
    message("No informative OTUs found for: ", phylum_name)
    return(NULL)
  }

  abundances$OTU <- otu_map$OTU[match(abundances$UniqueID, otu_map$UniqueID)]
  abundances_filtered <- dplyr::filter(abundances, OTU %in% informative_otus)

  if (nrow(abundances_filtered) == 0) {
    message("No informative OTUs after mapping for: ", phylum_name)
    return(NULL)
  }

  alignment <- Biostrings::readDNAStringSet(input_files$alignment_file)
  alignment_dnabin <- ape::as.DNAbin(alignment)

  haplotype_result <- assign_haplotypes(
    abundances_filtered,
    alignment_dnabin,
    max_sequences_per_otu = max_sequences_per_otu
  )
  abundances_filtered <- haplotype_result$abundances
  otus_with_network <- haplotype_result$otus_with_network

  if (length(haplotype_result$otus_skipped) > 0) {
    message(
      "Skipped OTUs with more than ", max_sequences_per_otu, " sequences: ",
      paste(haplotype_result$otus_skipped, collapse = ", ")
    )
  }

  if (nrow(abundances_filtered) == 0) {
    message("No OTUs with haplotype networks for: ", phylum_name)
    return(NULL)
  }

  abundances_grouped <- abundances_filtered %>%
    dplyr::group_by(OTU) %>%
    dplyr::group_modify(~ cluster_points_by_distance(.x, cluster_radius_m)) %>%
    dplyr::ungroup()

  abundance_outputs <- compute_abundance_metrics(abundances_grouped)
  connection_outputs <- compute_connections_metrics(
    abundances_grouped,
    alignment_dnabin,
    otus_with_network
  )
  diversity_outputs <- compute_diversity_metrics(
    abundances_grouped,
    alignment_dnabin,
    connection_outputs$connections_distances
  )

  diversity_for_merge <- diversity_outputs$diversity_distances %>%
    dplyr::rename(distance_km_diversity = distance_km)
  abundance_for_merge <- abundance_outputs$abundance_distances %>%
    dplyr::rename(distance_km_abundance = distance_km)

  diversity_abundance_combined <- dplyr::left_join(
    diversity_for_merge %>% dplyr::select(-seqs),
    abundance_for_merge %>%
      dplyr::select(OTU, cluster, total_abundance, log_abundance, distance_km_abundance),
    by = c("OTU", "cluster")
  )

  final_combined <- dplyr::left_join(
    diversity_abundance_combined,
    connection_outputs$connections_distances %>%
      dplyr::select(OTU, cluster, n_connections, log_connections, distance_km_connections),
    by = c("OTU", "cluster")
  )

  metrics_dir <- file.path(otus_root, phylum_name, "informative_otus_metrics")
  network_dir <- file.path(otus_root, phylum_name, "haplotype_network")
  dir.create(metrics_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(network_dir, showWarnings = FALSE, recursive = TRUE)

  write.csv(
    diversity_outputs$diversity_metrics,
    file.path(metrics_dir, "otu_metrics_diversity.csv"),
    row.names = FALSE
  )
  write.csv(
    diversity_outputs$diversity_distances %>% dplyr::select(-seqs),
    file.path(metrics_dir, "diversity_vs_distance.csv"),
    row.names = FALSE
  )
  write.csv(
    abundance_outputs$abundance_metrics,
    file.path(metrics_dir, "otu_metrics_abundance.csv"),
    row.names = FALSE
  )
  write.csv(
    abundance_outputs$abundance_distances,
    file.path(metrics_dir, "abundance_vs_distance.csv"),
    row.names = FALSE
  )
  write.csv(
    connection_outputs$connections_distances,
    file.path(metrics_dir, "otu_metrics_connections.csv"),
    row.names = FALSE
  )
  write.csv(
    final_combined,
    file.path(metrics_dir, "div_abun_conn_combined.csv"),
    row.names = FALSE
  )
  write.csv(
    connection_outputs$network_edges,
    file.path(network_dir, "edges_Mi_filo.csv"),
    row.names = FALSE
  )
  write.csv(
    connection_outputs$network_points,
    file.path(network_dir, "points_Mi_filo.csv"),
    row.names = FALSE
  )

  message("Haplotype network files written to: ", network_dir)

  if (write_maps) {
    build_maps(
      metrics_dir,
      diversity_outputs$diversity_results,
      abundance_outputs$abundance_distances,
      connection_outputs$connections_distances
    )
  }

  message("Completed analysis for: ", phylum_name)
}
