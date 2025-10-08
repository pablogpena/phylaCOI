# LIBRERÍAS NECESARIAS
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
})

# DIRECTORIO RAÍZ
base_dir <- ("C:/TEMP/FASE1/3_magia_nombre_carpeta/eKOI_metabarcoding_database_FILOS")

# LISTA DE CARPETAS DE FILOS
phyla_dirs <- list.dirs(base_dir, recursive = FALSE)

# FUNCIÓN PRINCIPAL PARA PROCESAR CADA FILO
process_phylum <- function(phylum_path) {
  phylum_name <- basename(phylum_path)
  message("\n========================")
  message("Procesando filo: ", phylum_name)
  message("========================")
  
  # Rutas de entrada
  abundances_file     <- file.path(phylum_path, "output", "abundances_unique_actualizado.csv")
  alignment_file      <- file.path(phylum_path, "output", "aligned_sequences_mafft.fasta")
  otu_mapping_file    <- file.path(phylum_path, "output", "otus", "otus_mapping.txt")
  filtered_otus_file  <- file.path(phylum_path, "output", "otus", "informative_OTU.txt")
  
  # Verificar existencia
  required_files <- c(abundances_file, alignment_file, otu_mapping_file, filtered_otus_file)
  if (!all(file.exists(required_files))) {
    warning("Faltan archivos en: ", phylum_name)
    return(NULL)
  }
  
  # Leer datos
  filtered_otus <- read.table(filtered_otus_file, stringsAsFactors = FALSE)
  abundances <- read.csv(abundances_file, stringsAsFactors = FALSE)
  otu_map <- read.table(otu_mapping_file, col.names = c("OTU", "UniqueID"))
  abundances$OTU <- otu_map$OTU[match(abundances$UniqueID, otu_map$UniqueID)]
  alignment <- readDNAStringSet(alignment_file)
  alignment_dnabin <- as.DNAbin(alignment)
  abundances_filtered <- dplyr::filter(abundances, OTU %in% filtered_otus$V1)
  abundances_filtered$Haplotype <- NA
 
  

  #----------------------------------------------------- 
  # ASIGNAR HAPLOTIPOS: descarta otus con >500 secuencias??
  #-----------------------------------------------------
  otus_con_red <- c()
  otus_descartados <- c()
  for (otu in unique(abundances_filtered$OTU)) {
    otu_ids <- abundances_filtered$UniqueID[abundances_filtered$OTU == otu]
    if (length(otu_ids) > 500) {
      otus_descartados <- c(otus_descartados, otu)
      next
    }
    otu_alignment <- alignment_dnabin[names(alignment_dnabin) %in% otu_ids]
    if (length(otu_alignment) > 1) {
      haps <- haplotype(otu_alignment)
      idx <- attr(haps, "index")
      seq_nm <- names(otu_alignment)
      for (h in seq_along(idx)) {
        abundances_filtered$Haplotype[abundances_filtered$UniqueID %in% seq_nm[idx[[h]]]] <- h
      }
      otus_con_red <- c(otus_con_red, otu)
    }
  }
  abundances_filtered <- dplyr::filter(abundances_filtered, OTU %in% otus_con_red)
  
  # Requiere columnas lat/lon. Si no están, salta.
  if (!all(c("lat", "lon") %in% colnames(abundances_filtered))) {
    warning("Faltan columnas lat/lon en: ", phylum_name)
    return(NULL)
  }
  
  #----------------------------------------------------- 
  # AGRUPACIÓN PUNTOS GEOGRÁFICOS: agrupa secuencias en 
  # un radio de 5KM por otu
  #-----------------------------------------------------
  group_points <- function(df, threshold = 5000) {
    coords <- df %>% dplyr::select(lon, lat) %>% as.matrix()
    dist_matrix <- distm(coords, fun = distHaversine)
    clusters <- rep(NA, nrow(coords))
    cluster_id <- 1
    for (i in seq_len(nrow(dist_matrix))) {
      if (is.na(clusters[i])) {
        nearby <- which(dist_matrix[i, ] <= threshold)
        clusters[nearby] <- cluster_id
        cluster_id <- cluster_id + 1
      }
    }
    df$cluster <- paste0(df$OTU[1], "_", clusters)
    return(df)
  }
  
  abundances_grouped <- abundances_filtered %>%
    group_by(OTU) %>%
    group_modify(~ group_points(.x)) %>%
    ungroup()
  
  ### ==========================
  ### ABUNDANCIA
  ### ==========================
  abundance_results <- abundances_grouped %>%
    group_by(OTU, cluster) %>%
    summarise(
      total_abundance = sum(Abundancia, na.rm = TRUE),
      lat = mean(lat),
      lon = mean(lon),
      .groups = "drop"
    ) %>%
    mutate(log_abundance = log1p(total_abundance))
  
  abundance_distances <- abundance_results %>%
    group_by(OTU) %>%
    mutate(
      max_log_abundance = max(log_abundance, na.rm = TRUE),
      center_lat = mean(lat[log_abundance == max_log_abundance], na.rm = TRUE),
      center_lon = mean(lon[log_abundance == max_log_abundance], na.rm = TRUE),
      distance_km = distHaversine(matrix(c(lon, lat), ncol = 2),
                                  matrix(c(center_lon, center_lat), ncol = 2)) / 1000
    ) %>%
    ungroup()
  
  abundance_model_data <- abundance_distances %>%
    group_by(OTU) %>%
    filter(n() >= 3) %>%
    ungroup()
  
  abundance_metrics <- abundance_model_data %>%
    group_by(OTU) %>%
    group_modify(~{
      df <- .x
      lm_model <- lm(log_abundance ~ distance_km, data = df)
      cor_spear <- suppressWarnings(cor(df$log_abundance, df$distance_km, method = "spearman"))
      nls_model <- tryCatch({
        nls(log_abundance ~ a * exp(-b * distance_km),
            data = df, start = list(a = max(df$log_abundance), b = 0.01))
      }, error = function(e) NULL)
      gam_model <- tryCatch({ gam(log_abundance ~ s(distance_km), data = df) }, error = function(e) NULL)
      
      tibble(
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
    ungroup()

  
  ### ==========================
  ### Numero de Conexiones redes de haplotipos
  ### ==========================
  
  conexiones_por_cluster <- list()
  
  for (otu in otus_con_red) {
    cat("\nCalculando conexiones para OTU:", otu, "\n")
    
    otu_data <- dplyr::filter(abundances_grouped, OTU == otu)
    ids <- otu_data$UniqueID
    aln <- alignment_dnabin[names(alignment_dnabin) %in% ids]
    
    if (length(aln) < 2) {
      cat("   -> Menos de 2 secuencias. Se omite.\n")
      next
    }
    
    haps <- pegas::haplotype(aln)
    if (nrow(haps) < 2) {
      cat("   -> Solo un haplotipo. Se omite.\n")
      next
    }
    
    net <- pegas::haploNet(haps)
    labels <- attr(net, "labels")
    haplo_to_seq <- lapply(attr(haps, "index"), function(idx) names(aln)[idx])
    names(haplo_to_seq) <- labels
    
    # Conexiones por diferencias de secuencia
    edges <- apply(as.matrix(net), 1, function(e) {
      from <- haplo_to_seq[[labels[e[1]]]]
      to   <- haplo_to_seq[[labels[e[2]]]]
      base::expand.grid(from = from, to = to, stringsAsFactors = FALSE)
    })
    seq_edges_df <- do.call(rbind, edges)
    
    # Conexiones por mismo haplotipo en distintos clusters
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
      dplyr::ungroup()
    
    # Unir ambas conexiones
    all_edges_df <- dplyr::bind_rows(seq_edges_df, haplo_edges_df) %>% dplyr::distinct()
    
    # Mapear coordenadas
    coord_map <- otu_data %>% dplyr::select(ID = UniqueID, cluster, lat, lon)
    
    get_attr <- function(id, attr) coord_map[[attr]][coord_map$ID == id][1]
    
    all_edges_df <- all_edges_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        from_cluster = get_attr(from, "cluster"),
        to_cluster   = get_attr(to, "cluster"),
        from_lat     = get_attr(from, "lat"),
        from_lon     = get_attr(from, "lon"),
        to_lat       = get_attr(to, "lat"),
        to_lon       = get_attr(to, "lon")
      ) %>%
      dplyr::ungroup()
    
    # Filtrar conexiones entre clusters diferentes
    all_edges_df <- dplyr::filter(all_edges_df, from_cluster != to_cluster)
    
    # Agrupar todos los clusters del OTU
    all_clusters <- otu_data %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(
        lat = mean(lat, na.rm = TRUE),
        lon = mean(lon, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Contar conexiones por cluster (desde "from_cluster")
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
    
    # Completar con clusters sin conexiones
    edge_counts_full <- dplyr::full_join(all_clusters, edge_counts,
                                         by = c("cluster", "lat", "lon")) %>%
      dplyr::mutate(
        n_connections = dplyr::coalesce(n_connections, 0),
        OTU = otu,
        log_connections = log1p(n_connections)
      )
    
    conexiones_por_cluster[[as.character(otu)]] <- edge_counts_full
  }
  
  
  # Unir todos los resultados
  connections_df <- dplyr::bind_rows(conexiones_por_cluster)
  
  
  connections_df2 <- connections_df %>%
    dplyr::group_by(OTU) %>%
    dplyr::mutate(
      max_conn = max(n_connections, na.rm = TRUE)
    ) %>%
    dplyr::ungroup()
  
  # Calcular centro como UNA coordenada (no el promedio)
  
  
  connections_df2 <- connections_df2 %>%
    dplyr::group_by(OTU, cluster) %>%
    dplyr::summarise(
      lat = mean(lat, na.rm = TRUE),
      lon = mean(lon, na.rm = TRUE),
      n_connections = sum(n_connections, na.rm = TRUE),
      log_connections = log1p(n_connections),
      .groups = "drop"
    )
  connections_df2 <- connections_df2 %>%
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
  
  ### ==========================
  ### DIVERSIDAD NUCLEOTÍDICA
  ### ==========================
  diversity_results <- abundances_grouped %>%
    group_by(OTU, cluster) %>%
    summarise(
      seqs = list(UniqueID),
      lat = mean(lat),
      lon = mean(lon),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      diversity = {
        selected <- names(alignment_dnabin) %in% seqs
        aln <- alignment_dnabin[selected]
        if (length(aln) > 1) nuc.div(aln) else 0
      }
    ) %>%
    ungroup()
  
  diversity_results <- diversity_results %>%
    left_join(
      connections_df2 %>% select(OTU, cluster, n_connections),
      by = c("OTU", "cluster")
    )
  
  diversity_distances <- diversity_results %>%
    group_by(OTU) %>%
    mutate(
      max_div = max(diversity, na.rm = TRUE),
      # Identificar índices con diversidad máxima
      center_cluster = {
        idx_max <- which(diversity == max_div)
        
        if (length(idx_max) == 1) {
          # Solo un cluster con diversidad máxima → lo usamos
          idx_max
        } else {
          # --- 1️⃣ Elegir el que tenga más conexiones ---
          conn_values <- n_connections[idx_max]  # columna con número de conexiones
          idx_best <- idx_max[which.max(conn_values)]
          
          # --- 2️⃣ Si hay empate en conexiones, elegir el más central geográficamente ---
          if (sum(conn_values == max(conn_values)) > 1) {
            idx_candidates <- idx_max[conn_values == max(conn_values)]
            
            mean_distances <- sapply(idx_candidates, function(i) {
              mean(
                distHaversine(
                  matrix(c(lon[i], lat[i]), ncol = 2),
                  matrix(c(lon, lat), ncol = 2)
                )
              )
            })
            idx_best <- idx_candidates[which.min(mean_distances)]
          }
          
          idx_best
        }
      },
      center_lat = lat[center_cluster],
      center_lon = lon[center_cluster],
      distance_km = distHaversine(
        matrix(c(lon, lat), ncol = 2),
        matrix(c(center_lon, center_lat), ncol = 2)
      ) / 1000
    ) %>%
    ungroup()
  
  diversity_model_data <- diversity_distances %>%
    group_by(OTU) %>%
    filter(n() >= 3) %>%
    ungroup()
  
  otu_metrics_div <- diversity_model_data %>%
    group_by(OTU) %>%
    group_modify(~{
      df <- .x
      lm_model <- lm(diversity ~ distance_km, data = df)
      lm_slope <- coef(lm_model)[2]
      cor_spear <- suppressWarnings(cor(df$diversity, df$distance_km, method = "spearman"))
      nls_model <- tryCatch({
        nls(diversity ~ a * exp(-b * distance_km),
            data = df, start = list(a = max(df$diversity), b = 0.01))
      }, error = function(e) NULL)
      gam_model <- tryCatch({ gam(diversity ~ s(distance_km), data = df) }, error = function(e) NULL)
      
      tibble(
        max_diversity = max(df$diversity, na.rm = TRUE),
        mean_diversity = mean(df$diversity, na.rm = TRUE),
        max_distance = max(df$distance_km, na.rm = TRUE),
        lm_slope = lm_slope,
        spearman_rho = cor_spear,
        nls_a = if (!is.null(nls_model)) coef(nls_model)["a"] else NA,
        nls_b = if (!is.null(nls_model)) coef(nls_model)["b"] else NA,
        gam_r2 = if (!is.null(gam_model)) summary(gam_model)$r.sq else NA
      )
    }) %>%
    ungroup()
  
  write.csv(otu_metrics_div, "otu_metrics_diversity.csv", row.names = FALSE)
  write.csv(diversity_distances %>% select(-seqs), "diversity_vs_distance.csv", row.names = FALSE)
  
  
  ### =====================================
  ### ARCHIVO COMBINADO: DIV, ABUN, CONEX
  ### =====================================
  
  ### Renombrar columnas de distancia para claridad (diversidad)
  diversity_distances <- diversity_results %>%
    group_by(OTU) %>%
    mutate(
      max_div = max(diversity, na.rm = TRUE),
      
      # Seleccionar el cluster central para cada OTU
      center_cluster = {
        idx_max <- which(diversity == max_div)
        
        if (length(idx_max) == 1) {
          idx_max
        } else {
          # Elegir el que tenga más conexiones
          conn_values <- n_connections[idx_max]
          idx_best <- idx_max[which.max(conn_values)]
          
          # Si hay empate, elegir el más central geográficamente
          if (sum(conn_values == max(conn_values)) > 1) {
            idx_candidates <- idx_max[conn_values == max(conn_values)]
            
            mean_distances <- sapply(idx_candidates, function(i) {
              # Distancia promedio al resto de clusters
              valid <- !is.na(lon) & !is.na(lat)
              mean(distHaversine(
                cbind(lon[valid], lat[valid]),
                cbind(lon[i], lat[i])
              ))
            })
            idx_best <- idx_candidates[which.min(mean_distances)]
          }
          idx_best
        }
      },
      
      center_lat = lat[center_cluster],
      center_lon = lon[center_cluster],
      
      # Distancia de cada cluster al cluster central (en km)
      distance_km_diversity = {
        valid <- !is.na(lat) & !is.na(lon)
        distHaversine(
          cbind(lon[valid], lat[valid]),
          cbind(center_lon, center_lat)
        ) / 1000
      }
    ) %>%
    ungroup()
  
  # Renombrar para abundancia
  abundance_distances <- abundance_results %>%
    group_by(OTU) %>%
    mutate(
      max_log_abundance = max(log_abundance, na.rm = TRUE),
      center_lat = mean(lat[log_abundance == max_log_abundance], na.rm = TRUE),
      center_lon = mean(lon[log_abundance == max_log_abundance], na.rm = TRUE),
      distance_km_abundance = distHaversine(matrix(c(lon, lat), ncol = 2),
                                            matrix(c(center_lon, center_lat), ncol = 2)) / 1000
    ) %>%
    ungroup()
  
  # Combinar diversidad y abundancia
  diversity_abundance_combined <- left_join(
    diversity_distances %>% select(-seqs),  # si existe
    abundance_distances %>%
      select(OTU, cluster, total_abundance, log_abundance, distance_km_abundance),
    by = c("OTU", "cluster")
  )
  
  ##########
  # 🔗 Añadir datos de conexiones de red (haplotipos)
  ##########
  final_combined <- left_join(
    diversity_abundance_combined,
    connections_df2 %>% 
      select(OTU, cluster, n_connections, log_connections, distance_km_connections),
    by = c("OTU", "cluster")
  )
  
  
  ### ==========================
  ### MAPAS INTERACTIVOS
  ### ==========================
  
  # MAPA DIVERSIDAD
  pal_div <- colorNumeric("viridis", domain = diversity_results$diversity)
  map_div <- leaflet() %>%
    addProviderTiles(providers$CartoDB.Positron)
  
  for (otu in unique(diversity_results$OTU)) {
    map_div <- map_div %>%
      addCircleMarkers(data = filter(diversity_results, OTU == otu),
                       lng = ~lon, lat = ~lat,
                       radius = ~rescale(diversity, to = c(4, 12)),
                       color = ~pal_div(diversity),
                       stroke = TRUE, fillOpacity = 0.8,
                       label = ~paste0("OTU: ", OTU,
                                       "<br>Cluster: ", cluster,
                                       "<br>Diversity: ", round(diversity, 5)),
                       group = paste0("OTU ", otu))
  }
  map_div <- map_div %>%
    addLegend("bottomright", pal = pal_div, values = diversity_results$diversity,
              title = "Nucleotide diversity") %>%
    addLayersControl(overlayGroups = paste0("OTU ", unique(diversity_results$OTU)),
                     options = layersControlOptions(collapsed = FALSE))
  
  saveWidget(map_div, file.path(phylum_path, "output", "otu_nucleotide_diversity_map.html"), selfcontained = TRUE)
  
  # MAPA ABUNDANCIA
  pal_abun <- colorNumeric("viridis", domain = abundance_distances$log_abundance)
  map_abun <- leaflet() %>%
    addProviderTiles(providers$CartoDB.Positron)
  
  for (otu in unique(abundance_distances$OTU)) {
    map_abun <- map_abun %>%
      addCircleMarkers(data = filter(abundance_distances, OTU == otu),
                       lng = ~lon, lat = ~lat,
                       radius = ~rescale(log_abundance, to = c(4, 12)),
                       color = ~pal_abun(log_abundance),
                       stroke = TRUE, fillOpacity = 0.8,
                       label = ~paste0("OTU: ", OTU,
                                       "<br>Cluster: ", cluster,
                                       "<br>Log(Abundancia): ", round(log_abundance, 3)),
                       group = paste0("OTU ", otu))
  }
  map_abun <- map_abun %>%
    addLegend("bottomright", pal = pal_abun, values = abundance_distances$log_abundance,
              title = "Log(Abundance)") %>%
    addLayersControl(overlayGroups = paste0("OTU ", unique(abundance_distances$OTU)),
                     options = layersControlOptions(collapsed = FALSE))
  
  saveWidget(map_abun, file.path(phylum_path, "output", "otu_log_abundance_map.html"), selfcontained = TRUE)
  
  ############conexiones
  
  # Paleta de colores para log(conexiones)
  pal_conn <- colorNumeric("viridis", domain = connections_df2$log_connections, na.color = "transparent")
  
  map_conn <- leaflet() %>%
    addProviderTiles(providers$CartoDB.Positron)
  
  # Añadir capas por OTU
  for (otu in unique(connections_df2$OTU)) {
    map_conn <- map_conn %>%
      addCircleMarkers(data = dplyr::filter(connections_df2, OTU == otu),
                       lng = ~lon, lat = ~lat,
                       radius = ~rescale(log_connections, to = c(4, 12)),
                       color = ~pal_conn(log_connections),
                       stroke = TRUE, fillOpacity = 0.8,
                       label = ~paste0("OTU: ", OTU,
                                       "<br>Cluster: ", cluster,
                                       "<br>Log(Conexiones): ", round(log_connections, 3)),
                       group = paste0("OTU ", otu))
  }
  
  # Añadir leyenda y control de capas
  map_conn <- map_conn %>%
    addLegend("bottomright", pal = pal_conn, values = connections_df2$log_connections,
              title = "Log(Connections)") %>%
    addLayersControl(overlayGroups = paste0("OTU ", unique(connections_df2$OTU)),
                     options = layersControlOptions(collapsed = FALSE))
  
  # Guardar mapa
  htmlwidgets::saveWidget(map_conn, file.path(phylum_path, "output", "otu_log_connections_map.html"), selfcontained = TRUE)
  
  
  
  
  
  ### =============================
  ### GUARDAR TODOS LOS RESULTADOS
  ### ============================= 
  output_dir <- file.path(phylum_path, "output")
  
  #DIVERSIDAD NUCLEOTÍDICA
  write.csv(otu_metrics_div, file.path(phylum_path, "output", "otu_metrics_diversity.csv"), row.names = FALSE)
  write.csv(diversity_distances %>% select(-seqs), file.path(phylum_path, "output", "diversity_vs_distance.csv"), row.names = FALSE)
  #ABUNDANCIA
  write.csv(abundance_metrics, file.path(phylum_path, "output", "otu_metrics_abundance.csv"), row.names = FALSE)
  write.csv(abundance_distances, file.path(phylum_path, "output", "abundance_vs_distance.csv"), row.names = FALSE)
  #CONEXIONES HAPLOTIPOS
  write.csv(connections_df2, file.path(phylum_path, "output", "otu_metrics_connections.csv"), row.names = FALSE)
  #COMBINADO DIVERSIDAD, ABUNDANCIA Y CONEXIONES
  write.csv(final_combined, file.path(phylum_path, "output", "div_abun_conn_combined.csv"), row.names = FALSE)
  
  
  message("Análisis completado para: ", phylum_name)
}

# APLICAR A TODOS LOS FILOS
for (phylum_path in phyla_dirs) {
  tryCatch({
    process_phylum(phylum_path)
  }, error = function(e) {
    message("Error en ", basename(phylum_path), ": ", e$message)
  })
}
