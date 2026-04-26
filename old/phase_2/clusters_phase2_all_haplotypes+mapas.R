#############################################
### Part 1: ALL HAPLOTYPES — 
### Mantiene 4 métodos: inverse, count, linear, RBF
#############################################

# Reproducibilidad global
set.seed(42)

# ============================
# Librerías
# ============================
library(dplyr)
library(igraph)
library(tidyr)
library(spdep)
library(aricode)
library(ggplot2)
library(mclust)
library(stringr)
library(stringi)
library(geosphere)   # distHaversine (m)
library(rstatix)

# ============================
# Rutas y carga de datos
# ============================
setwd("E:/TEMP/prueba_script_definitivos/_RESULTADOS_HAPLOTIPOS_SCRIPT_7/A_all_haplotypes/1_METAZOOS")

root_dir <- "E:/TEMP/prueba_script_definitivos/_RESULTADOS_HAPLOTIPOS_SCRIPT_7/A_all_haplotypes/1_METAZOOS"
filos <- list.dirs(root_dir, recursive = FALSE)

all_edges <- list()
all_points <- list()

for (filo_dir in filos) {
  filo_name <- basename(filo_dir)
  edges_file <- file.path(filo_dir, "edges_Mi_filo.csv")
  points_file <- file.path(filo_dir, "points_Mi_filo.csv")
  
  if (file.exists(edges_file) & file.exists(points_file)) {
    message("✅ Cargando datos de: ", filo_name)
    edges_df <- read.csv(edges_file)
    points_df <- read.csv(points_file)
    edges_df$Filo <- filo_name
    points_df$Filo <- filo_name
    all_edges[[filo_name]] <- edges_df
    all_points[[filo_name]] <- points_df
  } else {
    message("⚠️ No se encontraron archivos para ", filo_name)
  }
}

edges_df_all <- bind_rows(all_edges)
points_df_all <- bind_rows(all_points)

# ============================
# 1) Construcción de edges_with_info (filtrado por OTU y quitando intra-localidad)
# ============================
edges_with_info <- edges_df_all %>%
  left_join(points_df_all %>% select(UniqueID, Locality_from = Localities, OTU_ID_from = OTU_ID, Filo_from = Filo),
            by = c("from" = "UniqueID")) %>%
  left_join(points_df_all %>% select(UniqueID, Locality_to = Localities, OTU_ID_to = OTU_ID, Filo_to = Filo),
            by = c("to" = "UniqueID")) %>%
  filter(Locality_from != Locality_to) %>%
  filter(OTU_ID_from == OTU_ID_to) %>%
  mutate(OTU_ID = OTU_ID_from)

# ============================
# 2) Partición A: clusters ponderados por distancia genética (inversa)
# ============================
# --- A) Inversa (1 / distancia genética) ---
edges_locality <- edges_with_info %>%
  group_by(Locality_from, Locality_to) %>%
  summarise(
    mean_dist = mean(distancia_genetica, na.rm = TRUE),
    n_connections = n(),
    .groups = "drop"
  )

g_locality <- graph_from_data_frame(edges_locality, directed = FALSE)
E(g_locality)$weight <- 1 / (edges_locality$mean_dist + 1e-6)
set.seed(42)
clusters <- cluster_louvain(g_locality, weights = E(g_locality)$weight)

locality_clusters <- data.frame(
  Locality = names(membership(clusters)),
  Cluster  = as.integer(membership(clusters))
)

locality_coords <- edges_with_info %>%
  group_by(Locality_from) %>%
  summarise(
    lon = mean(x, na.rm = TRUE),
    lat = mean(y, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(Locality = Locality_from) %>%
  left_join(locality_clusters, by = "Locality")

# ============================
# 3) Partición B: clusters ponderados por Nº de conexiones (sin distancia)
# ============================
# --- B) Conteo (número de conexiones) ---
edges_locality_counts <- edges_with_info %>%
  group_by(Locality_from, Locality_to) %>%
  summarise(n_connections = n(), .groups = "drop") %>%
  filter(Locality_from != Locality_to)

g_locality_cnt <- graph_from_data_frame(edges_locality_counts, directed = FALSE)
E(g_locality_cnt)$weight <- edges_locality_counts$n_connections
set.seed(42)
cl_cnt <- cluster_louvain(g_locality_cnt, weights = E(g_locality_cnt)$weight)

locality_clusters_cnt <- data.frame(
  Locality = names(membership(cl_cnt)),
  Cluster_cnt = as.integer(membership(cl_cnt))
)

locality_coords_cnt <- edges_with_info %>%
  group_by(Locality_from) %>%
  summarise(
    lon = mean(x, na.rm = TRUE),
    lat = mean(y, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(Locality = Locality_from) %>%
  left_join(locality_clusters_cnt, by = "Locality")

# ============================
# 4) Partición C: clusters ponderados linealmente (1 - d/max)
# ============================

# Función para grafo con pesos
make_graph_with_weights <- function(edges_locality, mode = c("inverse","linear","rbf"), sigma = NULL) {
  mode <- match.arg(mode)
  g <- graph_from_data_frame(d = edges_locality, directed = FALSE)
  d <- edges_locality$mean_dist
  
  if (mode == "inverse") {
    w <- 1 / (d + 1e-6)
  } else if (mode == "linear") {
    maxd <- max(d, na.rm = TRUE)
    if (maxd == 0) maxd <- 1
    w <- 1 - d / maxd
  } else if (mode == "rbf") {
    if (is.null(sigma)) stop("Para 'rbf' debes proporcionar sigma")
    w <- exp(- d / sigma)
  }
  E(g)$weight <- w
  E(g)$dist <- d
  return(g)
}

# Función auxiliar: clustering y coordenadas
cluster_and_coords <- function(g, edges_with_info, label_name) {
  set.seed(42)
  cl <- cluster_louvain(g, weights = E(g)$weight)
  membership_df <- data.frame(
    Locality = names(membership(cl)),
    Cluster  = as.integer(membership(cl))
  )
  coords <- edges_with_info %>%
    group_by(Locality_from) %>%
    summarise(
      lon = mean(x, na.rm = TRUE),
      lat = mean(y, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(Locality = Locality_from) %>%
    left_join(membership_df, by = "Locality")
  names(coords)[names(coords)=="Cluster"] <- label_name
  list(cl = cl, coords = coords)
}

# --- Partición C: LINEAL ---
g_locality_linear <- make_graph_with_weights(edges_locality, mode = "linear")
res_linear <- cluster_and_coords(g_locality_linear, edges_with_info, "Cluster_linear")
cl_linear  <- res_linear$cl
locality_coords_linear <- res_linear$coords

# ============================
# 5) Partición D: clusters ponderados por función RBF (exp(-d/sigma))
# ============================

# Selección automática del mejor sigma
sigma_grid <- c(0.005, 0.010, 0.015)
rbf_grid <- lapply(sigma_grid, function(s) {
  gtmp <- make_graph_with_weights(edges_locality, mode = "rbf", sigma = s)
  cltmp <- cluster_louvain(gtmp, weights = E(gtmp)$weight)
  data.frame(sigma = s,
             modularidad = modularity(gtmp, membership(cltmp), weights = E(gtmp)$weight),
             n_comun = length(unique(membership(cltmp))))
})
rbf_grid <- do.call(rbind, rbf_grid)
print(rbf_grid)

# Escoger el sigma que maximiza la modularidad
best_sigma <- rbf_grid$sigma[which.max(rbf_grid$modularidad)]
cat("\nSigma RBF seleccionado (máxima modularidad):", best_sigma, "\n")

# Construir grafo y clustering final
g_locality_rbf <- make_graph_with_weights(edges_locality, mode = "rbf", sigma = best_sigma)
set.seed(42)
res_rbf <- cluster_and_coords(g_locality_rbf, edges_with_info, "Cluster_rbf")
cl_rbf  <- res_rbf$cl
locality_coords_rbf <- res_rbf$coords

# ======================================================
# 6) Comparaciones entre las cuatro particiones
#    (Inverse, Conexiones, Lineal, RBF)
# ======================================================

# --- Alinear todas las particiones por localidad ---
comp_all <- locality_coords %>%
  select(Locality, Cluster_inv = Cluster) %>%
  inner_join(locality_coords_cnt %>% select(Locality, Cluster_cnt), by = "Locality") %>%
  inner_join(locality_coords_linear %>% select(Locality, Cluster_linear), by = "Locality") %>%
  inner_join(locality_coords_rbf %>% select(Locality, Cluster_rbf), by = "Locality")

write.csv(comp_all, "clusters_metazoa.csv", row.names = FALSE)

# --- ARI entre pares de particiones ---
ari_table <- data.frame(
  Comparacion = c("Inv vs Conex", "Inv vs Lineal", "Inv vs RBF",
                  "Conex vs Lineal", "Conex vs RBF", "Lineal vs RBF"),
  ARI = c(
    adjustedRandIndex(comp_all$Cluster_inv,    comp_all$Cluster_cnt),
    adjustedRandIndex(comp_all$Cluster_inv,    comp_all$Cluster_linear),
    adjustedRandIndex(comp_all$Cluster_inv,    comp_all$Cluster_rbf),
    adjustedRandIndex(comp_all$Cluster_cnt,    comp_all$Cluster_linear),
    adjustedRandIndex(comp_all$Cluster_cnt,    comp_all$Cluster_rbf),
    adjustedRandIndex(comp_all$Cluster_linear, comp_all$Cluster_rbf)
  )
)

cat("\n=== Índice Rand Ajustado (ARI) entre particiones ===\n")
print(ari_table, row.names = FALSE)
write.csv(ari_table, "ari_metazoa.csv", row.names = FALSE)

#----------------------------------------------
# --- Modularidad y número de comunidades ---
#----------------------------------------------
mod_summary <- dplyr::bind_rows(
  data.frame(metodo = "Inverse_1_d",
             modularidad = modularity(g_locality, membership(clusters), weights = E(g_locality)$weight),
             n_comun = length(unique(membership(clusters)))),
  data.frame(metodo = "Conexiones",
             modularidad = modularity(g_locality_cnt, membership(cl_cnt), weights = E(g_locality_cnt)$weight),
             n_comun = length(unique(membership(cl_cnt)))),
  data.frame(metodo = "Lineal_1_minus_d_over_max",
             modularidad = modularity(g_locality_linear, membership(cl_linear), weights = E(g_locality_linear)$weight),
             n_comun = length(unique(membership(cl_linear)))),
  data.frame(metodo = paste0("RBF_sigma_", best_sigma),
             modularidad = modularity(g_locality_rbf, membership(cl_rbf), weights = E(g_locality_rbf)$weight),
             n_comun = length(unique(membership(cl_rbf))))
)

cat("\n=== Resumen modularidad y nº comunidades ===\n")
print(mod_summary[order(-mod_summary$modularidad), ], row.names = FALSE)
write.csv(mod_summary, "modularity_summary_metazoa.csv", row.names = FALSE)


# ======================================================
# 7) Autocorrelación espacial (Moran’s I por cluster)
# ======================================================

# Función general para calcular Moran por clúster
moran_table_by_cluster <- function(locality_coords_df, cluster_col,
                                   k = 5,
                                   jitter_amount = 1e-06,
                                   decimals_fixed = 8,
                                   tiny_thresh = 1e-12,
                                   round_I = 3) {
  stopifnot(all(c("Locality","lon","lat", cluster_col) %in% names(locality_coords_df)))
  
  coords_all <- locality_coords_df %>%
    select(Locality, lon, lat, !!rlang::sym(cluster_col)) %>%
    tidyr::drop_na(lon, lat)
  
  xy <- as.matrix(coords_all[, c("lon","lat")])
  if (jitter_amount > 0) {
    set.seed(1)
    xy <- jitter(xy, amount = jitter_amount)
  }
  
  nb <- knn2nb(knearneigh(xy, k = k))
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  labs <- factor(coords_all[[cluster_col]])
  X <- model.matrix(~ labs - 1)
  colnames(X) <- levels(labs)
  
  mor_mat <- t(sapply(1:ncol(X), function(j){
    mt <- moran.test(as.numeric(X[, j]), lw, zero.policy = TRUE)
    c(I = unname(mt$estimate["Moran I statistic"]), p = mt$p.value)
  }))
  
  mor_df <- as.data.frame(mor_mat)
  mor_df$Cluster <- colnames(X)
  
  freq_df <- data.frame(Cluster = as.character(labs)) %>%
    count(Cluster, name = "Freq")
  
  out <- mor_df %>%
    select(Cluster, I, p) %>%
    left_join(freq_df, by = "Cluster") %>%
    mutate(q_FDR = p.adjust(p, method = "BH")) %>%
    arrange(desc(I))
  
  format_small <- function(x) {
    ifelse(x < tiny_thresh,
           paste0("<", formatC(tiny_thresh, format = "e", digits = 0)),
           sprintf(paste0("%.", decimals_fixed, "f"), x))
  }
  
  out_formatted <- out %>%
    mutate(
      I     = round(I, round_I),
      p     = format_small(p),
      q_FDR = format_small(q_FDR)
    )
  
  list(raw = out, formatted = out_formatted)
}

# Calcular Moran para cada método
res_inv  <- moran_table_by_cluster(locality_coords,        "Cluster")
res_cnt  <- moran_table_by_cluster(locality_coords_cnt,    "Cluster_cnt")
res_lin  <- moran_table_by_cluster(locality_coords_linear, "Cluster_linear")
res_rbf  <- moran_table_by_cluster(locality_coords_rbf,    "Cluster_rbf")

# Mostrar resultados formateados
cat("\n=== Moran por cluster (Inverse 1/d) ===\n"); print(res_inv$formatted, row.names = FALSE)
cat("\n=== Moran por cluster (Conexiones) ===\n");  print(res_cnt$formatted, row.names = FALSE)
cat("\n=== Moran por cluster (Lineal) ===\n");      print(res_lin$formatted, row.names = FALSE)
cat("\n=== Moran por cluster (RBF) ===\n");         print(res_rbf$formatted, row.names = FALSE)

# Guardar tablas (opcional)
write.csv(res_inv$formatted, "moran_inverse_metazoa.csv", row.names = FALSE)
write.csv(res_cnt$formatted, "moran_count_metazoa.csv", row.names = FALSE)
write.csv(res_lin$formatted, "moran_linear_metazoa.csv", row.names = FALSE)
write.csv(res_rbf$formatted, "moran_rbf_metazoa.csv", row.names = FALSE)


# --- Resumen final de modularidad (para visualización rápida) ---
mod_results <- data.frame(
  Metodo = c("Inverse_1_d", "Conexiones", "Lineal", "RBF"),
  Modularity = c(
    modularity(g_locality, membership(clusters), weights = E(g_locality)$weight),
    modularity(g_locality_cnt, membership(cl_cnt), weights = E(g_locality_cnt)$weight),
    modularity(g_locality_linear, membership(cl_linear), weights = E(g_locality_linear)$weight),
    modularity(g_locality_rbf, membership(cl_rbf), weights = E(g_locality_rbf)$weight)
  )
)

cat("\n=== Resumen final de modularidad ===\n")
print(mod_results)
write.csv(mod_results, "modularity_results_metazoa.csv", row.names = FALSE)

# ======================================================
# 8) Porcentaje de localidades por océano (Mediterráneo / Atlántico)
# ======================================================

# Lectura y limpieza de ocean_df

ocean_csv <- "metadata_eKOI_ver2.csv"
ocean_df <- read.csv(ocean_csv, sep = ";", stringsAsFactors = FALSE) %>%
  transmute(
    Locality = id_sample,
    Oceans   = trimws(gsub("\\s+", " ", Oceans))
  )


# Normalización robusta de claves (quita espacios raros/acentos y pone MAYÚSCULAS)

norm_key <- function(x) {
  x %>%
    str_replace_all("\\s+", " ") %>%                       # colapsa espacios internos
    str_trim() %>%                                         # quita espacios extremos
    str_replace_all("\u00A0", " ") %>%                     # NBSP -> espacio normal
    stri_trans_general("NFD; [:Nonspacing Mark:] Remove; NFC") %>%  # sin diacríticos
    toupper()
}

# Normaliza claves en ocean_df

ocean_df <- ocean_df %>%
  mutate(Locality = norm_key(Locality),
         Oceans   = str_trim(Oceans)) %>%
  distinct(Locality, .keep_all = TRUE)   # 1 fila por Locality


# Estandariza tablas de caso (renombra y normaliza claves)

locality_coords <- locality_coords %>%
  rename(Locality = Locality) %>%
  mutate(
    Locality = norm_key(Locality),
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  )

locality_coords_cnt <- locality_coords_cnt %>%
  mutate(
    Locality = norm_key(Locality),
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  )

locality_coords_linear <- locality_coords_linear %>%
  mutate(
    Locality = norm_key(Locality),
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  )

locality_coords_rbf <- locality_coords_rbf %>%
  mutate(
    Locality = norm_key(Locality),
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  )


# Función: añade % Mediterráneo y % Atlántico Norte a tabla de Moran

# --- Helper: calcular % de océanos por cluster ---
ocean_share_by_cluster <- function(coords_df, cluster_col, ocean_df) {
  stopifnot(all(c("Locality", cluster_col) %in% names(coords_df)))
  stopifnot(all(c("Locality", "Oceans") %in% names(ocean_df)))
  
  missing_in_ocean <- coords_df %>%
    anti_join(ocean_df, by = "Locality") %>%
    distinct(Locality)
  if (nrow(missing_in_ocean) > 0) {
    message(sprintf("[%s] Aviso: %d localidades sin 'Oceans'.", cluster_col, nrow(missing_in_ocean)))
  }
  
  out <- coords_df %>%
    select(Locality, Cluster = all_of(cluster_col)) %>%
    inner_join(ocean_df, by = "Locality") %>%
    mutate(Cluster = as.integer(Cluster)) %>%
    group_by(Cluster, Oceans) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup() %>%
    select(Cluster, Oceans, pct) %>%
    pivot_wider(names_from = Oceans, values_from = pct, values_fill = 0) %>%
    mutate(across(-Cluster, ~round(.x, 1))) %>%
    arrange(Cluster)
  return(out)
}

# --- Ejecutar para cada método ---
tabla_coords        <- ocean_share_by_cluster(locality_coords,        "Cluster",        ocean_df)
tabla_coords_cnt    <- ocean_share_by_cluster(locality_coords_cnt,    "Cluster_cnt",    ocean_df)
tabla_coords_linear <- ocean_share_by_cluster(locality_coords_linear, "Cluster_linear", ocean_df)
tabla_coords_rbf    <- ocean_share_by_cluster(locality_coords_rbf,    "Cluster_rbf",    ocean_df)

# --- Guardar resultados ---
write.csv(tabla_coords,        "tabla_coords_inverse_metazoa.csv", row.names = FALSE)
write.csv(tabla_coords_cnt,    "tabla_coords_count_metazoa.csv",   row.names = FALSE)
write.csv(tabla_coords_linear, "tabla_coords_linear_metazoa.csv",  row.names = FALSE)
write.csv(tabla_coords_rbf,    "tabla_coords_rbf_metazoa.csv",     row.names = FALSE)

# ============================
# 9) Zonas oceanográficas por corrientes (ALL haplotypes)
# ============================

zones_df <- locality_coords %>%
  mutate(
    Zone = case_when(
      lon < 0  & lat >= 37 ~ "ATL_N_UPW",
      lon < -2 & lat < 37  ~ "ATL_S_UPW_ALM",
      lon >= -2            ~ "MED",
      TRUE                 ~ NA_character_
    )
  ) %>%
  select(Locality, Zone)

# Control rápido
cat("\n=== Control de asignación de zonas (ALL haplotypes) ===\n")
print(table(is.na(zones_df$Zone)))

write.csv(zones_df, "locality_zones_currents_all.csv", row.names = FALSE)

#Definimos la función:
zone_share_by_cluster <- function(coords_df, cluster_col, zones_df) {
  stopifnot(all(c("Locality", cluster_col) %in% names(coords_df)))
  
  coords_df %>%
    select(Locality, Cluster = all_of(cluster_col)) %>%
    inner_join(zones_df, by = "Locality") %>%
    count(Cluster, Zone) %>%
    group_by(Cluster) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup() %>%
    select(-n) %>%
    pivot_wider(
      names_from = Zone,
      values_from = pct,
      values_fill = 0
    ) %>%
    mutate(across(where(is.numeric), ~round(.x, 1))) %>%
    arrange(Cluster)
}


#Exportamos resultados:
zones_inv <- zone_share_by_cluster(
  locality_coords,
  "Cluster",
  zones_df
)

zones_cnt <- zone_share_by_cluster(
  locality_coords_cnt,
  "Cluster_cnt",
  zones_df
)

zones_linear <- zone_share_by_cluster(
  locality_coords_linear,
  "Cluster_linear",
  zones_df
)

zones_rbf <- zone_share_by_cluster(
  locality_coords_rbf,
  "Cluster_rbf",
  zones_df
)

write.csv(zones_inv,    "zone_share_inverse_metazoa.csv", row.names = FALSE)
write.csv(zones_cnt,    "zone_share_count_metazoa.csv",   row.names = FALSE)
write.csv(zones_linear, "zone_share_linear_metazoa.csv",  row.names = FALSE)
write.csv(zones_rbf,    "zone_share_rbf_metazoa.csv",     row.names = FALSE)

#-------------------------------------
# 10) Maps and results visualization
#-------------------------------------
#     - Paleta fija (16) + extensión automática si hay >16 clusters
#     - Colores consistentes entre los 4 mapas (mismo número = mismo color)
#     - Añade Francia (blanco) + Andorra para que Iberia no parezca isla
#     - Recorte norte controlado con ymax (p.ej. 44.5)
#     - Fix S2 (sf) para intersecciones/union en lon/lat
#     - Compatibilidad ggplot2: usa size= en vez de linewidth=
#     - Clusters SIN transparencia: voronoi_alpha = 1 (opaco)
#-------------------------------------

# ============================
# Librerías necesarias
# ============================
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(scatterpie)
library(rnaturalearth)
library(RColorBrewer)
library(ggspatial)
library(stringr)

# ============================
# 0) Paleta manual (16 colores) - EDITA AQUÍ si quieres
# ============================
pal_user16 <- c(
  "1"  = "#8DD3C7",
  "2"  = "#FDB462",
  "3"  = "#B3DE69",
  "4"  = "#FCCDE5",
  "5"  = "#D9D9D9",
  "6"  = "#BC80BD", 
  "7"  = "#CCEBC5",
  "8"  = "#FFED6F",
  "9"  = "#7F7F7F",
  "10" = "#FFFFB3",
  "11" = "#BEBADA",
  "12" = "#fb8072",
  "13" = "#80B1D3", 
  "14" = "#C49C94",
  "15" = "#AEC7E8",
  "16" = "#304847"
)

# ============================
# 1) Funciones base
# ============================

# 1.1) Construye paleta global: usa tus 16 y añade colores si faltan clusters (>16)
make_pal_from_user <- function(comp_all, pal_user, hcl_palette = "Dark 3") {
  levels_all <- sort(unique(unlist(lapply(
    comp_all[, c("Cluster_inv","Cluster_cnt","Cluster_linear","Cluster_rbf")],
    function(x) as.character(unique(x))
  ))))
  
  missing <- setdiff(levels_all, names(pal_user))
  if (length(missing) > 0) {
    extra_cols <- grDevices::hcl.colors(length(missing), palette = hcl_palette)
    names(extra_cols) <- missing
    pal_full <- c(pal_user, extra_cols)
  } else {
    pal_full <- pal_user
  }
  
  pal_full[levels_all]
}

# 1.2) Capas base (Iberia + Francia parcial + Andorra + bounding box)
#      - ymax controla cuánto de Francia entra
#      - FIX: desactiva S2 temporalmente para evitar errores de modelo (s2)
prepare_base_layers <- function(xmin = -10, ymin = 35.5, xmax = 4, ymax = 44.5, crs_planar = 3035) {
  
  old_s2 <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  
  # Iberia
  iberia_wgs <- world %>%
    dplyr::filter(admin %in% c("Spain", "Portugal")) %>%
    sf::st_union() %>% sf::st_as_sf()
  
  # Francia
  france_wgs <- world %>%
    dplyr::filter(admin == "France") %>%
    sf::st_union() %>% sf::st_as_sf()
  
  # Andorra (NO unir: mantener geometría original; más estable)
  andorra_wgs <- world %>% dplyr::filter(admin == "Andorra")
  if (nrow(andorra_wgs) == 0) {
    andorra_wgs <- world %>% dplyr::filter(sovereignt == "Andorra")
  }
  andorra_wgs <- sf::st_make_valid(andorra_wgs)
  
  # BBox en WGS84
  bbox_wgs <- sf::st_as_sfc(sf::st_bbox(c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), crs = 4326))
  
  # Recortes
  iberia_wgs_clip  <- suppressWarnings(sf::st_intersection(iberia_wgs,  bbox_wgs))
  france_wgs_clip  <- suppressWarnings(sf::st_intersection(france_wgs,  bbox_wgs))
  andorra_wgs_clip <- suppressWarnings(sf::st_intersection(andorra_wgs, bbox_wgs))
  
  # Land mask: Iberia + Francia + Andorra
  land_wgs_clip <- suppressWarnings(sf::st_union(iberia_wgs_clip, france_wgs_clip, andorra_wgs_clip)) %>% sf::st_as_sf()
  
  list(
    iberia_wgs_clip  = iberia_wgs_clip,
    france_wgs_clip  = france_wgs_clip,
    andorra_wgs_clip = andorra_wgs_clip,
    land_wgs_clip    = land_wgs_clip,
    bbox_wgs         = bbox_wgs,
    crs_planar       = crs_planar
  )
}

# 1.3) Generador de mapas
make_map_for_case <- function(comp_all, cluster_col, title_text,
                              base_layers = NULL,
                              cellsize = 60000,     # ~60 km
                              pal = NULL,
                              voronoi_alpha = 0.25,    
                              hcl_palette = "Dark 3") {
  
  stopifnot(all(c("Locality", "lon", "lat", cluster_col) %in% names(comp_all)))
  
  df_case <- comp_all %>%
    select(Locality, lon, lat, Cluster = all_of(cluster_col)) %>%
    drop_na(lon, lat)
  
  # Si no viene pal, crea una
  if (is.null(pal)) {
    levs <- sort(unique(as.character(df_case$Cluster)))
    nlev <- length(levs)
    if (nlev <= 12) {
      pal <- RColorBrewer::brewer.pal(nlev, "Set3")
      names(pal) <- levs
    } else {
      pal <- grDevices::hcl.colors(nlev, palette = hcl_palette)
      names(pal) <- levs
    }
  }
  
  # Mapping estable
  df_case$Cluster <- factor(as.character(df_case$Cluster), levels = names(pal))
  
  # 1) Puntos a sf
  pts_wgs <- st_as_sf(df_case, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  
  # 2) Capas base
  if (is.null(base_layers)) base_layers <- prepare_base_layers()
  iberia_wgs_clip  <- base_layers$iberia_wgs_clip
  france_wgs_clip  <- base_layers$france_wgs_clip
  andorra_wgs_clip <- base_layers$andorra_wgs_clip
  land_wgs_clip    <- base_layers$land_wgs_clip
  bbox_wgs         <- base_layers$bbox_wgs
  crs_planar       <- base_layers$crs_planar
  
  # 3) Transformaciones
  pts     <- st_transform(pts_wgs, crs_planar)
  iberia  <- st_transform(iberia_wgs_clip, crs_planar)
  france  <- st_transform(france_wgs_clip, crs_planar)
  andorra <- st_transform(andorra_wgs_clip, crs_planar)
  land    <- st_transform(land_wgs_clip, crs_planar)
  bbox_pr <- st_transform(bbox_wgs, crs_planar)
  
  # Andorra: fallback visible si queda vacía/demasiado pequeña
  andorra <- suppressWarnings(sf::st_make_valid(andorra))
  need_fallback <- (nrow(andorra) == 0) ||
    any(sf::st_is_empty(andorra)) ||
    (as.numeric(sf::st_area(sf::st_union(andorra))) < 5e7)  # ~50 km^2
  
  if (need_fallback) {
    a_pt <- sf::st_as_sf(
      data.frame(id = 1, lon = 1.6, lat = 42.55),
      coords = c("lon", "lat"), crs = 4326
    )
    a_pt <- sf::st_transform(a_pt, crs_planar)
    andorra <- sf::st_as_sf(sf::st_buffer(a_pt, dist = 12000)) # 12 km radio
  }
  
  # 4) Voronoi y recorte al mar (mar = bbox - tierra)
  voro <- st_voronoi(st_union(pts))
  voro_sf <- st_collection_extract(voro, "POLYGON") %>% st_as_sf()
  st_crs(voro_sf) <- st_crs(pts)
  
  idx <- st_nearest_feature(voro_sf, pts)
  voro_sf <- voro_sf %>%
    mutate(Cluster = as.character(pts$Cluster[idx])) %>%
    mutate(Cluster = factor(Cluster, levels = names(pal)))
  
  sea_mask <- st_difference(bbox_pr, st_buffer(land, 0))
  voro_sea <- suppressWarnings(st_intersection(voro_sf, sea_mask))
  
  # Robustez: normalizar sf y evitar capa vacía
  if (inherits(voro_sea, "sfc")) voro_sea <- st_as_sf(voro_sea)
  if (!inherits(voro_sea, "sf")) voro_sea <- tryCatch(st_as_sf(voro_sea), error = function(e) NULL)
  if (is.null(voro_sea) || nrow(voro_sea) == 0) voro_sea <- NULL
  
  # 5) Grid hex y proporciones por celda
  hex <- st_make_grid(bbox_pr, cellsize = cellsize, what = "polygons", square = FALSE)
  hex <- st_as_sf(hex)
  st_crs(hex) <- st_crs(bbox_pr)
  hex <- st_intersection(hex, bbox_pr) %>% mutate(id = row_number())
  hex_c <- st_centroid(hex) # aviso "attributes constant..." es normal
  
  join_idx <- st_join(pts, hex, join = st_within)
  tab <- join_idx %>% st_drop_geometry() %>%
    filter(!is.na(id)) %>%
    count(id, Cluster = as.character(Cluster), name = "n")
  
  totals <- tab %>% group_by(id) %>% summarise(total = sum(n), .groups = "drop")
  
  tab_w <- tab %>%
    group_by(id, Cluster) %>%
    summarise(n = sum(n), .groups = "drop_last") %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    pivot_wider(id_cols = id, names_from = Cluster, values_from = prop, values_fill = 0)
  
  pie_cols <- names(pal)[names(pal) %in% names(tab_w)]
  
  hex_dat <- hex_c %>%
    left_join(tab_w, by = "id") %>%
    left_join(totals, by = "id")
  
  coords <- st_coordinates(st_geometry(hex_dat))
  hex_df <- hex_dat %>% st_drop_geometry() %>%
    mutate(x = coords[, 1], y = coords[, 2]) %>%
    mutate(across(all_of(pie_cols), ~ replace_na(as.numeric(.), 0))) %>%
    mutate(total = replace_na(total, 0)) %>%
    mutate(radius = 15000 + sqrt(total) * 1500)
  
  hex_df_plot <- filter(hex_df, total > 0)
  lim <- st_bbox(bbox_pr)
  
  # --- Plot: orden correcto de capas (Francia -> Voronoi mar -> Iberia/Andorra -> pies) ---
  p <- ggplot() +
    # Francia (blanco)
    geom_sf(data = france, fill = "white", color = "grey75", size = 0.25)
  
  # Voronoi (OPACO) solo si existe
  if (!is.null(voro_sea)) {
    p <- p + geom_sf(data = voro_sea, aes(fill = Cluster), alpha = voronoi_alpha, color = NA)
  }
  
  # Continente encima
  p <- p +
    geom_sf(data = iberia, fill = "white", color = "grey30", size = 0.4) +
    geom_sf(data = andorra, fill = "white", color = "grey30", size = 0.35) +
    geom_scatterpie(
      data = hex_df_plot,
      aes(x = x, y = y, r = radius),
      cols = pie_cols,
      color = NA,
      alpha = 1
    ) +
    scale_fill_manual(values = pal, name = "Cluster", drop = FALSE) +
    guides(fill = guide_legend(override.aes = list(alpha = 1))) +
    coord_sf(
      xlim = c(lim["xmin"], lim["xmax"]),
      ylim = c(lim["ymin"], lim["ymax"]),
      expand = FALSE
    ) +
    labs(
      title = title_text,
      subtitle = "Sombreado Voronoi + proporciones por celda (hex grid)",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major = element_line(color = "grey85", size = 0.2))
  
  p
}

# ============================
# 2) Preparar coordenadas y unir a comp_all
# ============================
# (Si tienes norm_key definido antes, úsalo para evitar fallos de join por mayúsculas/tildes/espacios)
if (exists("norm_key")) {
  comp_all <- comp_all %>% mutate(Locality = norm_key(Locality))
}

coords_ref <- edges_with_info %>%
  group_by(Locality_from) %>%
  summarise(
    lon = mean(x, na.rm = TRUE),
    lat = mean(y, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(Locality = Locality_from)

if (exists("norm_key")) {
  coords_ref <- coords_ref %>% mutate(Locality = norm_key(Locality))
}

comp_all <- comp_all %>%
  left_join(coords_ref, by = "Locality")

# Si el join creó lon.x/lon.y etc., los unificamos
if (any(grepl("^lon\\.", names(comp_all))) || any(grepl("^lat\\.", names(comp_all)))) {
  comp_all <- comp_all %>%
    mutate(
      lon = coalesce(lon, lon.x, lon.y),
      lat = coalesce(lat, lat.x, lat.y)
    ) %>%
    select(-matches("^lon\\.|^lat\\."))
}

stopifnot(all(c("Locality","lon","lat",
                "Cluster_inv","Cluster_cnt","Cluster_linear","Cluster_rbf") %in% names(comp_all)))

# ============================
# 3) Paleta global (consistente entre mapas)
# ============================
pal_all <- make_pal_from_user(comp_all, pal_user16, hcl_palette = "Dark 3")

# Control del recorte norte (menos Francia -> baja ymax)
base_layers <- prepare_base_layers(
  xmin = -10, ymin = 35.5, xmax = 4, ymax = 44.5,
  crs_planar = 3035
)

# ============================
# 4) Generar y guardar mapas
# ============================
p_inv <- make_map_for_case(
  comp_all, "Cluster_inv",
  "Clusters por distancia genética (Inverse 1/d)",
  base_layers = base_layers, pal = pal_all,
  voronoi_alpha = 0.25
)

p_cnt <- make_map_for_case(
  comp_all, "Cluster_cnt",
  "Clusters por número de conexiones",
  base_layers = base_layers, pal = pal_all,
  voronoi_alpha = 0.25
)

p_lin <- make_map_for_case(
  comp_all, "Cluster_linear",
  "Clusters por similitud lineal (1 - d/max)",
  base_layers = base_layers, pal = pal_all,
  voronoi_alpha = 0.25
)

p_rbf <- make_map_for_case(
  comp_all, "Cluster_rbf",
  "Clusters por similitud RBF (exp(-d/sigma))",
  base_layers = base_layers, pal = pal_all,
  voronoi_alpha = 0.25
)

print(p_inv); print(p_cnt); print(p_lin); print(p_rbf)

outdir <- "figures"
dir.create(outdir, showWarnings = FALSE)

ggsave(file.path(outdir, "map_inverse_1d_metazoa.pdf"), plot = p_inv, width = 8, height = 7)
ggsave(file.path(outdir, "map_connections_metazoa.pdf"), plot = p_cnt, width = 8, height = 7)
ggsave(file.path(outdir, "map_linear_metazoa.pdf"), plot = p_lin, width = 8, height = 7)
ggsave(file.path(outdir, "map_rbf_metazoa.pdf"), plot = p_rbf, width = 8, height = 7)

