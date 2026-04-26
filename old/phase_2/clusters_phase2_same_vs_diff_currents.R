#############################################
### Part 2: SAME & DIFFERENT HAPLOTYPES
### A–D: inverse, count, linear, RBF
#############################################

# ============================
# Librerías
# ============================
library(dplyr)
library(tidyr)
library(igraph)
library(spdep)
library(aricode)
library(ggplot2)

# Reproducibilidad global
set.seed(42)

# ============================
# Supuestos de entrada
# ============================
# - edges_with_info ya existe y contiene:
#   Locality_from, Locality_to, distancia_genetica, x, y, xend, yend, group ("same"/"diff")
# - (Opcional) ocean_df con columnas: Locality, Oceans

setwd("E:/TEMP/prueba_script_definitivos/_RESULTADOS_HAPLOTIPOS_SCRIPT_7/B_different_same_haplotypes/1_same_diff_metazoa")
#locality_coords_df <- read.csv("locality_coords_metazoa.csv", header = TRUE)
locality_coords_df <- locality_coords
  
# ============================
# Helpers
# ============================

# Moran por clúster para una tabla Locality–lon–lat–Cluster
moran_table_by_cluster <- function(locality_coords_df, k = 5,
                                   jitter_amount = 1e-6,
                                   decimals_fixed = 8,
                                   tiny_thresh = 1e-12,
                                   round_I = 3) {
  stopifnot(all(c("Locality","lon","lat","Cluster") %in% names(locality_coords_df)))
  df <- locality_coords_df %>% tidyr::drop_na(lon, lat, Cluster)
  if (nrow(df) < 4 || length(unique(df$Cluster)) < 2) {
    return(data.frame(Cluster=character(), I=numeric(), p=numeric(),
                      Freq=integer(), q_FDR=numeric()))
  }
  # Vecindad kNN
  xy <- as.matrix(df[, c("lon","lat")])
  if (jitter_amount > 0) { set.seed(1); xy <- jitter(xy, amount = jitter_amount) }
  nb <- knn2nb(knearneigh(xy, k = max(1, min(k, nrow(df)-1))))
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  labs <- factor(df$Cluster)
  X <- model.matrix(~ labs - 1); colnames(X) <- levels(labs)
  
  mor_mat <- t(sapply(seq_len(ncol(X)), function(j){
    mt <- moran.test(as.numeric(X[, j]), lw, zero.policy = TRUE)
    c(I = unname(mt$estimate["Moran I statistic"]), p = mt$p.value)
  }))
  mor_df <- as.data.frame(mor_mat); mor_df$Cluster <- colnames(X)
  
  freq_df <- data.frame(Cluster = as.character(labs)) %>%
    count(Cluster, name = "Freq")
  
  out <- mor_df %>%
    select(Cluster, I, p) %>%
    left_join(freq_df, by = "Cluster") %>%
    mutate(q_FDR = p.adjust(p, method = "BH")) %>%
    arrange(desc(I))
  
  format_small <- function(x){
    ifelse(x < tiny_thresh, paste0("<", formatC(tiny_thresh, format="e", digits=0)),
           sprintf(paste0("%.", decimals_fixed, "f"), x))
  }
  
  out %>%
    mutate(I = round(I, round_I),
           p = format_small(p),
           q_FDR = format_small(q_FDR))
}

# Construye grafo localidad–localidad y clusteriza (4 modos)
analyze_case <- function(edges_sub,
                         mode = c("inv","count","linear","rbf"),
                         sigma_grid = c(0.005, 0.010, 0.015),
                         eps = 1e-6, k_moran = 5, seed = 42) {
  mode <- match.arg(mode)
  
  # 1) Colapsar a localidad–localidad
  edges_ll <- edges_sub %>%
    filter(Locality_from != Locality_to) %>%
    group_by(Locality_from, Locality_to) %>%
    summarise(
      mean_dist = mean(distancia_genetica, na.rm = TRUE),
      n_connections = n(),
      .groups = "drop"
    )
  
  if (nrow(edges_ll) == 0) {
    return(list(
      mode = mode, sigma = NA,
      g = make_empty_graph(),
      clusters = NULL, modularity = NA_real_,
      locality_coords = data.frame(Locality=character(), lon=double(), lat=double(), Cluster=integer()),
      moran = data.frame()
    ))
  }
  
  # 2) Pesos por modo
  if (mode == "count") {
    w <- edges_ll$n_connections
  } else if (mode == "inv") {
    w <- 1 / (edges_ll$mean_dist + eps)
  } else if (mode == "linear") {
    maxd <- max(edges_ll$mean_dist, na.rm = TRUE); if (maxd == 0) maxd <- 1
    w <- 1 - edges_ll$mean_dist / maxd
  } else if (mode == "rbf") {
    # Selección de sigma por modularidad
    try_grid <- lapply(sigma_grid, function(s){
      gtmp <- graph_from_data_frame(edges_ll, directed = FALSE)
      E(gtmp)$weight <- exp(- edges_ll$mean_dist / s)
      set.seed(seed)
      cltmp <- cluster_louvain(gtmp, weights = E(gtmp)$weight)
      data.frame(sigma = s,
                 modularidad = modularity(gtmp, membership(cltmp), weights = E(gtmp)$weight),
                 n_comun = length(unique(membership(cltmp))))
    })
    grid <- bind_rows(try_grid)
    best_sigma <- grid$sigma[which.max(grid$modularidad)]
    w <- exp(- edges_ll$mean_dist / best_sigma)
  }
  
  # 3) Grafo y Louvain
  g <- graph_from_data_frame(edges_ll, directed = FALSE)
  E(g)$weight <- w
  set.seed(seed)
  cl <- cluster_louvain(g, weights = E(g)$weight)
  
  # 4) Modularity
  mod <- modularity(g, membership(cl), weights = E(g)$weight)
  
  # 5) Coordenadas por localidad (usar from y to para no perder nodos)
  coords <- bind_rows(
    edges_sub %>% transmute(Locality = Locality_from, lon = x,    lat = y),
    edges_sub %>% transmute(Locality = Locality_to,   lon = xend, lat = yend)
  ) %>%
    group_by(Locality) %>%
    summarise(lon = mean(lon, na.rm = TRUE),
              lat = mean(lat, na.rm = TRUE),
              .groups = "drop")
  
  loc_clusters <- data.frame(
    Locality = names(membership(cl)),
    Cluster  = as.integer(membership(cl)),
    stringsAsFactors = FALSE
  )
  
  locality_coords <- loc_clusters %>%
    left_join(coords, by = "Locality") %>%
    relocate(lon, lat, .after = Locality)
  
  # 6) Moran por clúster
  moran_tbl <- moran_table_by_cluster(locality_coords, k = k_moran)
  
  list(
    mode = mode,
    sigma = ifelse(mode == "rbf", best_sigma, NA),
    g = g,
    clusters = cl,
    modularity = mod,
    locality_coords = locality_coords,
    moran = moran_tbl
  )
}

# ARI entre dos particiones dado locality_coords de cada una
pairwise_ari <- function(coordsA, coordsB, colA = "Cluster", colB = "Cluster") {
  A <- coordsA %>% select(Locality, !!colA) %>% rename(ClusterA = !!colA)
  B <- coordsB %>% select(Locality, !!colB) %>% rename(ClusterB = !!colB)
  comp <- inner_join(A, B, by = "Locality")
  if (nrow(comp) == 0) return(NA_real_)
  adjustedRandIndex(comp$ClusterA, comp$ClusterB)
}

# (opcional) % Mediterráneo/Atlántico por clúster
ocean_share_by_cluster <- function(coords_df, ocean_df) {
  if (!all(c("Locality","Cluster") %in% names(coords_df))) return(NULL)
  if (!all(c("Locality","Oceans") %in% names(ocean_df)))  return(NULL)
  coords_df %>%
    select(Locality, Cluster) %>%
    inner_join(ocean_df, by = "Locality") %>%
    group_by(Cluster, Oceans) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup() %>%
    select(Cluster, Oceans, pct) %>%
    tidyr::pivot_wider(names_from = Oceans, values_from = pct, values_fill = 0) %>%
    mutate(across(-Cluster, ~round(.x, 1))) %>%
    arrange(Cluster)
}
# ============================
# % por zonas de corrientes (NUEVO)
# ============================

zones_df <- locality_coords_df %>%
  mutate(
    Zone = case_when(
      lon < 0  & lat >= 37 ~ "ATL_N_UPW",
      lon < -2 & lat < 37  ~ "ATL_S_UPW_ALM",
      lon >= -2            ~ "MED",
      TRUE                 ~ NA_character_
    )
  ) %>%
  select(Locality, Zone)

# Control rápido: localidades sin asignación de zona
cat("\n=== Control de asignación de zonas ===\n")
print(table(is.na(zones_df$Zone)))

write.csv(zones_df, "locality_zones_currents.csv", row.names = FALSE)

zone_share_by_cluster <- function(coords_df, zones_df) {
  coords_df %>%
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

# ============================
# 1) Separación SAME / DIFF
# ============================
edges_same <- edges_with_info %>% filter(group == "same")
edges_diff <- edges_with_info %>% filter(group == "diff")

# ============================
# 2) Ejecutar las 4 particiones por grupo
# ============================
res_same_inv    <- analyze_case(edges_same, mode = "inv")
res_same_cnt    <- analyze_case(edges_same, mode = "count")
res_same_linear <- analyze_case(edges_same, mode = "linear")
res_same_rbf    <- analyze_case(edges_same, mode = "rbf")

res_diff_inv    <- analyze_case(edges_diff, mode = "inv")
res_diff_cnt    <- analyze_case(edges_diff, mode = "count")
res_diff_linear <- analyze_case(edges_diff, mode = "linear")
res_diff_rbf    <- analyze_case(edges_diff, mode = "rbf")

# ============================
# 3) Modularidad — resumen
# ============================
mod_summary <- bind_rows(
  data.frame(grupo = "SAME", metodo = "Inverse_1/d", modularidad = res_same_inv$modularity,
             n_comun = ifelse(is.null(res_same_inv$clusters), NA_integer_, length(unique(membership(res_same_inv$clusters))))),
  data.frame(grupo = "SAME", metodo = "Count",        modularidad = res_same_cnt$modularity,
             n_comun = ifelse(is.null(res_same_cnt$clusters), NA_integer_, length(unique(membership(res_same_cnt$clusters))))),
  data.frame(grupo = "SAME", metodo = "Linear",       modularidad = res_same_linear$modularity,
             n_comun = ifelse(is.null(res_same_linear$clusters), NA_integer_, length(unique(membership(res_same_linear$clusters))))),
  data.frame(grupo = "SAME", metodo = paste0("RBF (sigma=", sprintf("%.3f", res_same_rbf$sigma), ")"),
             modularidad = res_same_rbf$modularity,
             n_comun = ifelse(is.null(res_same_rbf$clusters), NA_integer_, length(unique(membership(res_same_rbf$clusters))))),
  data.frame(grupo = "DIFF", metodo = "Inverse_1/d", modularidad = res_diff_inv$modularity,
             n_comun = ifelse(is.null(res_diff_inv$clusters), NA_integer_, length(unique(membership(res_diff_inv$clusters))))),
  data.frame(grupo = "DIFF", metodo = "Count",        modularidad = res_diff_cnt$modularity,
             n_comun = ifelse(is.null(res_diff_cnt$clusters), NA_integer_, length(unique(membership(res_diff_cnt$clusters))))),
  data.frame(grupo = "DIFF", metodo = "Linear",       modularidad = res_diff_linear$modularity,
             n_comun = ifelse(is.null(res_diff_linear$clusters), NA_integer_, length(unique(membership(res_diff_linear$clusters))))),
  data.frame(grupo = "DIFF", metodo = paste0("RBF (sigma=", sprintf("%.3f", res_diff_rbf$sigma), ")"),
             modularidad = res_diff_rbf$modularity,
             n_comun = ifelse(is.null(res_diff_rbf$clusters), NA_integer_, length(unique(membership(res_diff_rbf$clusters)))))
)
cat("\n=== MODULARIDAD (SAME vs DIFF, 4 métodos) ===\n"); print(mod_summary, row.names = FALSE)
write.csv(mod_summary, "modularity_same_diff_summary_metazoa.csv", row.names = FALSE)

# ============================
# 4) ARI — (a) dentro de SAME; (b) dentro de DIFF; (c) SAME vs DIFF por método
# ============================
# (a) SAME: ARI entre sus 4 particiones
ari_same <- data.frame(
  Comparacion = c("Inv vs Count", "Inv vs Linear", "Inv vs RBF",
                  "Count vs Linear", "Count vs RBF", "Linear vs RBF"),
  ARI = c(
    pairwise_ari(res_same_inv$locality_coords,    res_same_cnt$locality_coords),
    pairwise_ari(res_same_inv$locality_coords,    res_same_linear$locality_coords),
    pairwise_ari(res_same_inv$locality_coords,    res_same_rbf$locality_coords),
    pairwise_ari(res_same_cnt$locality_coords,    res_same_linear$locality_coords),
    pairwise_ari(res_same_cnt$locality_coords,    res_same_rbf$locality_coords),
    pairwise_ari(res_same_linear$locality_coords, res_same_rbf$locality_coords)
  )
)

# (b) DIFF: ARI entre sus 4 particiones
ari_diff <- data.frame(
  Comparacion = c("Inv vs Count", "Inv vs Linear", "Inv vs RBF",
                  "Count vs Linear", "Count vs RBF", "Linear vs RBF"),
  ARI = c(
    pairwise_ari(res_diff_inv$locality_coords,    res_diff_cnt$locality_coords),
    pairwise_ari(res_diff_inv$locality_coords,    res_diff_linear$locality_coords),
    pairwise_ari(res_diff_inv$locality_coords,    res_diff_rbf$locality_coords),
    pairwise_ari(res_diff_cnt$locality_coords,    res_diff_linear$locality_coords),
    pairwise_ari(res_diff_cnt$locality_coords,    res_diff_rbf$locality_coords),
    pairwise_ari(res_diff_linear$locality_coords, res_diff_rbf$locality_coords)
  )
)

# (c) SAME vs DIFF — ARI por método
ari_same_vs_diff <- data.frame(
  Metodo = c("Inverse_1/d", "Count", "Linear", "RBF"),
  ARI = c(
    pairwise_ari(res_same_inv$locality_coords,    res_diff_inv$locality_coords),
    pairwise_ari(res_same_cnt$locality_coords,    res_diff_cnt$locality_coords),
    pairwise_ari(res_same_linear$locality_coords, res_diff_linear$locality_coords),
    pairwise_ari(res_same_rbf$locality_coords,    res_diff_rbf$locality_coords)
  )
)

cat("\n=== ARI (SAME — entre particiones) ===\n"); print(ari_same, row.names = FALSE)
cat("\n=== ARI (DIFF — entre particiones) ===\n"); print(ari_diff, row.names = FALSE)
cat("\n=== ARI (SAME vs DIFF — por método) ===\n"); print(ari_same_vs_diff, row.names = FALSE)

write.csv(ari_same,         "ari_same_within_methods_metazoa.csv", row.names = FALSE)
write.csv(ari_diff,         "ari_diff_within_methods_metazoa.csv", row.names = FALSE)
write.csv(ari_same_vs_diff, "ari_same_vs_diff_by_method_metazoa.csv", row.names = FALSE)

# ============================
# 5) Moran’s I — ya viene dentro de cada resultado
# ============================
cat("\n=== Moran por clúster — SAME (Inverse) ===\n"); print(res_same_inv$moran, row.names = FALSE)
cat("\n=== Moran por clúster — SAME (Count)   ===\n"); print(res_same_cnt$moran, row.names = FALSE)
cat("\n=== Moran por clúster — SAME (Linear)  ===\n"); print(res_same_linear$moran, row.names = FALSE)
cat("\n=== Moran por clúster — SAME (RBF)     ===\n"); print(res_same_rbf$moran, row.names = FALSE)

cat("\n=== Moran por clúster — DIFF (Inverse) ===\n"); print(res_diff_inv$moran, row.names = FALSE)
cat("\n=== Moran por clúster — DIFF (Count)   ===\n"); print(res_diff_cnt$moran, row.names = FALSE)
cat("\n=== Moran por clúster — DIFF (Linear)  ===\n"); print(res_diff_linear$moran, row.names = FALSE)
cat("\n=== Moran por clúster — DIFF (RBF)     ===\n"); print(res_diff_rbf$moran, row.names = FALSE)

# Guardar Moran (tablas)
write.csv(res_same_inv$moran,    "moran_same_inverse_metazoa.csv", row.names = FALSE)
write.csv(res_same_cnt$moran,    "moran_same_count_metazoa.csv",   row.names = FALSE)
write.csv(res_same_linear$moran, "moran_same_linear_metazoa.csv",  row.names = FALSE)
write.csv(res_same_rbf$moran,    "moran_same_rbf_metazoa.csv",     row.names = FALSE)

write.csv(res_diff_inv$moran,    "moran_diff_inverse_metazoa.csv", row.names = FALSE)
write.csv(res_diff_cnt$moran,    "moran_diff_count_metazoa.csv",   row.names = FALSE)
write.csv(res_diff_linear$moran, "moran_diff_linear_metazoa.csv",  row.names = FALSE)
write.csv(res_diff_rbf$moran,    "moran_diff_rbf_metazoa.csv",     row.names = FALSE)

# ============================
# 6) Exportar locality_coords (para mapas)
# ============================
#write.csv(res_same_inv$locality_coords,    "locality_coords_same_inverse.csv", row.names = FALSE)
#write.csv(res_same_cnt$locality_coords,    "locality_coords_same_count.csv",   row.names = FALSE)
#write.csv(res_same_linear$locality_coords, "locality_coords_same_linear.csv",  row.names = FALSE)
#write.csv(res_same_rbf$locality_coords,    "locality_coords_same_rbf.csv",     row.names = FALSE)

#write.csv(res_diff_inv$locality_coords,    "locality_coords_diff_inverse.csv", row.names = FALSE)
#write.csv(res_diff_cnt$locality_coords,    "locality_coords_diff_count.csv",   row.names = FALSE)
#write.csv(res_diff_linear$locality_coords, "locality_coords_diff_linear.csv",  row.names = FALSE)
#write.csv(res_diff_rbf$locality_coords,    "locality_coords_diff_rbf.csv",     row.names = FALSE)

# ============================
# 7) % Mediterráneo / Atlántico por clúster
# ============================
if (exists("ocean_df")) {
  same_inv_ocean    <- ocean_share_by_cluster(res_same_inv$locality_coords,    ocean_df)
  same_cnt_ocean    <- ocean_share_by_cluster(res_same_cnt$locality_coords,    ocean_df)
  same_linear_ocean <- ocean_share_by_cluster(res_same_linear$locality_coords, ocean_df)
  same_rbf_ocean    <- ocean_share_by_cluster(res_same_rbf$locality_coords,    ocean_df)
  
  diff_inv_ocean    <- ocean_share_by_cluster(res_diff_inv$locality_coords,    ocean_df)
  diff_cnt_ocean    <- ocean_share_by_cluster(res_diff_cnt$locality_coords,    ocean_df)
  diff_linear_ocean <- ocean_share_by_cluster(res_diff_linear$locality_coords, ocean_df)
  diff_rbf_ocean    <- ocean_share_by_cluster(res_diff_rbf$locality_coords,    ocean_df)
  
  write.csv(same_inv_ocean,    "ocean_share_same_inverse_metazoa.csv", row.names = FALSE)
  write.csv(same_cnt_ocean,    "ocean_share_same_count_metazoa.csv",   row.names = FALSE)
  write.csv(same_linear_ocean, "ocean_share_same_linear_metazoa.csv",  row.names = FALSE)
  write.csv(same_rbf_ocean,    "ocean_share_same_rbf_metazoa.csv",     row.names = FALSE)

  write.csv(diff_inv_ocean,    "ocean_share_diff_inverse_metazoa.csv", row.names = FALSE)
  write.csv(diff_cnt_ocean,    "ocean_share_diff_count_metazoa.csv",   row.names = FALSE)
  write.csv(diff_linear_ocean, "ocean_share_diff_linear_metazoa.csv",  row.names = FALSE)
  write.csv(diff_rbf_ocean,    "ocean_share_diff_rbf_metazoa.csv",     row.names = FALSE)
}

# ============================
# 8) % ZONAS DE CORRIENTES / FRENTES POR CLÚSTER  (NUEVO)
# ============================
same_inv_zone    <- zone_share_by_cluster(res_same_inv$locality_coords,    zones_df)
same_cnt_zone    <- zone_share_by_cluster(res_same_cnt$locality_coords,    zones_df)
same_linear_zone <- zone_share_by_cluster(res_same_linear$locality_coords, zones_df)
same_rbf_zone    <- zone_share_by_cluster(res_same_rbf$locality_coords,    zones_df)

diff_inv_zone    <- zone_share_by_cluster(res_diff_inv$locality_coords,    zones_df)
diff_cnt_zone    <- zone_share_by_cluster(res_diff_cnt$locality_coords,    zones_df)
diff_linear_zone <- zone_share_by_cluster(res_diff_linear$locality_coords, zones_df)
diff_rbf_zone    <- zone_share_by_cluster(res_diff_rbf$locality_coords,    zones_df)

write.csv(same_inv_zone,    "zone_share_same_inverse_metazoa.csv", row.names = FALSE)
write.csv(same_cnt_zone,    "zone_share_same_count_metazoa.csv",   row.names = FALSE)
write.csv(same_linear_zone, "zone_share_same_linear_metazoa.csv",  row.names = FALSE)
write.csv(same_rbf_zone,    "zone_share_same_rbf_metazoa.csv",     row.names = FALSE)

write.csv(diff_inv_zone,    "zone_share_diff_inverse_metazoa.csv", row.names = FALSE)
write.csv(diff_cnt_zone,    "zone_share_diff_count_metazoa.csv",   row.names = FALSE)
write.csv(diff_linear_zone, "zone_share_diff_linear_metazoa.csv",  row.names = FALSE)
write.csv(diff_rbf_zone,    "zone_share_diff_rbf_metazoa.csv",     row.names = FALSE)

################################################
### Part 3: GEOGRAPHICAL DISTANCES OF CONNECTIONS 
### All, same and different haplotypes + OTUs
################################################
library(dplyr)
library(geosphere)   # distHaversine (m)

# -----------------------------
# 0) Base: calcula distancia geográfica (km)
# -----------------------------
edges_base <- edges_with_info %>%
  mutate(
    geo_dist_m  = geosphere::distHaversine(cbind(x, y), cbind(xend, yend)),
    geo_dist_km = geo_dist_m / 1000
  )

# Excluir enlaces intra-localidad o d=0
edges_base <- edges_base %>% filter(Locality_from != Locality_to)
edges_base <- edges_base %>% filter(geo_dist_km > 0)

# -----------------------------
# 1) Medianas globales (km)
# -----------------------------
med_geo_all  <- median(edges_base$geo_dist_km, na.rm = TRUE)
med_geo_diff <- edges_base %>%
  filter(group == "diff") %>% summarise(m = median(geo_dist_km, na.rm = TRUE)) %>% pull(m)
med_geo_same <- edges_base %>%
  filter(group == "same") %>% summarise(m = median(geo_dist_km, na.rm = TRUE)) %>% pull(m)

cat("Medianas geográficas globales (km):\n")
print(c(all = med_geo_all, diff = med_geo_diff, same = med_geo_same))

# -----------------------------
# 2) Medianas por OTU_ID (km)
# -----------------------------
med_geo_by_otu <- edges_base %>%
  group_by(OTU_ID) %>%
  summarise(
    n_edges = n(),
    n_diff  = sum(group == "diff"),
    n_same  = sum(group == "same"),
    med_km_all  = median(geo_dist_km, na.rm = TRUE),
    med_km_diff = ifelse(n_diff > 0,
                         median(geo_dist_km[group == "diff"], na.rm = TRUE),
                         NA_real_),
    med_km_same = ifelse(n_same > 0,
                         median(geo_dist_km[group == "same"], na.rm = TRUE),
                         NA_real_),
    .groups = "drop"
  )

cat("\nPrimeras filas (medianas geográficas por OTU_ID, km):\n")
print(head(med_geo_by_otu, 10))

# -----------------------------
# 3) “Mediana de las medianas” a nivel OTU (km)
# -----------------------------
resumen_geo_otu <- med_geo_by_otu %>%
  summarise(
    n_OTUs             = n(),
    n_OTUs_con_diff    = sum(!is.na(med_km_diff)),
    n_OTUs_con_same    = sum(!is.na(med_km_same)),
    med_of_meds_all_km  = median(med_km_all,  na.rm = TRUE),
    med_of_meds_diff_km = median(med_km_diff, na.rm = TRUE),
    med_of_meds_same_km = median(med_km_same, na.rm = TRUE)
  )

cat("\nResumen por OTU (mediana de las medianas, km):\n")
print(resumen_geo_otu)

# (opcional) exportar
write.csv(med_geo_by_otu, "medianas_geo_por_OTU_metazoa.csv", row.names = FALSE)
write.csv(resumen_geo_otu, "resumen_medianas_geo_por_OTU_metazoa.csv", row.names = FALSE)

# -----------------------------
# 4) Preparar datos para el análisis (crear dist_long)
# -----------------------------

# Seleccionamos columnas relevantes
edges_geo <- edges_base %>%
  select(group, geo_dist_km)

# Crear subconjuntos: all / same / diff
df_all  <- edges_geo %>% mutate(tipo = "all")
df_same <- edges_geo %>% filter(group == "same") %>% mutate(tipo = "same")
df_diff <- edges_geo %>% filter(group == "diff") %>% mutate(tipo = "diff")

# Combinar todo en formato largo
dist_long <- bind_rows(df_all, df_same, df_diff) %>%
  mutate(tipo = factor(tipo, levels = c("all", "same", "diff")))

# (opcional) Resumen para gráficos (mediana + IQR)
sum_stats <- dist_long %>%
  group_by(tipo) %>%
  summarise(
    n = n(),
    median_km = median(geo_dist_km, na.rm = TRUE),
    q1 = quantile(geo_dist_km, 0.25, na.rm = TRUE),
    q3 = quantile(geo_dist_km, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# 5) Análisis estadísticos de distancias geográficas
# Comparación entre grupos same y diff (t de Welch)
# -----------------------------
library(dplyr)
library(rstatix)    # cohens_d / levene_test

# Datos para el test
dist_test <- dist_long %>%
  filter(tipo %in% c("same", "diff")) %>%
  droplevels()

# (opcional) Comprobar homogeneidad de varianzas
var_test <- rstatix::levene_test(geo_dist_km ~ tipo, data = dist_test)
print(var_test)
# Nota: Usaremos t de Welch igualmente (var.equal = FALSE)

# t de Welch en escala original
t_res <- t.test(geo_dist_km ~ tipo, data = dist_test, var.equal = FALSE)
print(t_res)

# Tamaño del efecto (Cohen's d con corrección de Hedges)
d_res <- rstatix::cohens_d(geo_dist_km ~ tipo, data = dist_test,
                           var.equal = FALSE, hedges.correction = TRUE)
print(d_res)

# (opcional) t de Welch sobre log1p si la distribución está muy sesgada
dist_test <- dist_test %>% mutate(geo_log1p = log1p(geo_dist_km))
t_res_log <- t.test(geo_log1p ~ tipo, data = dist_test, var.equal = FALSE)
d_res_log <- rstatix::cohens_d(geo_log1p ~ tipo, data = dist_test,
                               var.equal = FALSE, hedges.correction = TRUE)

# Resumen y guardado
stats_summary <- dist_test %>%
  group_by(tipo) %>%
  summarise(
    n = n(),
    median_km = median(geo_dist_km, na.rm = TRUE),
    IQR_km = IQR(geo_dist_km, na.rm = TRUE),
    mean_km = mean(geo_dist_km, na.rm = TRUE)
  )
print(stats_summary)

write.csv(stats_summary, "estadisticos_dist_geograficas.csv", row.names = FALSE)
#saveRDS(list(t_res = t_res, d_res = d_res,
#             t_res_log = t_res_log, d_res_log = d_res_log,
#             var_test = var_test),
#        "tests_t_student_results.rds")

# -----------------------------
# 6) Geographic distance visualization
# -----------------------------
library(ggplot2)

#Paleta de colores
pal <- c(
  all  = "#bdbdbd",  # gris medio neutro
  same = "#756bb1",  # violeta/morado suave
  diff = "#de2d26"   # rojo oscuro
)

# Etiqueta de p y estrellas según t de Welch (escala original)
p_value_t <- t_res$p.value
signif_label_t <- dplyr::case_when(
  p_value_t < 0.001 ~ "***",
  p_value_t < 0.01  ~ "**",
  p_value_t < 0.05  ~ "*",
  TRUE              ~ "ns"
)

# Altura para la etiqueta
y_top <- max(dist_long$geo_dist_km, na.rm = TRUE)
y_lab <- y_top * 1.05

# ---------------------------------------
# Gráfico principal
# ---------------------------------------
p_dist <- ggplot() +
  geom_jitter(
    data = dist_long,
    aes(x = tipo, y = geo_dist_km, color = tipo),
    width = 0.15, alpha = 0.10, size = 1, show.legend = FALSE
  ) +
  geom_col(
    data = sum_stats,
    aes(x = tipo, y = median_km, fill = tipo),
    width = 0.6, alpha = 0.8, color = NA
  ) +
  geom_errorbar(
    data = sum_stats,
    aes(x = tipo, ymin = q1, ymax = q3),
    width = 0.16, linewidth = 0.7
  ) +
  # Anotación entre same (x=2) y diff (x=3)
  annotate(
    "text",
    x = 2.5, y = y_lab,
    label = paste0("Welch t-test: p = ", signif(p_value_t, 3),
                   " (", signif_label_t, ")"),
    size = 4.2
  ) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_color_manual(values = pal, guide = "none") +
  labs(
    title = "Distancias geográficas de las conexiones",
    subtitle = "Barras = mediana; líneas = IQR (Q1–Q3); puntos = todas las conexiones\nSignificancia: Welch t-test (same vs diff)",
    x = "Conjunto de conexiones",
    y = "Distancia geográfica (km)"
  ) +
  expand_limits(y = y_lab * 1.05) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank())

# Mostrar y guardar gráfico
print(p_dist)

ggsave(
  filename = "p_dist_metazoa.pdf",
  plot = p_dist,
  width = 8, height = 6,
  device = cairo_pdf
)

################################################################################
#############################################
### FIGURE: OCEANOGRAPHIC ZONES × ALL METHODS
### SAME / DIFF 
#############################################
#############################################
### FIGURE: OCEANOGRAPHIC ZONES × ALL METHODS
### ALL / SAME / DIFF  + ALL GROUPS
#############################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(grid)   # unit()

# ============================
# 0) Directorio de trabajo
# ============================
setwd("E:/TEMP/prueba_script_definitivos/_RESULTADOS_HAPLOTIPOS_SCRIPT_7/B_different_same_haplotypes/FIGURA_CORRIENTES")

# ============================
# 1) Leer todos los CSV
# ============================
files <- list.files(
  pattern = "^zone_share_.*\\.csv$",
  full.names = TRUE
)
stopifnot(length(files) > 0)

# ============================
# 2) Función robusta para parsear nombre
#    Espera: zone_share_<type>_<method>_<group>.csv
#    donde <group> puede contener underscores (se pegan)
# ============================
parse_meta_from_filename <- function(f) {
  base <- basename(f) %>%
    str_remove("^zone_share_") %>%
    str_remove("\\.csv$")
  
  parts <- str_split(base, "_", simplify = TRUE)
  
  # Si no cumple mínimo, devolvemos NA para no romper
  if (ncol(parts) < 3) {
    return(list(type = NA_character_, method = NA_character_, group = NA_character_))
  }
  
  type   <- parts[1]
  method <- parts[2]
  group  <- paste(parts[1, 3:ncol(parts)], collapse = "_")  # todo lo demás
  
  list(type = type, method = method, group = group)
}

# ============================
# 3) Cargar y anotar metadatos
# ============================
all_data <- lapply(files, function(f) {
  
  df <- read.csv(f, check.names = FALSE)
  
  info <- parse_meta_from_filename(f)
  
  df %>%
    mutate(
      type   = info$type,    # all / same / diff
      method = info$method,  # inverse / count / linear / rbf
      group  = info$group    # metazoa / arqueaplastidia / protists / multicellular / unicellular / ...
    )
}) %>% bind_rows()

# (Opcional) si quieres asegurarte de que las columnas existen
needed_cols <- c("Cluster", "ATL_N_UPW", "ATL_S_UPW_ALM", "MED", "type", "method", "group")
missing_cols <- setdiff(needed_cols, colnames(all_data))
if (length(missing_cols) > 0) {
  stop("Faltan columnas en algún CSV: ", paste(missing_cols, collapse = ", "))
}

# ============================
# 4) Pasar a formato largo (zonas)
# ============================
long_data <- all_data %>%
  pivot_longer(
    cols = c("ATL_N_UPW", "ATL_S_UPW_ALM", "MED"),
    names_to = "Zone",
    values_to = "pct"
  )

# ============================
# 5) Ordenar factores
#    (incluimos los grupos nuevos + 'all')
# ============================
type_levels   <- c("all", "same", "diff")
method_levels <- c("inverse", "count", "linear", "rbf")
group_levels  <- c("metazoa", "arqueaplastidia", "protists", "multicellular", "unicellular")

long_data <- long_data %>%
  mutate(
    method = factor(method, levels = method_levels),
    type   = factor(type,   levels = type_levels),
    # Mantiene los conocidos primero, y deja el resto al final si existieran
    group  = factor(group, levels = c(group_levels, setdiff(sort(unique(group)), group_levels))),
    Zone   = factor(Zone, levels = c("ATL_N_UPW", "ATL_S_UPW_ALM", "MED")),
    Cluster = factor(Cluster)
  )

# ============================
# 6) Figura: barras apiladas por clúster
# ============================
p <- ggplot(long_data, aes(x = Cluster, y = pct, fill = Zone)) +
  geom_col(width = 0.9) +
  facet_grid(
    group ~ method + type,
    scales = "free_x",
    space  = "free_x"
  ) +
  scale_fill_manual(
    values = c(
      ATL_N_UPW     = "#4575b4",
      ATL_S_UPW_ALM = "#91bfdb",
      MED           = "#8DD3C7"
    ),
    name = "Oceanographic zone"
  ) +
  labs(
    x = "Genetic cluster",
    y = "% of localities"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(size = 10),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    panel.spacing    = unit(0.4, "lines"),
    legend.position  = "right"
  )

print(p)

# ============================
# 7) Guardar figura
# ============================
ggsave(
  filename = "Figure_currents_all_methods_ALL_same_diff_ALLgroups.png",
  plot     = p,
  width    = 18,
  height   = 11,
  dpi      = 300
)
