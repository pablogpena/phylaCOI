# Dependencias
library(sf)
library(dplyr)
library(clue)   # para solve_LSAP (algoritmo húngaro)
library(units)

# ---------------- 1) Cargar datos ----------------
setwd("C:/TEMP/FASE1/3_magia_nombre_carpeta/eKOI_metabarcoding_database_FILOS/_RESULTADOS_HAPLOTIPOS_SCRIPT_7")

# --- 1) Rutas a tus CSVs (ajusta) ---
files <- c(
  metazoos        = "coordenadas_kmeans_metazoos.csv",
  archaeaplastidia = "coordenadas_kmeans_arqueaplastidia.csv",
  protistas       = "coordenadas_kmeans_protists_fungi.csv"
)

# --- 2) Leer dataframes ---
dfs <- lapply(files, function(f) read.csv(f, stringsAsFactors = FALSE))
names(dfs) <- names(files)
dfs <- lapply(dfs, function(df){ df$Cluster <- as.character(df$Cluster); df })

# --- 3) CRS métrico ---
crs_planar <- 3035

# --- 4) Centroides de clusters ---
centroids_list <- lapply(dfs, function(df){
  pts <- st_as_sf(df, coords = c("lon","lat"), crs = 4326, remove = FALSE)
  pts <- st_transform(pts, crs_planar)
  cent <- pts %>%
    group_by(Cluster) %>%
    summarise(geometry = st_centroid(st_combine(geometry)), .groups = "drop") %>%
    st_as_sf()
  cent$Cluster <- as.character(cent$Cluster)
  cent
})

# --- 5) Usar metazoos como referencia ---
ref_name <- "metazoos"
ref_centroids <- centroids_list[[ref_name]]


# --- 6) Asignar clusters 1:1 con Hungarian algorithm ---
for(nm in names(dfs)){
  if(nm == ref_name){
    dfs[[nm]]$Cluster_geo <- dfs[[nm]]$Cluster
    next
  }
  this_cent <- centroids_list[[nm]]
  # matriz de distancias (rows: this, cols: ref) SIN unidades
  dmat <- st_distance(this_cent, ref_centroids) %>% drop_units() %>% as.matrix()
  assignment <- solve_LSAP(dmat)  # Hungarian: minimiza distancias totales
  mapped_ref <- ref_centroids$Cluster[assignment]
  map_df <- data.frame(Cluster = this_cent$Cluster,
                       Cluster_geo = mapped_ref,
                       stringsAsFactors = FALSE)
  dfs[[nm]] <- left_join(dfs[[nm]], map_df, by = "Cluster")
}


# --- 7) Guardar CSVs alineados ---
for(nm in names(dfs)){
  out <- dfs[[nm]]
  out$Cluster <- out$Cluster_geo
  out$Cluster_geo <- NULL
  outfile <- paste0("aligned_", basename(files[nm]))
  write.csv(out, outfile, row.names = FALSE)
  message("Escrito: ", outfile)
}
