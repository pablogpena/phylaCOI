library(dplyr)
library(igraph)
library(tidygraph)
library(ggraph)
library(leaflet)
library(RColorBrewer)
library(scales)
library(htmlwidgets)
library(factoextra)
library(geosphere)
library(ClustGeo)
library(flexclust)
library(fpc)
library(mclust)
library(mapview)
library(webshot2)
library(cluster)
library(ggplot2)
library(ggspatial)   
library(sf)
library(networkD3)
library(aricode)


setwd("C:/TEMP/FASE1/3_magia_nombre_carpeta/eKOI_metabarcoding_database_FILOS/_RESULTADOS_HAPLOTIPOS_SCRIPT_7/3_PROTISTAS_FUNGI")

# Carpeta raíz donde están las carpetas de los filos
root_dir <- "C:/TEMP/FASE1/3_magia_nombre_carpeta/eKOI_metabarcoding_database_FILOS/_RESULTADOS_HAPLOTIPOS_SCRIPT_7/3_PROTISTAS_FUNGI"

# Listar carpetas de filos
filos <- list.dirs(root_dir, recursive = FALSE)

# Inicializar listas para acumular datos
all_edges <- list()
all_points <- list()

for (filo_dir in filos) {
  filo_name <- basename(filo_dir)
  
  # Archivos CSV dentro de cada carpeta de filo
  edges_file <- file.path(filo_dir, "edges_Mi_filo.csv")
  points_file <- file.path(filo_dir, "points_Mi_filo.csv")
  
  if (file.exists(edges_file) & file.exists(points_file)) {
    message("✅ Cargando datos de: ", filo_name)
    
    edges_df <- read.csv(edges_file)
    points_df <- read.csv(points_file)
    
    # Añadir columna con el nombre del filo
    edges_df$Filo <- filo_name
    points_df$Filo <- filo_name
    
    all_edges[[filo_name]] <- edges_df
    all_points[[filo_name]] <- points_df
  } else {
    message("⚠️ No se encontraron archivos para ", filo_name)
  }
}

# Combinar todos los datos en un único dataframe
edges_df_all <- bind_rows(all_edges)
points_df_all <- bind_rows(all_points)

#
#
#
#
#-----------------------------------------------------------------------------
### 1. Network of Localities with Adjusted Weights (OTUs and Genetic Distance)
#-----------------------------------------------------------------------------

# --- Procesamiento conjunto corregido ---
edges_with_info <- edges_df_all %>%
  left_join(points_df_all %>% select(UniqueID, Locality_from = Localities, OTU_ID_from = OTU_ID, Filo_from = Filo),
            by = c("from" = "UniqueID")) %>%
  left_join(points_df_all %>% select(UniqueID, Locality_to = Localities, OTU_ID_to = OTU_ID, Filo_to = Filo),
            by = c("to" = "UniqueID")) %>%
  filter(Locality_from != Locality_to) %>%
  filter(OTU_ID_from == OTU_ID_to) %>%     # Mantener solo coincidencias
  mutate(OTU_ID = OTU_ID_from)             # Crear columna OTU_ID para agrupar

edges_locality_net <- edges_with_info %>%
  group_by(Locality_from, Locality_to, OTU_ID, Filo_from) %>%
  summarise(dist_gen_prom = mean(distancia_genetica, na.rm = TRUE), .groups = "drop") %>%
  group_by(Locality_from, Locality_to) %>%
  summarise(
    n_otus = n(),
    dist_gen_prom = mean(dist_gen_prom, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(peso = n_otus / (dist_gen_prom + 1e-6))

# Crear grafo global
g_loc <- graph_from_data_frame(
  d = edges_locality_net,
  vertices = data.frame(name = unique(c(edges_locality_net$Locality_from, edges_locality_net$Locality_to))),
  directed = FALSE
)

clusters <- cluster_louvain(g_loc, weights = E(g_loc)$peso)

locality_clusters <- data.frame(
  Locality = V(g_loc)$name,
  Cluster = clusters$membership
)

# Coordenadas globales
locality_coords <- points_df_all %>%
  group_by(Localities) %>%
  summarise(lat = mean(lat), lon = mean(lon), .groups = "drop") %>%
  left_join(locality_clusters, by = c("Localities" = "Locality"))

# Añadir coordenadas a edges
edges_mapa <- edges_locality_net %>%
  left_join(locality_coords, by = c("Locality_from" = "Localities")) %>%
  rename(lon_from = lon, lat_from = lat) %>%
  left_join(locality_coords, by = c("Locality_to" = "Localities")) %>%
  rename(lon_to = lon, lat_to = lat)

#Exportar las coordenadas de las localidades y el cluster al que pertenecen para representación gráfica en mapa
write.csv(locality_coords,
          file = "locality_coords_louvain_protist.csv",
          row.names = FALSE)
#-------------------------------------------------------------

#
#
#
#
#---------------------
### 2. K-means max=10
#---------------------
library(ggplot2)
library(dplyr)
library(scales)  # para rescale
library(RColorBrewer)
library(ggrepel)
library(ggspatial)  # para fondos de mapas


# Preparar datos para K-means
data_cluster <- locality_coords %>%
  select(Locality = Localities, lon, lat) %>%
  left_join(
    edges_locality_net %>%
      group_by(Locality_from) %>%
      summarise(
        mean_genetic_dist = mean(dist_gen_prom, na.rm = TRUE),
        total_otus = sum(n_otus, na.rm = TRUE),
        total_peso = sum(peso, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("Locality" = "Locality_from")
  ) %>%
  #na.omit() #remove all localities without connections
  mutate(
    mean_genetic_dist = ifelse(is.na(mean_genetic_dist), 0, mean_genetic_dist),
    total_otus = ifelse(is.na(total_otus), 0, total_otus),
    total_peso = ifelse(is.na(total_peso), 0, total_peso),
  )

# Escalar datos
kmeans_data <- scale(data_cluster[, c("lon", "lat", "mean_genetic_dist", "total_otus", "total_peso")])

# Elegir k automáticamente (silhouette)
# Función robusta para número óptimo de clusters con silueta
get_best_k <- function(data) {
  n <- nrow(data)
  n_unique <- nrow(unique(data))
  
  if (n_unique < 3) {
    return(1)  # menos de 3 puntos distintos -> 1 cluster
  }
  
  # límite superior = 10
  max_k <- min(n - 1, n_unique - 1, 10)
  
  sil_width <- numeric(max_k)
  
  for (k in 2:max_k) {
    km_res <- kmeans(data, centers = k, nstart = 25)
    ss <- silhouette(km_res$cluster, dist(data))
    sil_width[k] <- mean(ss[, 3])
  }
  
  best_k <- which.max(sil_width)
  return(best_k)
}

# --- Uso ---
cat("Calculando número óptimo de clusters globales con método de silueta...\n")
#best_k <- get_best_k(kmeans_data) #chose the best option in a range k=(1:10)
best_k <- 10 #to force k=10 to make the clusters comparable between groups
cat("Número óptimo de clusters (k) sugerido:", best_k, "\n")

# Ejecutar K-means
set.seed(42)
kmeans_res <- kmeans(kmeans_data, centers = best_k, nstart = 25)


# Añadir clusters a localidad global
locality_coords_kmeans <- locality_coords %>%
  select(Localities, lon, lat) %>%
  left_join(
    data_cluster %>%
      select(Locality, mean_genetic_dist) %>%
      mutate(Cluster = kmeans_res$cluster),
    by = c("Localities" = "Locality")
  )

# Preparar edges con coordenadas de origen y destino
edges_mapa_kmeans <- edges_locality_net %>%
  left_join(locality_coords_kmeans, by = c("Locality_from" = "Localities")) %>%
  rename(lon_from = lon, lat_from = lat) %>%
  left_join(locality_coords_kmeans, by = c("Locality_to" = "Localities")) %>%
  rename(lon_to = lon, lat_to = lat)

# Guardamos las coordenadas de las localidades de los clusters de Kmeans para hacer mapas a partir de este archivo
write.csv(locality_coords_kmeans, "coordenadas_kmeans_proists_fungi.csv", row.names = FALSE) #necesario para hacer el mapa con mejor calidad desde otro programa



#
#
#
#----------------------------------
### 3. ClustGeo estático kmax=10
#----------------------------------
library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)
library(sf)
library(ggspatial)
library(ggrepel)

# -----------------------------
# Preparar datos con clusters ClustGeo
# -----------------------------
# Localidades únicas globales
localidades <- locality_coords_kmeans %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  distinct(Localities, lat, lon)

# Matriz de distancias geográficas
D_geo <- distm(localidades[, c("lon", "lat")], fun = distHaversine)

# Matriz de distancias genéticas global
gen_dist <- edges_locality_net %>%
  filter(Locality_from %in% localidades$Localities,
         Locality_to %in% localidades$Localities) %>%
  select(Locality_from, Locality_to, dist_gen_prom)

locality_names <- localidades$Localities
n <- length(locality_names)
D_gen <- matrix(NA, n, n, dimnames = list(locality_names, locality_names))
for (i in 1:nrow(gen_dist)) {
  D_gen[gen_dist$Locality_from[i], gen_dist$Locality_to[i]] <- gen_dist$dist_gen_prom[i]
  D_gen[gen_dist$Locality_to[i], gen_dist$Locality_from[i]] <- gen_dist$dist_gen_prom[i]
}
D_gen[is.na(D_gen)] <- max(D_gen, na.rm = TRUE)

# Convertir a dist
D0 <- as.dist(D_gen)
D1 <- as.dist(D_geo)

# ClustGeo: selection of the optimum alpha value
range_alpha <- seq(0, 1, 0.1)
tree <- hclustgeo(D0) #Realiza un clustering jerárquico solo basado en distancias genéticas
stab <- choicealpha(D0, D1, range.alpha = range_alpha, K = 6, graph = TRUE) #Selección de alpha (peso geográfico vs genético) y gráfico
Qnorm <- stab$Qnorm
range_alpha <- stab$range.alpha
sum_qnorm <- rowSums(Qnorm)
best_idx <- which.max(sum_qnorm)
alpha_optimo <- range_alpha[best_idx]
cat("Alpha óptimo seleccionado:", alpha_optimo, "\n")
tree_opt <- hclustgeo(D0, D1, alpha = alpha_optimo)

# Seleccionar k (máx. 10)
range_k <- 3:min(10, nrow(localidades) - 1)

stab_k <- sapply(range_k, function(k){
  part <- cutree(tree_opt, k = k)
  cluster.stats(D0, part)$avg.silwidth
})

best_k_clustgeo <- best_k #choice to force the same best_k for Kmeans and make possible the comparison btw methods
cat("Número óptimo de clusters (ClustGeo):", best_k_clustgeo, "\n")

# Cortar el árbol con el k óptimo
clust_opt <- cutree(tree_opt, k = best_k_clustgeo)

# Añadir columna Cluster a las localidades
localidades_clust <- localidades %>%
  mutate(Cluster = clust_opt)

write.csv(localidades_clust, 
          file = "localidades_clustgeo.csv", 
          row.names = FALSE)
#_______________________________________________________________________________

#
#
#
#----------------------------------------
### Estadísticas intra-cluster kmeans
#----------------------------------------

# Nombre del análisis
nombre_especie <- "Global_metazoos"

# Distancias geográficas promedio intra-cluster
geo_dist_prom <- locality_coords_kmeans %>%
  select(Localities, lon, lat, Cluster) %>%
  group_by(Cluster) %>%
  summarise(
    dist_geo_promedio = {
      if (n() > 1) {
        locs <- cbind(lon, lat)
        mean(distm(locs, locs, fun = distHaversine)[upper.tri(locs, diag = FALSE)]) / 1000
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  )

# Número de OTUs, haplotipos y filos por cluster
cluster_info_extra <- points_df_all %>%
  left_join(locality_coords_kmeans %>% select(Localities, Cluster),
            by = c("Localities" = "Localities")) %>%
  filter(!is.na(Cluster)) %>%
  group_by(Cluster) %>%
  summarise(
    n_otus = n_distinct(OTU_ID),
    n_haplotipos = n_distinct(UniqueID),
    n_filos = n_distinct(Filo),
    .groups = "drop"
  )

# Estadísticas resumidas por cluster
resumen_clusters_global <- locality_coords_kmeans %>%
  group_by(Cluster) %>%
  summarise(
    n_localidades = n(),
    dist_gen_promedio = mean(mean_genetic_dist, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(geo_dist_prom, by = "Cluster") %>%
  left_join(cluster_info_extra, by = "Cluster") %>%
  mutate(
    k_optimo = n_clusters_km,
    tipo_cluster = "K-means global",
    especie = nombre_especie
  )


# Ver o guardar
print(resumen_clusters_global)

write.csv(resumen_clusters_global,
          paste0("resumen_clusters_kmeans_metazoos.csv"),
          row.names = FALSE)

#_______________________________________________________________________________

#
#
#
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ESTADÍSTICA COMPARATIVA ENTRE KMEANS vs CLUSTALGEO Kmax=10
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Asignar clusters
clusters_clustgeo <- cutree(tree_opt, k = best_k_clustgeo)
# Guardar en el dataframe
localidades$Cluster <- clusters_clustgeo

# Crear vectores de cluster por localidad
vec_kmeans   <- locality_coords_kmeans$Cluster[match(localidades$Localities, locality_coords_kmeans$Localities)]
vec_clustgeo <- localidades$Cluster

# Calcular índice de Rand ajustado
adjustedRandIndex(vec_kmeans, vec_clustgeo)

# Crear la tabla de contingencia
tabla_confusion <- table(Kmeans = vec_kmeans, ClustGeo = vec_clustgeo)

# Convertir a data.frame
tabla_df <- as.data.frame.matrix(tabla_confusion)

# Añadir nombres de fila como columna si lo deseas
tabla_df <- cbind(Kmeans = rownames(tabla_df), tabla_df)
tabla_df
# 4. Guardar como CSV
write.csv(tabla_df, "clusters_Kmeans_vs_ClustGeo_metazoos_k10.csv", row.names = FALSE)


#-----------------------------
### --- 1. Sankey Diagram ---
#-----------------------------
df_sankey <- data.frame(
  source = paste0("Kmeans_", vec_kmeans),
  target = paste0("ClustGeo_", vec_clustgeo)
)

links <- as.data.frame(table(df_sankey))
names(links) <- c("source", "target", "value")
nodes <- data.frame(name = unique(c(links$source, links$target)))

links$source <- match(links$source, nodes$name) - 1
links$target <- match(links$target, nodes$name) - 1

sankey_plot <- sankeyNetwork(Links = links, Nodes = nodes,
                             Source = "source", Target = "target",
                             Value = "value", NodeID = "name", 
                             fontSize = 12, nodeWidth = 30)
sankey_plot

# Guardar como archivo HTML interactivo
saveWidget(sankey_plot, file = "sankey_clusters_metazoos.html", selfcontained = TRUE)
# convertir a PDF
webshot("sankey_clusters_metazoos.html", "sankey_clusters_metazoos.pdf", vwidth = 1200, vheight = 900)

#-----------------------------
### 2. Métricas de evaluación
#-----------------------------
# Asegurar que no hay NA y que están alineados
valid_idx <- which(!is.na(vec_kmeans) & !is.na(vec_clustgeo))
vec_kmeans_clean <- vec_kmeans[valid_idx]
vec_clustgeo_clean <- vec_clustgeo[valid_idx]
nmi <- NMI(vec_kmeans_clean, vec_clustgeo_clean)
cat("Normalized Mutual Information (NMI):", round(nmi, 3), "\n")

# Definir los valores que ya calculaste
ari_value <- round(adjustedRandIndex(vec_kmeans_clean, vec_clustgeo_clean), 3)
nmi_value <- round(nmi, 3)
ari_value
nmi_value

# Nombre de la especie o grupo
nombre_especie <- "metazoos"

# Crear dataframe con métricas
metricas_comparacion <- data.frame(
  especie = nombre_especie,
  NMI = nmi_value,
  ARI = ari_value
)

# Guardar como CSV (puedes usar append = TRUE si haces un loop por especies)
write.csv(metricas_comparacion, "metrics_nmi_ari_metazoos.csv", row.names = FALSE)

#____________________________________________________________________________________________



