#Cargar paquetes

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
library(geosphere)
library(leaflet)
library(RColorBrewer)
library(flexclust) # para adjustedRandIndex
library(fpc)
library(mclust)
library(mapview)
library(webshot2)  # instalar con install.packages("webshot2")
library(cluster)
library(networkD3)
library(mclust)
library(aricode)


#Establecer directorio de trabajo: input file
setwd("C:/TEMP/FASE1/3_magia_nombre_carpeta/eKOI_metabarcoding_database_FILOS/_RESULTADOS_HAPLOTIPOS_SCRIPT_7/1_METAZOOS")
# Carpeta raíz donde están las carpetas de los filos
root_dir <- "C:/TEMP/FASE1/3_magia_nombre_carpeta/eKOI_metabarcoding_database_FILOS/_RESULTADOS_HAPLOTIPOS_SCRIPT_7/1_METAZOOS"

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


#Establecer directorio de trabajo:output files
setwd("C:/TEMP/FASE1/3_magia_nombre_carpeta/eKOI_metabarcoding_database_FILOS/_RESULTADOS_HAPLOTIPOS_SCRIPT_7/1_METAZOOS/prueba")


#-------------------------------------
### 1.clustering method: K-means
# (does not include geographical data) 
#-------------------------------------

# 1. Preparar datos
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

# 2. Escalar datos
kmeans_data <- scale(data_cluster[, c("lon", "lat", "mean_genetic_dist", "total_otus", "total_peso")])

# 3. Elegir k automáticamente (silhouette)
# Función robusta para número óptimo de clusters con silueta
get_best_k <- function(data) {
  n <- nrow(data)
  n_unique <- nrow(unique(data))
  
  if (n_unique < 3) {
    return(1)  # menos de 3 puntos distintos -> 1 cluster
  }
  
  max_k <- min(n - 1, n_unique - 1)
  
  sil_width <- numeric(max_k)
  
  for (k in 2:max_k) {
    km_res <- kmeans(data, centers = k, nstart = 25)
    ss <- silhouette(km_res$cluster, dist(data))
    sil_width[k] <- mean(ss[, 3])
  }
  
  best_k <- which.max(sil_width)
  return(best_k)
}

cat("Calculando número óptimo de clusters globales con método de silueta...\n")
best_k <- get_best_k(kmeans_data)
cat("Número óptimo de clusters (k) sugerido:", best_k, "\n")

# 4. Ejecutar K-means
set.seed(42)
kmeans_res <- kmeans(kmeans_data, centers = best_k, nstart = 25)

### Si hay demasiados clusters para ver patrones de distribución, 
# acotamos el número de K óptimo máximo: 

###+++++++++++
# Kmeans
# K MAX = 10
###+++++++++++

# 1. Preparar datos para K-means a nivel global
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

# 2. Escalar datos
kmeans_data <- scale(data_cluster[, c("lon", "lat", "mean_genetic_dist", "total_otus", "total_peso")])

# 3. Elegir k automáticamente (silhouette)
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

cat("Calculando número óptimo de clusters globales con método de silueta...\n")
best_k <- get_best_k(kmeans_data) 
cat("Número óptimo de clusters (k) sugerido:", best_k, "\n")

# 4. Ejecutar K-means
set.seed(42)
kmeans_res <- kmeans(kmeans_data, centers = best_k, nstart = 25)

# 5. Añadir clusters a localidad global
locality_coords_kmeans <- locality_coords %>%
  select(Localities, lon, lat) %>%
  left_join(
    data_cluster %>%
      select(Locality, mean_genetic_dist) %>%
      mutate(Cluster = kmeans_res$cluster),
    by = c("Localities" = "Locality")
  )

# 6. Preparar edges con coordenadas de origen y destino
edges_mapa_kmeans <- edges_locality_net %>%
  left_join(locality_coords_kmeans, by = c("Locality_from" = "Localities")) %>%
  rename(lon_from = lon, lat_from = lat) %>%
  left_join(locality_coords_kmeans, by = c("Locality_to" = "Localities")) %>%
  rename(lon_to = lon, lat_to = lat)

### Mapa interactivo
# 7. Crear paleta por cluster
n_clusters_km <- length(unique(na.omit(locality_coords_kmeans$Cluster)))
pal_km <- colorFactor(brewer.pal(max(3, min(n_clusters_km, 8)), "Set2"),
                      domain = locality_coords_kmeans$Cluster)

# 8. Mapa interactivo global K-means
map_kmeans <- leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolylines(data = edges_mapa_kmeans,
               lng = ~c(lon_from, lon_to),
               lat = ~c(lat_from, lat_to),
               group = "Connections",
               color = "gray40",
               weight = ~rescale(peso, to = c(0.5, 4)),
               opacity = 0.6,
               label = ~paste0("From: ", Locality_from, " ↔ To: ", Locality_to,
                               "<br>Connected by ", n_otus, " OTUs",
                               "<br>Mean genetic distance: ", round(dist_gen_prom, 4),
                               "<br>Weight: ", round(peso, 3)))

for (cl in sort(unique(na.omit(locality_coords_kmeans$Cluster)))) {
  cluster_data <- locality_coords_kmeans %>% filter(Cluster == cl)
  
  map_kmeans <- map_kmeans %>%
    addCircleMarkers(data = cluster_data,
                     lng = ~lon, lat = ~lat,
                     radius = 6,
                     fillColor = ~pal_km(Cluster),
                     color = "black", weight = 1,
                     fillOpacity = 0.8,
                     stroke = TRUE,
                     label = ~paste0("Locality: ", Localities, "<br>Cluster (KM): ", Cluster),
                     group = paste0("Cluster ", cl))
}

map_kmeans <- map_kmeans %>%
  addLegend("bottomright", pal = pal_km, values = locality_coords_kmeans$Cluster,
            title = "Clusters K-means", opacity = 1) #%>%
#addLayersControl(
#overlayGroups = c("Connections", paste0("Cluster ", sort(unique(na.omit(locality_coords_kmeans$Cluster))))),
#options = layersControlOptions(collapsed = FALSE)
#)

# Guardar como HTML
saveWidget(map_kmeans, file = "mapa_kmeans_metazoos_k4.html", selfcontained = TRUE)
##Guardar el mapa en HTML para .png
saveWidget(map_kmeans, "mapa_kmeans_metazoos_k4_2.html", selfcontained = TRUE)
#Generar ruta absoluta correcta
ruta_html <- normalizePath("mapa_kmeans_metazoos_k4_2.html", winslash = "/")
#Exportar a PNG con webshot2 (forzando el paquete)
webshot2::webshot(
  ruta_html,
  "mapa_kmeans_metazoos_k4.png",
  vwidth = 1600,    # ancho más grande
  vheight = 1200,   # alto más grande
  delay = 120        # espera a que todo se cargue
)


#----------------------------------------
### Estadísticas intra-cluster (multi-filo global)
#----------------------------------------

# Nombre del análisis
nombre_especie <- "Global_metazoa"

# 1. Distancias geográficas promedio intra-cluster
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

# 2. Número de OTUs, haplotipos y filos por cluster
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

# 3. Estadísticas resumidas por cluster
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


# 4. Ver o guardar
print(resumen_clusters_global)

write.csv(resumen_clusters_global,
          paste0("resumen_clusters_kmeans_metazoa.csv"),
          row.names = FALSE)



#
#
#
#
#-------------------------------------
### 2. Clustering method: ClustGeo
# (does include geographical data)
#-------------------------------------

# 1. Localidades únicas globales
localidades <- locality_coords_kmeans %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  distinct(Localities, lat, lon)

# 2. Matriz de distancias geográficas
D_geo <- distm(localidades[, c("lon", "lat")], fun = distHaversine)

# 3. Matriz de distancias genéticas global
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

# 4. Convertir a dist
D0 <- as.dist(D_gen)
D1 <- as.dist(D_geo)

# 5. ClustGeo: selection of the optimum alpha value
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

# 6. Seleccionar k
#range_k <- 2:(nrow(localidades) - 1)
#stab_k <- sapply(range_k, function(k){
#  part <- cutree(tree_opt, k = k)
#  cluster.stats(D0, part)$avg.silwidth
#})
#best_k_clustgeo <- range_k[which.max(stab_k)]
best_k_clustgeo <- best_k #forzar a que tenga el mismo numero de clusters que el método Kmeans para poder comparar los resultados
cat("Número óptimo de clusters (ClustGeo):", best_k_clustgeo, "\n")

# 7. Asignar clusters
localidades$Cluster <- cutree(tree_opt, k = best_k_clustgeo)

###Mapa interactivo
# 8. Mapa interactivo ClustGeo
pal <- colorFactor(brewer.pal(best_k_clustgeo, "Set2"), domain = localidades$Cluster)
mapa <- leaflet(localidades) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircleMarkers(
    ~lon, ~lat,
    color = "black",
    fillColor = ~pal(Cluster),
    fillOpacity = 0.8,
    weight = 1,
    radius = 6,
    label = ~paste0(Localities, "<br>Cluster: ", Cluster)
  ) %>%
  addLegend("bottomright", pal = pal, values = localidades$Cluster,
            title = "Cluster ClustGeo", opacity = 1)
saveWidget(mapa, "mapa_clusters_clustgeo_fungi.html", selfcontained = TRUE)


#Generar ruta absoluta correcta
ruta_html <- normalizePath("mapa_clusters_clustgeo_fungi.html", winslash = "/")
#Exportar a PNG con webshot2 (forzando el paquete)
webshot2::webshot(
  ruta_html,
  "mapa_clustgeo_fungi.png",
  vwidth = 1600,    # ancho más grande
  vheight = 1200,   # alto más grande
  delay = 15        # espera a que todo se cargue
)

# 9. Comparación con K-means
vec_kmeans   <- locality_coords_kmeans$Cluster[match(localidades$Localities, locality_coords_kmeans$Localities)]
vec_clustgeo <- localidades$Cluster
cat("Índice de Rand ajustado K-means vs ClustGeo:", adjustedRandIndex(vec_kmeans, vec_clustgeo), "\n")

# 10. Tabla de contingencia
tabla_confusion <- table(Kmeans = vec_kmeans, ClustGeo = vec_clustgeo)
tabla_df <- as.data.frame.matrix(tabla_confusion)
tabla_df <- cbind(Kmeans = rownames(tabla_df), tabla_df)
write.csv(tabla_df, "comparacion_clusters_tabla_Kmeans_vs_ClustGeo_metazoa.csv", row.names = FALSE)


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

# 1. Crear la tabla de contingencia
tabla_confusion <- table(Kmeans = vec_kmeans, ClustGeo = vec_clustgeo)

# 2. Convertir a data.frame
tabla_df <- as.data.frame.matrix(tabla_confusion)

# 3. Añadir nombres de fila como columna si lo deseas
tabla_df <- cbind(Kmeans = rownames(tabla_df), tabla_df)
tabla_df
# 4. Guardar como CSV
write.csv(tabla_df, "clusters_Kmeans_vs_ClustGeo_metazoa.csv", row.names = FALSE)


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
saveWidget(sankey_plot, file = "sankey_clusters.html", selfcontained = TRUE)
# convertir a PNG
webshot("sankey_clusters.html", "sankey_clusters.png", vwidth = 1200, vheight = 900)
# convertir a PDF
webshot("sankey_clusters.html", "sankey_clusters.pdf", vwidth = 1200, vheight = 900)


### --- 2. Métricas de evaluación ---
# Asegurar que no hay NA y que están alineados
valid_idx <- which(!is.na(vec_kmeans) & !is.na(vec_clustgeo))
vec_kmeans_clean <- vec_kmeans[valid_idx]
vec_clustgeo_clean <- vec_clustgeo[valid_idx]
nmi <- NMI(vec_kmeans_clean, vec_clustgeo_clean)
cat("Normalized Mutual Information (NMI):", round(nmi, 3), "\n")

# Definir los valores que ya calculaste
nmi_value <- round(nmi, 3)
ari_value <- round(adjustedRandIndex(vec_kmeans_clean, vec_clustgeo_clean), 3)

# Nombre de la especie o grupo
nombre_especie <- "Metazoa"

# Crear dataframe con métricas
metricas_comparacion <- data.frame(
  especie = nombre_especie,
  NMI = nmi_value,
  ARI = ari_value
)

# Guardar como CSV (puedes usar append = TRUE si haces un loop por especies)
write.csv(metricas_comparacion, "metrics_nmi_ari_metazoa.csv", row.names = FALSE)

#____________________________________________________________________________________________