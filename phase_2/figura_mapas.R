# =========================================================
# Mapa con pies por zona + sombreado en el MAR (Voronoi)
# =========================================================

# Paquetes (descomenta si necesitas instalar)
# install.packages(c("sf","ggplot2","dplyr","tidyr","rnaturalearth","rnaturalearthdata","scatterpie","RColorBrewer"))

suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(scatterpie)
  library(RColorBrewer)
})

# ---------------- 1) Cargar datos ----------------
setwd("C:/TEMP/FASE1/3_magia_nombre_carpeta/eKOI_metabarcoding_database_FILOS/_RESULTADOS_HAPLOTIPOS_SCRIPT_7")
file_csv <- "aligned_coordenadas_kmeans_protists_fungi.csv"  # ajuta si está en otra ruta
df <- read.csv(file_csv, stringsAsFactors = FALSE)

# Niveles de clúster (aseguramos orden estable)
levs <- sort(unique(as.character(df$Cluster)))

# 🎨 Colores manuales para cada clúster (ajusta a gusto)
# Ejemplo con 10 clústeres
pal <- c(
  "1"  = "#1f77b4",   # azul
  "2"  = "#ff7f0e",   # naranja
  "3"  = "#2ca02c",   # verde
  "4"  = "#d62728",   # rojo
  "5"  = "#9467bd",   # morado
  "6"  = "#8c564b",   # marrón
  "7"  = "#e377c2",   # rosa
  "8"  = "#7f7f7f",   # gris
  "9"  = "#17becf",    # cian
  "10" = "#bcbd22"   # oliva
)

# asegurar columnas mínimas
stopifnot(all(c("lon","lat","Cluster") %in% names(df)))

# Niveles/leyenda para clústeres (consistentes en todo el plot)
levs <- sort(unique(as.character(df$Cluster)))
pal  <- brewer.pal(max(3, length(levs)), "Set3")[seq_along(levs)]
names(pal) <- levs

# ---------------- 2) Puntos sf ----------------
pts_wgs <- st_as_sf(df, coords = c("lon","lat"), crs = 4326, remove = FALSE)

# ---------------- 3) Máscara Península y bbox ----------------
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
iberia_wgs <- world %>% filter(admin %in% c("Spain","Portugal")) %>% st_union() %>% st_as_sf()

# BBox de trabajo (solo Península y alrededores para mar)
bbox_wgs <- st_as_sfc(st_bbox(c(xmin = -10, ymin = 35.5, xmax = 4, ymax = 44.5), crs = 4326))

iberia_wgs_clip <- st_intersection(iberia_wgs, bbox_wgs)
# ---------------- 4) Proyección a métrico ----------------
crs_planar <- 3035  # ETRS89 / LAEA Europe
pts     <- st_transform(pts_wgs,   crs_planar)
iberia  <- st_transform(iberia_wgs_clip, crs_planar)
bbox_pr <- st_transform(bbox_wgs,   crs_planar)

# ---------------- 5) Voronoi y recorte al MAR ----------------
# Voronoi de puntos
voro      <- st_voronoi(st_union(pts))
voro_sf   <- st_collection_extract(voro, "POLYGON") %>% st_as_sf()
st_crs(voro_sf) <- st_crs(pts)  # reasignar CRS

# heredar cluster del punto más cercano
idx <- st_nearest_feature(voro_sf, pts)
voro_sf <- voro_sf %>%
  mutate(Cluster = as.character(pts$Cluster[idx])) %>%
  mutate(Cluster = factor(Cluster, levels = levs))

# máscara de MAR: bbox - tierra
# buffer(0) para evitar problemas topológicos
sea_mask <- st_difference(bbox_pr, st_buffer(iberia, 0))

# recorte Voronoi al mar
voro_sea <- suppressWarnings(st_intersection(voro_sf, sea_mask))

# ---------------- 6) Hex grid y proporciones por celda ----------------
# Malla hexagonal (tamaño ~ 60 km; ajusta 'cellsize' si quieres más/menos detalle)
cellsize <- 60000
hex <- st_make_grid(bbox_pr, cellsize = cellsize, what = "polygons", square = FALSE)
hex <- st_as_sf(hex)
st_crs(hex) <- st_crs(bbox_pr)
hex <- st_intersection(hex, bbox_pr) %>% mutate(id = row_number())

# centroides para ubicar las tartas
hex_c <- st_centroid(hex)

# asignar puntos a celdas
join_idx <- st_join(pts, hex, join = st_within)

# tabla de recuentos por celda y clúster
tab <- join_idx %>%
  st_drop_geometry() %>%
  filter(!is.na(id)) %>%
  count(id, Cluster = as.character(Cluster), name = "n")

# totales por celda
totals <- tab %>% group_by(id) %>% summarise(total = sum(n), .groups = "drop")

# proporciones por celda y clúster en ancho
tab_w <- tab %>%
  group_by(id, Cluster) %>%
  summarise(n = sum(n), .groups = "drop_last") %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  select(id, Cluster, prop) %>%
  pivot_wider(id_cols = id, names_from = Cluster, values_from = prop, values_fill = 0)

# columnas para las porciones (aseguramos el mismo orden que 'levs')
# si faltan clusters en algunas celdas, pivot_wider ya mete 0
pie_cols <- setdiff(names(tab_w), "id")
# (opcional) reordenar columnas a 'levs' si estaban como números/strings
pie_cols <- unique(c(levs, pie_cols))
pie_cols <- pie_cols[pie_cols %in% names(tab_w)]

# unir a centroides y totales
hex_dat <- hex_c %>%
  left_join(tab_w, by = "id") %>%
  left_join(totals, by = "id")

# coordenadas x/y y dataframe plano
coords <- st_coordinates(st_geometry(hex_dat))
hex_df <- hex_dat %>%
  st_drop_geometry() %>%
  mutate(x = coords[,1], y = coords[,2])

# asegurar numérico y sin NA
hex_df <- hex_df %>%
  mutate(across(all_of(pie_cols), ~ replace_na(as.numeric(.), 0))) %>%
  mutate(total = replace_na(total, 0))

# radio de las tartas (ajusta estos números a gusto)
hex_df <- hex_df %>%
  mutate(radius = 15000 + sqrt(total) * 1500)

# solo celdas con puntos
hex_df_plot <- dplyr::filter(hex_df, total > 0)

# ---------------- 7) Gráfico final p_pies ----------------
lim <- st_bbox(bbox_pr)

p_pies <- ggplot() +
  geom_sf(data = voro_sea, aes(fill = Cluster), alpha = 0.18, color = NA) +
  geom_sf(data = iberia,   fill = "white",       color = "grey30", linewidth = 0.4) +
  geom_scatterpie(
    data = hex_df_plot,
    aes(x = x, y = y, r = radius),
    cols = pie_cols, color = NA
  ) +
  scale_fill_manual(values = pal, name = "Cluster", drop = FALSE) +
  coord_sf(xlim = c(lim["xmin"], lim["xmax"]),
           ylim = c(lim["ymin"], lim["ymax"]),
           expand = FALSE) +                                 # ⬅️ fija el encuadre
  labs(
    title = "Proporción de clusters por zona (hex grid)",
    subtitle = "Sombreado en el mar (Voronoi) + pies por celda",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major = element_line(color = "grey85", linewidth = 0.2))

print(p_pies)

# (opcional) guardar
# ggsave("mapa_pies_metazoa.png", p_pies, width = 9, height = 8, dpi = 300)

# Guardar en PDF
ggsave(
  filename = "mapa_pies_protists_fungi2.pdf",  # nombre del archivo de salida
  plot = p_pies,                   # el objeto del gráfico
  device = cairo_pdf,              # mejor soporte de fuentes/vectorial
  width = 8, height = 6, units = "in"  # tamaño en pulgadas
)

