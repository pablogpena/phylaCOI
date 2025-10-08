###======================================================
# HEATMAPS POR VARIABLE PARA MODELO EXPONENCIAL (nls)
###======================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

setwd("C:/TEMP/FASE1/3_magia_nombre_carpeta/eKOI_metabarcoding_database_FILOS/_RESULTADOS FINALES FILOS_script6_2/_maria(finales)/2_heatmap_metric")

df <- read_csv2("comparacion_phylum_dentro_metricas_10otus_modifserv_orden.csv")
df$observed_diff <- as.numeric(df$observed_diff)
df$p_value <- as.numeric(df$p_value)

# Filtrar solo modelo nls
df_nls <- df %>% filter(model == "nls")

# 🔧 Función para forzar todo a mitad superior
prepare_half_matrix <- function(data) {
  
  data_fixed <- data %>%
    rowwise() %>%
    mutate(
      # Ordenamos alfabéticamente los pares (para decidir qué va a x y qué a y)
      phylum_pair = list(sort(c(phylum1, phylum2))),
      # Definimos las nuevas coordenadas
      phylum1_new = phylum_pair[[1]],
      phylum2_new = phylum_pair[[2]],
      # Si cambiamos el orden, invertimos el signo
      observed_diff_new = ifelse(phylum1 == phylum1_new,
                                 observed_diff,
                                 -observed_diff),
      p_value_new = p_value
    ) %>%
    ungroup() %>%
    select(phylum1 = phylum1_new,
           phylum2 = phylum2_new,
           observed_diff = observed_diff_new,
           p_value = p_value_new)
  
  return(data_fixed)
}

plot_heatmap <- function(data, metric_name) {
  
  df_metric <- data %>%
    filter(metric == metric_name) %>%
    prepare_half_matrix() %>%
    mutate(
      signif_label = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE ~ ""
      )
    )
  
  p <- ggplot(df_metric, aes(x = phylum1, y = phylum2, fill = observed_diff)) +
    geom_tile(color = "white") +
    geom_text(aes(label = signif_label), color = "black", size = 3) +
    scale_fill_gradient2(
      low = "red", mid = "white", high = "blue", midpoint = 0
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = paste("Heatmap -", metric_name, "(nls)"),
      x = "Phylum 1", y = "Phylum 2", fill = "Observed Diff",
      caption = "* p < 0.05   ** p < 0.01   *** p < 0.001"
    )
  
  return(p)
}

# ✅ Generar heatmap corregido
p1 <- plot_heatmap(df_nls, "scaled_diversity")
p2 <- plot_heatmap(df_nls, "scaled_log_abundance")
p3 <- plot_heatmap(df_nls, "scaled_log_connections")

print(p1)
print(p2)
print(p3)

# Guardar cada heatmap en un PDF individual
pdf("heatmap_diversity.pdf", width = 8, height = 6)
print(p1)
dev.off()

#pdf("heatmap_abundance.pdf", width = 8, height = 6)
#print(p2)
#dev.off()

#pdf("heatmap_connections.pdf", width = 8, height = 6)
#print(p3)
#dev.off()

###======================================================
# HEATMAPS CON DENDROGRAMA Y ESCALERA (pheatmap)
###======================================================

#### 1. PREPARAR LA MATRIZ CUADRADA A PARTIR DEL CSV

# Cargar librerías
library(dplyr)

# 1. Cargar CSV
setwd("C:/TEMP/FASE1/3_magia_nombre_carpeta/eKOI_metabarcoding_database_FILOS/_RESULTADOS FINALES FILOS_script6_2/_maria(finales)/2_heatmap_metric")
df <- read.csv("comparacion_phylum_dentro_metricas_10otus_abundance_nls.csv", 
               sep = ";", stringsAsFactors = FALSE)

# Organizar datos para compatibilidad
df$phylum1 <- trimws(as.character(df$phylum1))
df$phylum2 <- trimws(as.character(df$phylum2))

# 2. Lista de todos los filos únicos
phylos <- sort(unique(c(df$phylum1, df$phylum2)))

# 3. Crear matrices cuadradas inicializadas en 0
matrix_diff <- matrix(0, nrow = length(phylos), ncol = length(phylos),
                      dimnames = list(phylos, phylos))
matrix_pval <- matrix(0, nrow = length(phylos), ncol = length(phylos),
                      dimnames = list(phylos, phylos))

# 4. Rellenar las matrices
for(i in 1:nrow(df)) {
  p1 <- df$phylum1[i]
  p2 <- df$phylum2[i]
  val_diff <- df$observed_diff[i]
  val_pval <- df$p_value[i]
  
  i1 <- match(p1, phylos)
  i2 <- match(p2, phylos)
  
  matrix_diff[i1, i2] <- val_diff
  matrix_diff[i2, i1] <- -val_diff
  
  matrix_pval[i1, i2] <- val_pval
  matrix_pval[i2, i1] <- val_pval
}

# 5. Guardar resultados
write.csv(matrix_diff, "matriz_abundance_observed_diff.csv", row.names = TRUE)
write.csv(matrix_pval, "matriz_abundance_p_value.csv", row.names = TRUE)


#-------------------------------------------------
#### 2. HEATMAP CON DENDROGRAMA (PAQUETE PHEATMAP)
#-------------------------------------------------

library(pheatmap)
library(RColorBrewer)

# Leer las matrices correctas
observed_diff <- as.matrix(read.csv("matriz_abundance_observed_diff.csv", row.names = 1, check.names = FALSE))
p_values      <- as.matrix(read.csv("matriz_abundance_p_value.csv", row.names = 1, check.names = FALSE))

mode(observed_diff) <- "numeric"

# Crear etiquetas de significación
signif_labels <- matrix("", nrow = nrow(p_values), ncol = ncol(p_values))
signif_labels[p_values < 0.05] <- "*"
signif_labels[p_values < 0.01] <- "**"
signif_labels[p_values < 0.001] <- "***"

rownames(signif_labels) <- rownames(observed_diff)
colnames(signif_labels) <- colnames(observed_diff)

# Escala de colores centrada en 0
lim <- max(abs(observed_diff), na.rm = TRUE)
paleta_div <- colorRampPalette(brewer.pal(11, "RdBu"))(100)   # sin rev()
breaks_div <- seq(-lim, lim, length.out = length(paleta_div) + 1)

# Heatmap con dendrograma
heatmap <- pheatmap(
  observed_diff,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = signif_labels,
  number_color = "black",
  color = paleta_div,
  breaks = breaks_div,
  fontsize_row = 8,
  fontsize_col = 8,
  main = "Heatmap de Observed Diff con dendrograma",
  filename = "prueba_abundance.pdf",   # <<--- clave
  width = 9,
  height = 7
)

heatmap
