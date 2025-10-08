# =============================
# Comparación de tres grupos
# =============================

# install.packages(c("dplyr","tidyr","ggplot2","aricode","mclust","ggalluvial"))
library(dplyr)
library(tidyr)
library(ggplot2)
library(aricode)
library(mclust)
library(ggalluvial)
library(purrr)
library(patchwork)

setwd("C:/TEMP/FASE1/3_magia_nombre_carpeta/eKOI_metabarcoding_database_FILOS/_RESULTADOS_HAPLOTIPOS_SCRIPT_7")

# ---------- 1) Cargar datos ----------
meta <- read.csv("aligned_coordenadas_kmeans_metazoos.csv", stringsAsFactors = FALSE)  %>%
  select(Localities, Cluster_meta = Cluster)
arq  <- read.csv("aligned_coordenadas_kmeans_arqueaplastidia.csv", stringsAsFactors = FALSE) %>%
  select(Localities, Cluster_arq  = Cluster)
pf   <- read.csv("aligned_coordenadas_kmeans_protists_fungi.csv", stringsAsFactors = FALSE) %>%
  select(Localities, Cluster_pf   = Cluster)

# ---------- 2) Uniones ----------
# Pares
ov_meta_arq <- inner_join(meta, arq, by = "Localities")
ov_meta_pf  <- inner_join(meta, pf,  by = "Localities")
ov_arq_pf   <- inner_join(arq,  pf,  by = "Localities")

# Triple solapamiento
ov_all3 <- reduce(list(meta, arq, pf), inner_join, by = "Localities")

# ---------- 3) Funciones aux ----------
pair_metrics <- function(x, y, nperm = 1000, seed = 123){
  set.seed(seed)
  ari  <- ARI(x, y)
  nmi  <- NMI(x, y)
  perm <- replicate(nperm, ARI(x, sample(y)))
  pval <- mean(perm >= ari)
  list(ARI = ari, NMI = nmi, pval = pval, perm = perm)
}

mk_heatmap <- function(tab, title){
  # normalizar por filas (proporciones), mantener también el conteo
  df <- as.data.frame(tab)
  names(df) <- c("row","col","n")
  row_tot <- df %>% group_by(row) %>% summarise(t = sum(n), .groups="drop")
  df <- df %>% left_join(row_tot, by="row") %>% mutate(prop = ifelse(t>0, n/t, 0))
  ggplot(df, aes(x = col, y = row, fill = prop)) +
    geom_tile() +
    geom_text(aes(label = n), size = 3) +
    scale_fill_viridis_c(labels = scales::percent_format(accuracy = 1)) +
    labs(title = title, x = "Clusters (2º grupo)", y = "Clusters (1º grupo)", fill = "% fila") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ---------- 4) Métricas por par ----------
m_meta_arq <- pair_metrics(ov_meta_arq$Cluster_meta, ov_meta_arq$Cluster_arq)
m_meta_pf  <- pair_metrics(ov_meta_pf$Cluster_meta,   ov_meta_pf$Cluster_pf)
m_arq_pf   <- pair_metrics(ov_arq_pf$Cluster_arq,     ov_arq_pf$Cluster_pf)

metrics_tbl <- tibble(
  pares = c("Metazoos vs Arqueaplastidia", "Metazoos vs Protists/Fungi", "Arqueaplastidia vs Protists/Fungi"),
  ARI   = c(m_meta_arq$ARI, m_meta_pf$ARI, m_arq_pf$ARI),
  NMI   = c(m_meta_arq$NMI, m_meta_pf$NMI, m_arq_pf$NMI),
  pval_perm = c(m_meta_arq$pval, m_meta_pf$pval, m_arq_pf$pval)
)

print(metrics_tbl)

# Exportar la tabla de métricas a CSV
write.csv(metrics_tbl, "metrics_clusters.csv", row.names = FALSE)

# ---------- 5) Heatmaps de tablas (sintético) ----------
# Tablas de contingencia para cada par
tab_meta_arq <- with(ov_meta_arq, table(Cluster_meta, Cluster_arq))
tab_meta_pf  <- with(ov_meta_pf,  table(Cluster_meta, Cluster_pf))
tab_arq_pf   <- with(ov_arq_pf,   table(Cluster_arq,  Cluster_pf))

p_hm_ma <- mk_heatmap(tab_meta_arq, "Metazoos vs Arqueaplastidia")
p_hm_ma
p_hm_mp <- mk_heatmap(tab_meta_pf,  "Metazoos vs Protists/Fungi")
p_hm_mp
p_hm_ap <- mk_heatmap(tab_arq_pf,   "Arqueaplastidia vs Protists/Fungi")
p_hm_ap

# Combinar los tres heatmaps en una sola figura
combined_plot <- p_hm_ma | p_hm_mp | p_hm_ap   # en 3 filas
print(combined_plot)
# también podrías usar p_hm_ma | p_hm_mp | p_hm_ap para columnas

# Exportar a un único PDF
ggsave("heatmaps_comparacion.pdf", plot = combined_plot,
       device = cairo_pdf, width = 8, height = 12, units = "in")

# ---------- 6) Aluvial a 3 vías (solo localidades presentes en los tres) ----------
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)

# alluv_df: una fila por combinación de 3 clústeres con su Freq
# (si no lo tienes aún)
# alluv_df <- ov_all3 %>%
#   count(Cluster_meta, Cluster_arq, Cluster_pf, name = "Freq")
alluv_df <- ov_all3 %>%
  count(Cluster_meta, Cluster_arq, Cluster_pf, name = "Freq") %>%
  mutate(
    Cluster_meta = as.character(Cluster_meta),
    Cluster_arq  = as.character(Cluster_arq),
    Cluster_pf   = as.character(Cluster_pf)
  ) %>%
  arrange(Cluster_meta, Cluster_arq, Cluster_pf)
alluv_df <- alluv_df %>%
  mutate(
    Cluster_meta = as.character(Cluster_meta),
    Cluster_arq  = as.character(Cluster_arq),
    Cluster_pf   = as.character(Cluster_pf)
  )

# 1) id único por combinación + color por cluster de metazoos
df_wide <- alluv_df %>%
  arrange(Cluster_meta, Cluster_arq, Cluster_pf) %>%
  mutate(
    id        = row_number(),
    fill_meta = Cluster_meta
  ) %>%
  select(id, Freq, fill_meta, Cluster_meta, Cluster_arq, Cluster_pf)

# 2) Formato largo (¡una fila por (id, eje)!)
ps <- df_wide %>%
  pivot_longer(
    cols = c(Cluster_meta, Cluster_arq, Cluster_pf),
    names_to  = "x",
    values_to = "stratum"
  )

# 3) Ejes bonitos
ps$x <- factor(ps$x,
               levels = c("Cluster_meta","Cluster_arq","Cluster_pf"),
               labels = c("Metazoos","Arqueaplastidia","Protists/Fungi"))

# 4) Paleta (asegúrate de tenerla; claves deben ser strings de Cluster_meta)
if (!exists("pal")) {
  pal <- c("1"="#1f77b4","2"="#ff7f0e","3"="#2ca02c","4"="#d62728","5"="#9467bd",
           "6"="#8c564b","7"="#e377c2","8"="#7f7f7f","9"="#bcbd22","10"="#17becf")
}
ps$fill_meta <- as.character(ps$fill_meta)

# (Comprobación opcional: un único id por eje)
# stopifnot(!any(duplicated(ps[c("x","id")])))

# 5) Gráfico tipo aluvial/paralelo (estable con ggforce)
p_alluvial_alt <- ggplot(
  ps,
  aes(x = x, id = id, split = stratum, value = Freq, fill = fill_meta)
) +
  geom_parallel_sets(alpha = 0.7, axis.width = 0.15) +
  geom_parallel_sets_axes(axis.width = 0.15, fill = "grey90", color = "grey30") +
  geom_parallel_sets_labels(angle = 0, size = 3, color = "black") +
  scale_fill_manual(values = pal, name = "Cluster (Metazoos)", drop = FALSE) +
  labs(
    title = "Correspondencia de clústeres (localidades presentes en los 3 grupos)",
    x = NULL, y = "Número de localidades"
  ) +
  theme_minimal()

print(p_alluvial_alt)

# Exportar a un único PDF
ggsave("comparacion_clusters.pdf", plot = p_alluvial_alt,
       device = cairo_pdf, width = 8, height = 12, units = "in")


# ---------- 7) (Opcional) Histogramas nulos ARI por par ----------
par(mfrow = c(1,3))
hist(m_meta_arq$perm, breaks=30, main="ARI nulo: Meta vs Arq", xlab="ARI permutado"); abline(v=m_meta_arq$ARI, col="red", lwd=2)
hist(m_meta_pf$perm,  breaks=30, main="ARI nulo: Meta vs PF",  xlab="ARI permutado"); abline(v=m_meta_pf$ARI,  col="red", lwd=2)
hist(m_arq_pf$perm,   breaks=30, main="ARI nulo: Arq vs PF",   xlab="ARI permutado"); abline(v=m_arq_pf$ARI,   col="red", lwd=2)
par(mfrow = c(1,1))
