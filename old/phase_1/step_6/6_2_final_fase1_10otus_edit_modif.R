# #-------------------
# # ANÁLISIS
# #-------------------
# Cargar y preparar datos
library(readr)
library(dplyr)
library(mgcv)
library(tibble)
library(tidyr)
library(ggplot2)
library(vegan)
library(geosphere)
library(mgcv)
library(broom)
library(purrr)
library(igraph)
library(ggraph)
library(RColorBrewer)

setwd("/workspace/maria")
data <- read_delim("div_abun_conn_10_OTUs_v2.csv", delim = ",", locale = locale(decimal_mark = ".")) 

data <- data %>%
  mutate(
    scaled_log_abundance = (log_abundance - min(log_abundance, na.rm = TRUE)) /
      (max(log_abundance, na.rm = TRUE) - min(log_abundance, na.rm = TRUE)),
    scaled_diversity = (diversity - min(diversity, na.rm = TRUE)) /
      (max(diversity, na.rm = TRUE) - min(diversity, na.rm = TRUE)),
    scaled_log_connections = (log_connections - min(log_connections, na.rm = TRUE)) /
      (max(log_connections, na.rm = TRUE) - min(log_connections, na.rm = TRUE))
  )


ajustar_modelos <- function(data, var_y, nombre_metric) {
  
  # Determinar columna de distancia según variable de interés
  if (grepl("diversity", var_y)) {
    var_dist <- "distance_km_diversity"
  } else if (grepl("abundance", var_y)) {
    var_dist <- "distance_km_abundance"
  } else if (grepl("connections", var_y)) {
    var_dist <- "distance_km_connections"
  } else {
    stop("Variable no reconocida. Usa una que contenga 'diversity', 'abundance' o 'connections'")
  }
  
  # Filtrar datos válidos
  data_filtered <- data %>%
    filter(!is.na(.data[[var_y]]), !is.na(.data[[var_dist]]))
  
  otu_weights <- data_filtered %>%
    group_by(phylum, OTU) %>%
    summarise(weight = 1 / n(), .groups = "drop")
  
  data_filtered_weighted <- data_filtered %>%
    left_join(otu_weights, by = c("phylum", "OTU"))
  
  phyla <- unique(data_filtered_weighted$phylum)
  
  resultados_modelos <- list()
  all_predictions <- list()
  
  for (p in phyla) {
    datos_phylum <- data_filtered_weighted %>% filter(phylum == p)
    
    # Fórmulas para modelos
    f_lm <- as.formula(paste0(var_y, " ~ ", var_dist))
    
    modelo_lm <- lm(f_lm, data = datos_phylum, weights = weight)
    pendiente <- coef(modelo_lm)[2]
    r2_lm <- summary(modelo_lm)$r.squared
    aic_lm <- AIC(modelo_lm)
    bic_lm <- BIC(modelo_lm)
    
    cor_spearman <- cor(datos_phylum[[var_y]], datos_phylum[[var_dist]], method = "spearman", use = "complete.obs")
    
    # EXPONENCIAL
    f_nls <- as.formula(paste0(var_y, " ~ a * exp(-b * ", var_dist, ")"))
    modelo_nls <- tryCatch({
      datos_validos <- datos_phylum %>%
        filter(
          is.finite(.data[[var_y]]),
          !is.na(.data[[var_y]]),
          .data[[var_dist]] > 0,
          is.finite(.data[[var_dist]])
        )
      if (nrow(datos_validos) < 5) stop("Muy pocos datos válidos para NLS")
      
      nls(f_nls,
          data = datos_validos,
          start = list(
            a = max(datos_validos[[var_y]], na.rm = TRUE),
            b = 1 / mean(datos_validos[[var_dist]], na.rm = TRUE)
          ),
          control = nls.control(maxiter = 200))
    }, error = function(e) {
      message("Error en modelo NLS para phylum ", p, ": ", e$message)
      NULL
    })
    
    nls_a <- if (!is.null(modelo_nls)) coef(modelo_nls)["a"] else NA
    nls_b <- if (!is.null(modelo_nls)) coef(modelo_nls)["b"] else NA
    aic_nls <- if (!is.null(modelo_nls)) AIC(modelo_nls) else NA
    bic_nls <- if (!is.null(modelo_nls)) BIC(modelo_nls) else NA
    
    # POWER
    epsilon <- 1e-6
    datos_power <- datos_phylum %>%
      mutate(
        log_dist = log(.data[[var_dist]] + epsilon),
        log_response = log(.data[[var_y]] + epsilon)
      ) %>%
      filter(is.finite(log_dist), is.finite(log_response))
    
    modelo_power_loglog <- tryCatch({
      lm(log_response ~ log_dist, data = datos_power)
    }, error = function(e) NULL)
    
    power_a <- if (!is.null(modelo_power_loglog)) exp(coef(modelo_power_loglog)[1]) else NA
    power_b <- if (!is.null(modelo_power_loglog)) coef(modelo_power_loglog)[2] else NA
    aic_power <- if (!is.null(modelo_power_loglog)) AIC(modelo_power_loglog) else NA
    bic_power <- if (!is.null(modelo_power_loglog)) BIC(modelo_power_loglog) else NA
    
    # GOMPERTZ
    f_gomp <- as.formula(paste0(var_y, " ~ a * exp(-b * exp(-c * ", var_dist, "))"))
    modelo_gompertz <- tryCatch({
      datos_validos <- datos_phylum %>%
        filter(
          is.finite(.data[[var_y]]),
          !is.na(.data[[var_y]]),
          .data[[var_dist]] > 0,
          is.finite(.data[[var_dist]])
        )
      if (nrow(datos_validos) < 5) stop("Muy pocos datos válidos para Gompertz")
      
      nls(f_gomp,
          data = datos_validos,
          start = list(
            a = max(datos_validos[[var_y]], na.rm = TRUE),
            b = 1,
            c = 0.01
          ),
          control = nls.control(maxiter = 200))
    }, error = function(e) {
      message("Error en modelo Gompertz para phylum ", p, ": ", e$message)
      NULL
    })
    
    gompertz_a <- if (!is.null(modelo_gompertz)) coef(modelo_gompertz)["a"] else NA
    gompertz_b <- if (!is.null(modelo_gompertz)) coef(modelo_gompertz)["b"] else NA
    gompertz_c <- if (!is.null(modelo_gompertz)) coef(modelo_gompertz)["c"] else NA
    aic_gompertz <- if (!is.null(modelo_gompertz)) AIC(modelo_gompertz) else NA
    bic_gompertz <- if (!is.null(modelo_gompertz)) BIC(modelo_gompertz) else NA
    
    # GAM
    f_gam <- as.formula(paste0(var_y, " ~ s(", var_dist, ")"))
    modelo_gam <- tryCatch({
      gam(f_gam, data = datos_phylum, weights = weight)
    }, error = function(e) NULL)
    
    gam_r2 <- if (!is.null(modelo_gam)) summary(modelo_gam)$r.sq else NA
    aic_gam <- if (!is.null(modelo_gam)) AIC(modelo_gam) else NA
    bic_gam <- if (!is.null(modelo_gam)) BIC(modelo_gam) else NA
    
    # PREDICCIONES
    dist_pred <- seq(min(datos_phylum[[var_dist]]), max(datos_phylum[[var_dist]]), length.out = 100)
    df_pred <- tibble(!!var_dist := dist_pred)
    df_pred$linear <- predict(modelo_lm, newdata = df_pred)
    df_pred$nls <- if (!is.null(modelo_nls)) predict(modelo_nls, newdata = df_pred) else NA
    df_pred$gam <- if (!is.null(modelo_gam)) predict(modelo_gam, newdata = df_pred) else NA
    df_pred$power <- if (!is.null(modelo_power_loglog)) {
      df_pred$log_dist <- log(df_pred[[var_dist]] + epsilon)
      exp(predict(modelo_power_loglog, newdata = df_pred))
    } else { NA }
    df_pred$gompertz <- if (!is.null(modelo_gompertz)) predict(modelo_gompertz, newdata = df_pred) else NA
    df_pred$phylum <- p
    df_pred$metric <- nombre_metric
    all_predictions[[p]] <- df_pred
    
    # Mantel test usando matriz de distancia ya presente
    mantel_r <- NA
    mantel_p <- NA
    try({
      # Matriz de distancia geográfica
      geo_dist_matrix <- datos_phylum %>%
        select(site_id = cluster, dist = !!sym(var_dist)) %>%
        distinct(site_id, dist)
      
      # Revisamos que haya varios sitios únicos
      site_data <- datos_phylum %>%
        group_by(cluster) %>%
        summarise(response = mean(.data[[var_y]], na.rm = TRUE), .groups = "drop")
      
      if (nrow(site_data) > 3) {
        # Matrices de disimilitud
        dist_response <- dist(site_data$response)
        
        # Para geo_dist usamos distancias entre pares de sitios
        coords <- datos_phylum %>%
          select(cluster, !!sym(var_dist)) %>%
          distinct(cluster, .keep_all = TRUE) %>%
          arrange(cluster)
        
        mat_geo <- as.matrix(dist(coords[[2]]))
        rownames(mat_geo) <- coords$cluster
        colnames(mat_geo) <- coords$cluster
        
        mat_geo <- as.dist(mat_geo)
        
        if (length(mat_geo) == length(dist_response)) {
          mantel_res <- vegan::mantel(mat_geo, dist_response, method = "spearman", permutations = 999)
          mantel_r <- mantel_res$statistic
          mantel_p <- mantel_res$signif
        }
      }
    }, silent = TRUE)
    
    
    # GUARDAR MÉTRICAS
    resultados_modelos[[p]] <- tibble(
      phylum = p,
      metric = nombre_metric,
      valor_max = max(datos_phylum[[var_y]], na.rm = TRUE),
      valor_media = mean(datos_phylum[[var_y]], na.rm = TRUE),
      distancia_max = max(datos_phylum[[var_dist]], na.rm = TRUE),
      pendiente_lineal = pendiente,
      r2_lineal = r2_lm,
      cor_spearman = cor_spearman,
      nls_a = nls_a,
      nls_b = nls_b,
      power_a = power_a,
      power_b = power_b,
      gompertz_a = gompertz_a,
      gompertz_b = gompertz_b,
      gompertz_c = gompertz_c,
      gam_r2 = gam_r2,
      aic_lm = aic_lm, bic_lm = bic_lm,
      aic_nls = aic_nls, bic_nls = bic_nls,
      aic_power = aic_power, bic_power = bic_power,
      aic_gompertz = aic_gompertz, bic_gompertz = bic_gompertz,
      aic_gam = aic_gam, bic_gam = bic_gam,
      mantel_r = mantel_r,
      mantel_p = mantel_p
    )
  }
  
  list(
    resultados = bind_rows(resultados_modelos),
    predicciones = bind_rows(all_predictions)
  )
}

res_abundancia <- ajustar_modelos(data, var_y = "scaled_log_abundance", nombre_metric = "log_abundance")
res_diversidad <- ajustar_modelos(data, var_y = "scaled_diversity", nombre_metric = "diversity")
res_conexiones <- ajustar_modelos(data, var_y = "scaled_log_connections", nombre_metric = "log_connections")

# Unir resultados
resultados_modelos_totales <- dplyr::bind_rows(
  res_abundancia$resultados,
  res_diversidad$resultados,
  res_conexiones$resultados
)
predicciones_totales <- dplyr::bind_rows(
  res_abundancia$predicciones,
  res_diversidad$predicciones,
  res_conexiones$predicciones
)

predicciones_largo <- predicciones_totales %>%
  tidyr::pivot_longer(
    cols = c("linear", "nls", "gam", "power", "gompertz"),
    names_to = "model",
    values_to = "predicted"
  )

predicciones_largo <- predicciones_largo %>%
  dplyr::mutate(
    distance = dplyr::case_when(
      metric == "diversity" ~ distance_km_diversity,
      metric == "log_abundance" ~ distance_km_abundance,
      metric == "log_connections" ~ distance_km_connections,
      TRUE ~ NA_real_
    )
  )


# +
write.csv(resultados_modelos_totales, "resultados_modelos_combinados_10otus_modif.csv", row.names = FALSE)
write.csv(predicciones_largo, "predicciones_modelos_combinadas_10otus_modif.csv", row.names = FALSE)

# +
#PLOT: Model Comparison by Phylum and Metric"
model_comparision <- ggplot(predicciones_largo, aes(x = distance, y = predicted, color = model)) +
  geom_line() +
  facet_grid(phylum ~ metric) +
  labs(
    title = "Model Comparison by Phylum and Metric",
    y = "Predicted Value (scaled)",
    x = "Distance (km)"
  ) +
  theme_minimal()

ggsave("model_comparision_all_modif.pdf", plot = model_comparision, width = 8, height = 5)
model_comparision
# +
#PLOT: version solo 200km
model_comparision_200km <- ggplot(predicciones_largo, aes(x = distance, y = predicted, color = model)) +
  geom_line() +
  facet_grid(phylum ~ metric) +
  labs(
    title = "Model Comparison by Phylum and Metric",
    y = "Predicted Value (scaled)",
    x = "Distance (km)"
  ) +
  xlim(0, 202) +   # Limita el eje x entre 0 y 200 km para poder ver bien las correlaciones entre los distintos modelos
  theme_minimal()
model_comparision_200km

ggsave("model_comparision_200km_modif.pdf", plot = model_comparision_200km, width = 8, height = 5)

# +
#PLOT: lm, power, gam solo 200km
#predicciones_filtradas <- predicciones_largo %>%
#  filter(model %in% c("linear", "gam", "nls"))
#model_comparision_lm_power_gam_nls <- ggplot(predicciones_filtradas, aes(x = distance, y = predicted, color = model)) +
#  geom_line() +
#  facet_grid(phylum ~ metric) +
#  labs(
#    title = "Model Comparison by Phylum and Metric",
#    y = "Predicted Value (scaled)",
#    x = "Distance (km)"
#  ) +
#  xlim(0, 200) +   # si quieres limitar a 200 km como antes
#  theme_minimal()

#ggsave("model_comparision_lm_nls_gam.pdf", plot = model_comparision_lm_power_gam_nls, width = 8, height = 5)
# -

#-----------------------------------------------------
#####PARTE 2: CORRELACION SPEARMAN
#-----------------------------------------------------

######correlacion entre abundancia y diversidad nucleotidica
# Variables a comparar
var_pairs <- list(
  c("log_abundance", "diversity"),
  c("log_connections", "diversity"),
  c("log_abundance", "log_connections")
)

# ---- Correlación GLOBAL ----
cor_global_df <- purrr::map_dfr(var_pairs, function(vars) {
  var1 <- vars[1]
  var2 <- vars[2]
  
  tmp <- data %>%
    dplyr::filter(!is.na(.data[[var1]]) & !is.na(.data[[var2]]))
  
  test <- cor.test(tmp[[var1]], tmp[[var2]], method = "spearman")
  
  tibble::tibble(
    phylum = "Global",
    variable_x = var1,
    variable_y = var2,
    rho = test$estimate,
    p_value = test$p.value,
    n = nrow(tmp)
  )
})

# ---- Correlación por PHYLUM ----
cor_phylum_df <- purrr::map_dfr(var_pairs, function(vars) {
  var1 <- vars[1]
  var2 <- vars[2]
  
  data %>%
    dplyr::filter(!is.na(.data[[var1]]) & !is.na(.data[[var2]])) %>%
    dplyr::group_by(phylum) %>%
    dplyr::summarise(
      variable_x = var1,
      variable_y = var2,
      rho = cor(.data[[var1]], .data[[var2]], method = "spearman"),
      p_value = cor.test(.data[[var1]], .data[[var2]], method = "spearman")$p.value,
      n = dplyr::n(),
      .groups = "drop"
    )
})

# ---- Unir todo ----
cor_total <- dplyr::bind_rows(cor_global_df, cor_phylum_df)

# ---- Guardar ----
readr::write_csv(cor_total, "correlaciones_multivariadas.csv")

comparisons <- list(
  c("log_abundance", "diversity"),
  c("log_connections", "diversity"),
  c("log_abundance", "log_connections")
)

#PLOT: Spearman correlation
plot_list <- purrr::map(comparisons, function(pair) {
  var1 <- pair[1]
  var2 <- pair[2]
  
  cor_phylum_df %>%
    filter(variable_x == var1, variable_y == var2) %>%
    ggplot(aes(x = reorder(phylum, -rho), y = rho, fill = p_value < 0.05)) +
    geom_col() +
    geom_text(aes(label = paste0("ρ = ", round(rho, 2))), vjust = -0.5, size = 3) +
    scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray70"),
                      name = "p < 0.05") +
    labs(
      title = paste("Spearman Correlation:", var1, "vs", var2, "by phylum"),
      x = "Phylum", y = "ρ (Spearman)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
})

# +
print(plot_list[[1]])
print(plot_list[[2]])
print(plot_list[[3]])

# Guardar cada plot en un archivo PDF separado
ggsave("plot_1.pdf", plot = plot_list[[1]], device = "pdf", width = 8, height = 6)
ggsave("plot_2.pdf", plot = plot_list[[2]], device = "pdf", width = 8, height = 6)
ggsave("plot_3.pdf", plot = plot_list[[3]], device = "pdf", width = 8, height = 6)



#PLOT: Heatmap of Correlations by Phylum
heatmap_plot <- ggplot(cor_total %>% filter(phylum != "Global"),
                       aes(x = variable_x, y = variable_y, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0("ρ=", round(rho, 2), 
                               ifelse(p_value < 0.05, "*", ""))),
            size = 3) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0,
                       name = "Spearman ρ") +
  facet_wrap(~ phylum) +
  labs(title = "Heatmap of Correlations by Phylum",
       x = "Variable X", y = "Variable Y") +
  theme_minimal()

ggsave("heatmap_correlations.pdf", plot = heatmap_plot,
       device = "pdf", width = 10, height = 8)
# -


# Preparar data global
edges <- cor_total %>% 
  filter(phylum == "Global") %>%
  select(from = variable_x, to = variable_y, rho)

# +
g <- graph_from_data_frame(edges, directed = FALSE)

#PLOT: Global Correlation Network (Spearman ρ)
network_plot <- ggraph(g, layout = "circle") +
  geom_edge_link(aes(width = abs(rho), color = rho), alpha = 0.8) +
  geom_node_point(size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4) +
  scale_edge_color_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0) +
  theme_void() +
  labs(title = "Global Correlation Network (Spearman ρ)")

ggsave("network_correlation.pdf", plot = network_plot,
       device = "pdf", width = 8, height = 8)



# =======================================================
# PARTE 3: Calcular pendientes y comparaciones por modelo
# =======================================================

library(tidyverse)
library(mgcv)

# --- operador auxiliar ---
`%||%` <- function(a, b) if (!is.null(a)) a else b

# --- FUNCIÓN BASE: get_slope (con 5 modelos) ---
get_slope <- function(df, yvar, model_type = "lm", distance_col,
                      x_ref = NULL, gam_k = NULL,
                      eps_log = 1e-8, h_frac = 50) {
  y <- df[[yvar]]
  x <- df[[distance_col]]
  ok <- is.finite(y) & is.finite(x)
  y <- y[ok]; x <- x[ok]
  n <- length(y)
  if (n < 5) return(NA_real_)
  
  if (is.null(x_ref)) x_ref <- stats::median(x, na.rm = TRUE)
  h <- stats::sd(x, na.rm = TRUE) / h_frac
  if (!is.finite(h) || h == 0) h <- (max(x) - min(x)) / (2*h_frac)
  
  y_pos <- ifelse(y <= 0, eps_log, y)
  x_pos <- ifelse(x <= 0, eps_log, x)
  
  tryCatch({
    if (model_type == "lm") {
      fit <- lm(y ~ x)
      unname(coef(fit)[["x"]])
      
    } else if (model_type == "power") {
      fit <- lm(log(y_pos) ~ log(x_pos))
      unname(coef(fit)[["log(x_pos)"]])
      
    } else if (model_type == "gam") {
      k_use <- gam_k %||% max(5, min(10, floor(n/3)))
      fit <- mgcv::gam(y ~ s(x, k = k_use), method = "REML")
      pred <- function(xx) predict(fit, newdata = data.frame(x = xx), type = "response")
      (pred(x_ref + h) - pred(x_ref - h)) / (2*h)
      
    } else if (model_type == "nls") {
      if (length(unique(x)) < 3) return(NA_real_)
      keep <- y > 0 & is.finite(y) & is.finite(x)
      y_raw <- y[keep]; x_raw <- x[keep]
      if (length(y_raw) < 5) return(NA_real_)
      
      mu_x <- mean(x_raw); sd_x <- sd(x_raw)
      if (!is.finite(sd_x) || sd_x == 0) return(NA_real_)
      
      z <- (x_raw - mu_x)/sd_x
      z_ref <- (x_ref - mu_x)/sd_x
      df_nls <- data.frame(z = z, y = y_raw)
      
      fit0 <- lm(log(y) ~ z, data = df_nls)
      b0 <- unname(coef(fit0)[["z"]])
      a0 <- unname(exp(coef(fit0)[["(Intercept)"]]))
      
      fit <- nls(
        y ~ a * exp(b * z),
        data = df_nls,
        start = list(a = a0, b = b0),
        algorithm = "port",
        lower = c(a = .Machine$double.eps, b = -50),
        upper = c(a = Inf, b = 50),
        control = nls.control(warnOnly = TRUE, maxiter = 500)
      )
      
      co <- coef(fit)
      slope <- (co[["a"]] * co[["b"]] * exp(co[["b"]] * z_ref)) * (1/sd_x)
      unname(as.numeric(slope))
      
    } else if (model_type == "gompertz") {
      a_start <- max(y_pos, na.rm = TRUE)
      c_start <- 0.05
      b_start <- 1
      df_g <- data.frame(x = x, y = y_pos)
      fit <- nls(y ~ a * exp(-b * exp(-c * x)),
                 data = df_g,
                 start = list(a = a_start, b = b_start, c = c_start),
                 control = nls.control(warnOnly = TRUE, maxiter = 500))
      co <- coef(fit)
      yx <- co[["a"]] * exp(-co[["b"]] * exp(-co[["c"]] * x_ref))
      unname(yx * co[["b"]] * co[["c"]] * exp(-co[["c"]] * x_ref))
      
    } else {
      stop("model_type no reconocido: ", model_type)
    }
  }, error = function(e) {
    warning(sprintf("get_slope(%s) falló: %s", model_type, e$message))
    NA_real_
  })
}

# --- FUNCIÓN: permutation_test_slope ---
permutation_test_slope <- function(data, phylum, model_type, 
                                   metric1 = "scaled_log_abundance", 
                                   metric2 = "scaled_diversity", 
                                   B = 999) {
  df <- data %>% 
    filter(phylum == !!phylum) %>%
    filter(!is.na(.data[[metric1]]), !is.na(.data[[metric2]]))
  
  n <- nrow(df)
  if (n < 5) return(tibble(phylum = phylum, model = model_type, observed_diff = NA, p_value = NA, n = n))
  
  dist_for <- function(m) case_when(
    m == "scaled_log_abundance"   ~ "distance_km_abundance",
    m == "scaled_diversity"       ~ "distance_km_diversity",
    m == "scaled_log_connections" ~ "distance_km_connections",
    TRUE ~ NA_character_
  )
  dist1 <- dist_for(metric1)
  dist2 <- dist_for(metric2)
  
  slope1 <- get_slope(df, metric1, model_type, dist1)
  slope2 <- get_slope(df, metric2, model_type, dist2)
  observed_diff <- as.numeric(slope1 - slope2)
  
  diffs <- numeric(B)
  for (b in seq_len(B)) {
    swap <- sample(c(TRUE, FALSE), n, replace = TRUE)
    temp1 <- ifelse(swap, df[[metric2]], df[[metric1]])
    temp2 <- ifelse(swap, df[[metric1]], df[[metric2]])
    dist_temp1 <- ifelse(swap, df[[dist2]], df[[dist1]])
    dist_temp2 <- ifelse(swap, df[[dist1]], df[[dist2]])
    
    df_b1 <- tibble(metric = temp1, distance_km = dist_temp1)
    df_b2 <- tibble(metric = temp2, distance_km = dist_temp2)
    
    slope1_b <- get_slope(df_b1, "metric", model_type, "distance_km")
    slope2_b <- get_slope(df_b2, "metric", model_type, "distance_km")
    
    diffs[b] <- slope1_b - slope2_b
  }
  
  valid_diffs <- diffs[is.finite(diffs)]
  p_val <- if (!is.finite(observed_diff) || length(valid_diffs) == 0) NA_real_
  else mean(abs(valid_diffs) >= abs(observed_diff))
  
  tibble(
    phylum = phylum,
    model = model_type,
    observed_diff = observed_diff,
    p_value = p_val,
    n = n
  )
}

# --- COMBINACIONES DE MÉTRICAS ---
metric_pairs <- list(
  c("scaled_log_abundance", "scaled_diversity"),
  c("scaled_log_abundance", "scaled_log_connections"),
  c("scaled_diversity", "scaled_log_connections")
)

# --- PERMUTACIONES PARA TODOS LOS MODELOS ---
resultados_perm_metricas <- purrr::pmap_dfr(
  tidyr::expand_grid(
    phylum = unique(data$phylum),
    model_type = c("lm", "power", "gam", "nls", "gompertz"),
    pair = metric_pairs
  ),
  function(phylum, model_type, pair) {
    cat("Analizando phylum:", phylum,
        "| Modelo:", model_type,
        "| Métrica 1:", pair[[1]],
        "| Métrica 2:", pair[[2]], "\n")
    
    permutation_test_slope(
      data,
      phylum = phylum,
      model_type = model_type,
      metric1 = pair[[1]],
      metric2 = pair[[2]],
      B = 999
    ) %>%
      mutate(metric1 = pair[[1]], metric2 = pair[[2]])
  }
)

# --- GUARDAR RESULTADOS DE PERMUTACIONES ---
write_csv(resultados_perm_metricas, "comparaciones_metricas_permutaciones_10otus_modif.csv")

# --- FUNCIÓN: obtener pendientes por modelo ---
get_all_slopes <- function(data, model_type = "lm") {
  phyla <- unique(data$phylum)
  metrics <- c("scaled_log_abundance", "scaled_diversity", "scaled_log_connections")
  
  combinaciones <- tidyr::expand_grid(phylum = phyla, metric = metrics)
  
  purrr::pmap_dfr(
    list(combinaciones$phylum, combinaciones$metric),
    function(ph, met) {
      df_sub <- dplyr::filter(data, phylum == ph)
      
      dist_col <- dplyr::case_when(
        met == "scaled_log_abundance" ~ "distance_km_abundance",
        met == "scaled_diversity" ~ "distance_km_diversity",
        met == "scaled_log_connections" ~ "distance_km_connections"
      )
      
      slope_val <- get_slope(df_sub, yvar = met, model_type = model_type, distance_col = dist_col)
      
      # 🔹 Convertir siempre a numérico simple (NA si falla)
      slope_num <- suppressWarnings(as.numeric(slope_val[1]))
      
      tibble(phylum = ph, metric = met, model = model_type, slope = slope_num)
    }
  )
}


# --- PENDIENTES PARA TODOS LOS MODELOS ---
slopes_all <- bind_rows(
  get_all_slopes(data, "lm"),
  get_all_slopes(data, "power"),
  get_all_slopes(data, "gam"),
  get_all_slopes(data, "nls"),
  get_all_slopes(data, "gompertz")
)

# --- GUARDAR PENDIENTES ---
write_csv(slopes_all, "pendientes_metricas_por_modelo_10otus_modif.csv")

# Guardar todo el entorno actual
save.image(file = "mi_entorno_6_2.RData")
# =====================================================
# FIN DEL SCRIPT
# =====================================================