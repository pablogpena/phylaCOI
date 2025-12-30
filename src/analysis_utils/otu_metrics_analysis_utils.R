# -*- coding: utf-8 -*-
# Utilities for analyzing combined OTU metrics across phyla.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(mgcv)
  library(vegan)
  library(purrr)
  library(igraph)
  library(ggraph)
  library(RColorBrewer)
  library(pheatmap)
})

validate_metrics_columns <- function(df) {
  required <- c(
    "phylum",
    "OTU",
    "cluster",
    "log_abundance",
    "diversity",
    "log_connections",
    "distance_km_abundance",
    "distance_km_diversity",
    "distance_km_connections"
  )
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
}

read_metrics_table <- function(input_file) {
  df <- readr::read_csv(input_file, show_col_types = FALSE)
  validate_metrics_columns(df)
  df
}

scale_to_unit <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }

  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)

  if (!is.finite(x_min) || !is.finite(x_max) || x_max == x_min) {
    return(rep(NA_real_, length(x)))
  }

  (x - x_min) / (x_max - x_min)
}

add_scaled_metrics <- function(df) {
  df %>%
    mutate(
      scaled_log_abundance = scale_to_unit(log_abundance),
      scaled_diversity = scale_to_unit(diversity),
      scaled_log_connections = scale_to_unit(log_connections)
    )
}

distance_column_for <- function(response_col) {
  if (grepl("diversity", response_col)) {
    return("distance_km_diversity")
  }
  if (grepl("abundance", response_col)) {
    return("distance_km_abundance")
  }
  if (grepl("connections", response_col)) {
    return("distance_km_connections")
  }

  stop("Unrecognized response column: ", response_col)
}

fit_models_for_metric <- function(data, response_col, metric_label) {
  distance_col <- distance_column_for(response_col)

  data_filtered <- data %>%
    filter(!is.na(.data[[response_col]]), !is.na(.data[[distance_col]]))

  otu_weights <- data_filtered %>%
    group_by(phylum, OTU) %>%
    summarise(weight = 1 / n(), .groups = "drop")

  data_weighted <- data_filtered %>%
    left_join(otu_weights, by = c("phylum", "OTU"))

  phyla <- unique(data_weighted$phylum)
  model_results <- list()
  predictions <- list()

  for (phylum_name in phyla) {
    df <- data_weighted %>% filter(phylum == phylum_name)

    # Linear model
    lm_formula <- as.formula(paste0(response_col, " ~ ", distance_col))
    model_lm <- lm(lm_formula, data = df, weights = weight)
    lm_slope <- coef(model_lm)[2]
    r2_lm <- summary(model_lm)$r.squared
    aic_lm <- AIC(model_lm)
    bic_lm <- BIC(model_lm)

    spearman_rho <- cor(
      df[[response_col]],
      df[[distance_col]],
      method = "spearman",
      use = "complete.obs"
    )

    # Exponential (NLS)
    nls_formula <- as.formula(paste0(response_col, " ~ a * exp(-b * ", distance_col, ")"))
    model_nls <- tryCatch({
      df_valid <- df %>%
        filter(
          is.finite(.data[[response_col]]),
          !is.na(.data[[response_col]]),
          .data[[distance_col]] > 0,
          is.finite(.data[[distance_col]])
        )
      if (nrow(df_valid) < 5) {
        stop("Too few valid points for NLS")
      }
      nls(
        nls_formula,
        data = df_valid,
        start = list(
          a = max(df_valid[[response_col]], na.rm = TRUE),
          b = 1 / mean(df_valid[[distance_col]], na.rm = TRUE)
        ),
        control = nls.control(maxiter = 200)
      )
    }, error = function(e) NULL)

    nls_a <- if (!is.null(model_nls)) coef(model_nls)["a"] else NA
    nls_b <- if (!is.null(model_nls)) coef(model_nls)["b"] else NA
    aic_nls <- if (!is.null(model_nls)) AIC(model_nls) else NA
    bic_nls <- if (!is.null(model_nls)) BIC(model_nls) else NA

    # Power (log-log)
    epsilon <- 1e-6
    data_power <- df %>%
      mutate(
        log_dist = log(.data[[distance_col]] + epsilon),
        log_response = log(.data[[response_col]] + epsilon)
      ) %>%
      filter(is.finite(log_dist), is.finite(log_response))

    model_power <- tryCatch({
      lm(log_response ~ log_dist, data = data_power)
    }, error = function(e) NULL)

    power_a <- if (!is.null(model_power)) exp(coef(model_power)[1]) else NA
    power_b <- if (!is.null(model_power)) coef(model_power)[2] else NA
    aic_power <- if (!is.null(model_power)) AIC(model_power) else NA
    bic_power <- if (!is.null(model_power)) BIC(model_power) else NA

    # Gompertz
    gompertz_formula <- as.formula(
      paste0(response_col, " ~ a * exp(-b * exp(-c * ", distance_col, "))")
    )
    model_gompertz <- tryCatch({
      df_valid <- df %>%
        filter(
          is.finite(.data[[response_col]]),
          !is.na(.data[[response_col]]),
          .data[[distance_col]] > 0,
          is.finite(.data[[distance_col]])
        )
      if (nrow(df_valid) < 5) {
        stop("Too few valid points for Gompertz")
      }
      nls(
        gompertz_formula,
        data = df_valid,
        start = list(
          a = max(df_valid[[response_col]], na.rm = TRUE),
          b = 1,
          c = 0.01
        ),
        control = nls.control(maxiter = 200)
      )
    }, error = function(e) NULL)

    gompertz_a <- if (!is.null(model_gompertz)) coef(model_gompertz)["a"] else NA
    gompertz_b <- if (!is.null(model_gompertz)) coef(model_gompertz)["b"] else NA
    gompertz_c <- if (!is.null(model_gompertz)) coef(model_gompertz)["c"] else NA
    aic_gompertz <- if (!is.null(model_gompertz)) AIC(model_gompertz) else NA
    bic_gompertz <- if (!is.null(model_gompertz)) BIC(model_gompertz) else NA

    # GAM
    gam_formula <- as.formula(paste0(response_col, " ~ s(", distance_col, ")"))
    model_gam <- tryCatch({
      mgcv::gam(gam_formula, data = df, weights = weight)
    }, error = function(e) NULL)

    gam_r2 <- if (!is.null(model_gam)) summary(model_gam)$r.sq else NA
    aic_gam <- if (!is.null(model_gam)) AIC(model_gam) else NA
    bic_gam <- if (!is.null(model_gam)) BIC(model_gam) else NA

    # Predictions across the distance range
    dist_pred <- seq(
      min(df[[distance_col]]),
      max(df[[distance_col]]),
      length.out = 100
    )
    pred_df <- tibble(!!distance_col := dist_pred)
    pred_df$linear <- predict(model_lm, newdata = pred_df)
    pred_df$nls <- if (!is.null(model_nls)) predict(model_nls, newdata = pred_df) else NA
    pred_df$gam <- if (!is.null(model_gam)) predict(model_gam, newdata = pred_df) else NA
    pred_df$power <- if (!is.null(model_power)) {
      pred_df$log_dist <- log(pred_df[[distance_col]] + epsilon)
      exp(predict(model_power, newdata = pred_df))
    } else {
      NA
    }
    pred_df$gompertz <- if (!is.null(model_gompertz)) {
      predict(model_gompertz, newdata = pred_df)
    } else {
      NA
    }
    pred_df$phylum <- phylum_name
    pred_df$metric <- metric_label
    predictions[[phylum_name]] <- pred_df

    # Mantel test using per-cluster distances from the combined metrics.
    mantel_r <- NA
    mantel_p <- NA
    try({
      site_data <- df %>%
        group_by(cluster) %>%
        summarise(response = mean(.data[[response_col]], na.rm = TRUE), .groups = "drop")

      if (nrow(site_data) > 3) {
        dist_response <- dist(site_data$response)
        coords <- df %>%
          select(cluster, !!sym(distance_col)) %>%
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

    # Collect summary metrics per phylum and metric.
    model_results[[phylum_name]] <- tibble(
      phylum = phylum_name,
      metric = metric_label,
      max_value = max(df[[response_col]], na.rm = TRUE),
      mean_value = mean(df[[response_col]], na.rm = TRUE),
      max_distance = max(df[[distance_col]], na.rm = TRUE),
      lm_slope = lm_slope,
      r2_lm = r2_lm,
      spearman_rho = spearman_rho,
      nls_a = nls_a,
      nls_b = nls_b,
      power_a = power_a,
      power_b = power_b,
      gompertz_a = gompertz_a,
      gompertz_b = gompertz_b,
      gompertz_c = gompertz_c,
      gam_r2 = gam_r2,
      aic_lm = aic_lm,
      bic_lm = bic_lm,
      aic_nls = aic_nls,
      bic_nls = bic_nls,
      aic_power = aic_power,
      bic_power = bic_power,
      aic_gompertz = aic_gompertz,
      bic_gompertz = bic_gompertz,
      aic_gam = aic_gam,
      bic_gam = bic_gam,
      mantel_r = mantel_r,
      mantel_p = mantel_p
    )
  }

  list(
    results = bind_rows(model_results),
    predictions = bind_rows(predictions)
  )
}

build_predictions_long <- function(predictions) {
  predictions %>%
    tidyr::pivot_longer(
      cols = c("linear", "nls", "gam", "power", "gompertz"),
      names_to = "model",
      values_to = "predicted"
    ) %>%
    mutate(
      distance = case_when(
        metric == "diversity" ~ distance_km_diversity,
        metric == "log_abundance" ~ distance_km_abundance,
        metric == "log_connections" ~ distance_km_connections,
        TRUE ~ NA_real_
      )
    )
}

plot_model_comparison <- function(predictions_long, output_path) {
  plot <- ggplot(predictions_long, aes(x = distance, y = predicted, color = model)) +
    geom_line() +
    facet_grid(phylum ~ metric) +
    labs(
      title = "Model Comparison by Phylum and Metric",
      y = "Predicted Value (scaled)",
      x = "Distance (km)"
    ) +
    theme_minimal()

  ggsave(output_path, plot = plot, width = 8, height = 5)
  plot
}

compute_spearman_correlations <- function(data) {
  comparisons <- list(
    c("log_abundance", "diversity"),
    c("log_connections", "diversity"),
    c("log_abundance", "log_connections")
  )

  cor_global_df <- purrr::map_dfr(comparisons, function(vars) {
    var1 <- vars[1]
    var2 <- vars[2]

    tmp <- data %>%
      filter(!is.na(.data[[var1]]) & !is.na(.data[[var2]]))

    test <- cor.test(tmp[[var1]], tmp[[var2]], method = "spearman")

    tibble(
      phylum = "Global",
      variable_x = var1,
      variable_y = var2,
      rho = test$estimate,
      p_value = test$p.value,
      n = nrow(tmp)
    )
  })

  cor_phylum_df <- purrr::map_dfr(comparisons, function(vars) {
    var1 <- vars[1]
    var2 <- vars[2]

    data %>%
      filter(!is.na(.data[[var1]]) & !is.na(.data[[var2]])) %>%
      group_by(phylum) %>%
      summarise(
        variable_x = var1,
        variable_y = var2,
        rho = cor(.data[[var1]], .data[[var2]], method = "spearman"),
        p_value = cor.test(.data[[var1]], .data[[var2]], method = "spearman")$p.value,
        n = n(),
        .groups = "drop"
      )
  })

  cor_total <- bind_rows(cor_global_df, cor_phylum_df)

  list(
    comparisons = comparisons,
    cor_global_df = cor_global_df,
    cor_phylum_df = cor_phylum_df,
    cor_total = cor_total
  )
}

plot_spearman_bars <- function(cor_phylum_df, comparisons, output_dir) {
  for (pair in comparisons) {
    var1 <- pair[1]
    var2 <- pair[2]

    plot <- cor_phylum_df %>%
      filter(variable_x == var1, variable_y == var2) %>%
      ggplot(aes(x = reorder(phylum, -rho), y = rho, fill = p_value < 0.05)) +
      geom_col() +
      geom_text(aes(label = paste0("rho = ", round(rho, 2))), vjust = -0.5, size = 3) +
      scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray70"), name = "p < 0.05") +
      labs(
        title = paste("Spearman Correlation:", var1, "vs", var2, "by Phylum"),
        x = "Phylum",
        y = "rho (Spearman)"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    output_path <- file.path(
      output_dir,
      paste0("spearman_correlation_", var1, "_vs_", var2, ".pdf")
    )
    ggsave(output_path, plot = plot, device = "pdf", width = 8, height = 6)
  }
}

plot_spearman_heatmap <- function(cor_total, output_path) {
  plot <- ggplot(
    cor_total %>% filter(phylum != "Global"),
    aes(x = variable_x, y = variable_y, fill = rho)
  ) +
    geom_tile(color = "white") +
    geom_text(
      aes(label = paste0("rho=", round(rho, 2), ifelse(p_value < 0.05, "*", ""))),
      size = 3
    ) +
    scale_fill_gradient2(
      low = "red",
      mid = "white",
      high = "blue",
      midpoint = 0,
      name = "Spearman rho"
    ) +
    facet_wrap(~ phylum) +
    labs(
      title = "Heatmap of Correlations by Phylum",
      x = "Variable X",
      y = "Variable Y"
    ) +
    theme_minimal()

  ggsave(output_path, plot = plot, device = "pdf", width = 10, height = 8)
  plot
}

plot_spearman_network <- function(cor_global_df, output_path) {
  edges <- cor_global_df %>%
    select(from = variable_x, to = variable_y, rho)

  graph <- igraph::graph_from_data_frame(edges, directed = FALSE)

  plot <- ggraph(graph, layout = "circle") +
    geom_edge_link(aes(width = abs(rho), color = rho), alpha = 0.8) +
    geom_node_point(size = 5) +
    geom_node_text(aes(label = name), repel = TRUE, size = 4) +
    scale_edge_color_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0) +
    theme_void() +
    labs(title = "Global Correlation Network (Spearman rho)")

  ggsave(output_path, plot = plot, device = "pdf", width = 8, height = 8)
  plot
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

get_slope <- function(df, yvar, model_type = "lm", distance_col,
                      x_ref = NULL, gam_k = NULL,
                      eps_log = 1e-8, h_frac = 50) {
  y <- df[[yvar]]
  x <- df[[distance_col]]
  ok <- is.finite(y) & is.finite(x)
  y <- y[ok]
  x <- x[ok]

  n <- length(y)
  if (n < 5) {
    return(NA_real_)
  }

  if (is.null(x_ref)) {
    x_ref <- stats::median(x, na.rm = TRUE)
  }
  h <- stats::sd(x, na.rm = TRUE) / h_frac
  if (!is.finite(h) || h == 0) {
    h <- (max(x) - min(x)) / (2 * h_frac)
  }

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
      k_use <- gam_k %||% max(5, min(10, floor(n / 3)))
      fit <- mgcv::gam(y ~ s(x, k = k_use), method = "REML")
      pred <- function(xx) predict(fit, newdata = data.frame(x = xx), type = "response")
      (pred(x_ref + h) - pred(x_ref - h)) / (2 * h)
    } else if (model_type == "nls") {
      if (length(unique(x)) < 3) {
        return(NA_real_)
      }
      keep <- y > 0 & is.finite(y) & is.finite(x)
      y_raw <- y[keep]
      x_raw <- x[keep]
      if (length(y_raw) < 5) {
        return(NA_real_)
      }

      mu_x <- mean(x_raw)
      sd_x <- sd(x_raw)
      if (!is.finite(sd_x) || sd_x == 0) {
        return(NA_real_)
      }

      z <- (x_raw - mu_x) / sd_x
      z_ref <- (x_ref - mu_x) / sd_x
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
      slope <- (co[["a"]] * co[["b"]] * exp(co[["b"]] * z_ref)) * (1 / sd_x)
      unname(as.numeric(slope))
    } else if (model_type == "gompertz") {
      a_start <- max(y_pos, na.rm = TRUE)
      c_start <- 0.05
      b_start <- 1
      df_g <- data.frame(x = x, y = y_pos)

      fit <- nls(
        y ~ a * exp(-b * exp(-c * x)),
        data = df_g,
        start = list(a = a_start, b = b_start, c = c_start),
        control = nls.control(warnOnly = TRUE, maxiter = 500)
      )

      co <- coef(fit)
      yx <- co[["a"]] * exp(-co[["b"]] * exp(-co[["c"]] * x_ref))
      unname(yx * co[["b"]] * co[["c"]] * exp(-co[["c"]] * x_ref))
    } else {
      stop("Unknown model_type: ", model_type)
    }
  }, error = function(e) {
    warning(sprintf("get_slope(%s) failed: %s", model_type, e$message))
    NA_real_
  })
}

permutation_test_slope <- function(data, phylum, model_type,
                                   metric1 = "scaled_log_abundance",
                                   metric2 = "scaled_diversity",
                                   B = 999) {
  df <- data %>%
    filter(phylum == !!phylum) %>%
    filter(!is.na(.data[[metric1]]), !is.na(.data[[metric2]]))

  n <- nrow(df)
  if (n < 5) {
    return(tibble(phylum = phylum, model = model_type, observed_diff = NA, p_value = NA, n = n))
  }

  distance_for <- function(metric) case_when(
    metric == "scaled_log_abundance" ~ "distance_km_abundance",
    metric == "scaled_diversity" ~ "distance_km_diversity",
    metric == "scaled_log_connections" ~ "distance_km_connections",
    TRUE ~ NA_character_
  )

  dist1 <- distance_for(metric1)
  dist2 <- distance_for(metric2)

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
  p_val <- if (!is.finite(observed_diff) || length(valid_diffs) == 0) {
    NA_real_
  } else {
    mean(abs(valid_diffs) >= abs(observed_diff))
  }

  tibble(
    phylum = phylum,
    model = model_type,
    observed_diff = observed_diff,
    p_value = p_val,
    n = n
  )
}

run_metric_permutation_tests <- function(data, metric_pairs, models, B = 999) {
  purrr::pmap_dfr(
    tidyr::expand_grid(
      phylum = unique(data$phylum),
      model_type = models,
      pair = metric_pairs
    ),
    function(phylum, model_type, pair) {
      permutation_test_slope(
        data,
        phylum = phylum,
        model_type = model_type,
        metric1 = pair[[1]],
        metric2 = pair[[2]],
        B = B
      ) %>%
        mutate(metric1 = pair[[1]], metric2 = pair[[2]])
    }
  )
}

get_all_slopes <- function(data, model_type = "lm") {
  phyla <- unique(data$phylum)
  metrics <- c("scaled_log_abundance", "scaled_diversity", "scaled_log_connections")

  combos <- tidyr::expand_grid(phylum = phyla, metric = metrics)

  purrr::pmap_dfr(
    list(combos$phylum, combos$metric),
    function(phylum_name, metric_name) {
      df_sub <- data %>% filter(phylum == phylum_name)

      distance_col <- case_when(
        metric_name == "scaled_log_abundance" ~ "distance_km_abundance",
        metric_name == "scaled_diversity" ~ "distance_km_diversity",
        metric_name == "scaled_log_connections" ~ "distance_km_connections"
      )

      slope_val <- get_slope(df_sub, yvar = metric_name, model_type = model_type, distance_col = distance_col)
      slope_num <- suppressWarnings(as.numeric(slope_val[1]))

      tibble(
        phylum = phylum_name,
        metric = metric_name,
        model = model_type,
        slope = slope_num
      )
    }
  )
}

read_phylum_comparisons <- function(input_file) {
  required <- c("phylum1", "phylum2", "observed_diff", "p_value", "model", "metric")
  df <- readr::read_delim(input_file, delim = ",", show_col_types = FALSE)

  if (!all(required %in% names(df))) {
    df <- readr::read_delim(input_file, delim = ";", show_col_types = FALSE)
  }

  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Missing required columns in comparison file: ", paste(missing, collapse = ", "))
  }

  df$observed_diff <- as.numeric(df$observed_diff)
  df$p_value <- as.numeric(df$p_value)
  df
}

prepare_half_matrix <- function(data) {
  data %>%
    rowwise() %>%
    mutate(
      phylum_pair = list(sort(c(phylum1, phylum2))),
      phylum1_new = phylum_pair[[1]],
      phylum2_new = phylum_pair[[2]],
      observed_diff_new = ifelse(phylum1 == phylum1_new, observed_diff, -observed_diff),
      p_value_new = p_value
    ) %>%
    ungroup() %>%
    select(
      phylum1 = phylum1_new,
      phylum2 = phylum2_new,
      observed_diff = observed_diff_new,
      p_value = p_value_new
    )
}

plot_phylum_heatmap <- function(data, metric_name, model_name, output_path) {
  df_metric <- data %>%
    filter(metric == metric_name, model == model_name) %>%
    prepare_half_matrix() %>%
    mutate(
      signif_label = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ ""
      )
    )

  plot <- ggplot(df_metric, aes(x = phylum1, y = phylum2, fill = observed_diff)) +
    geom_tile(color = "white") +
    geom_text(aes(label = signif_label), color = "black", size = 3) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = paste("Heatmap -", metric_name, "(", model_name, ")"),
      x = "Phylum 1",
      y = "Phylum 2",
      fill = "Observed Diff",
      caption = "* p < 0.05   ** p < 0.01   *** p < 0.001"
    )

  ggsave(output_path, plot = plot, width = 8, height = 6)
  plot
}

write_comparison_matrices <- function(data, metric_name, model_name, output_dir) {
  df_metric <- data %>%
    filter(metric == metric_name, model == model_name)

  phyla <- sort(unique(c(df_metric$phylum1, df_metric$phylum2)))
  matrix_diff <- matrix(
    0,
    nrow = length(phyla),
    ncol = length(phyla),
    dimnames = list(phyla, phyla)
  )
  matrix_pval <- matrix(
    0,
    nrow = length(phyla),
    ncol = length(phyla),
    dimnames = list(phyla, phyla)
  )

  for (i in seq_len(nrow(df_metric))) {
    p1 <- df_metric$phylum1[i]
    p2 <- df_metric$phylum2[i]
    val_diff <- df_metric$observed_diff[i]
    val_pval <- df_metric$p_value[i]

    i1 <- match(p1, phyla)
    i2 <- match(p2, phyla)

    matrix_diff[i1, i2] <- val_diff
    matrix_diff[i2, i1] <- -val_diff
    matrix_pval[i1, i2] <- val_pval
    matrix_pval[i2, i1] <- val_pval
  }

  diff_path <- file.path(
    output_dir,
    paste0("phylum_comparisons_", model_name, "_", metric_name, "_observed_diff.csv")
  )
  pval_path <- file.path(
    output_dir,
    paste0("phylum_comparisons_", model_name, "_", metric_name, "_p_value.csv")
  )

  write.csv(matrix_diff, diff_path, row.names = TRUE)
  write.csv(matrix_pval, pval_path, row.names = TRUE)

  list(diff_path = diff_path, pval_path = pval_path)
}

plot_phylum_dendrogram <- function(diff_path, pval_path, metric_name, model_name, output_path) {
  observed_diff <- as.matrix(read.csv(diff_path, row.names = 1, check.names = FALSE))
  p_values <- as.matrix(read.csv(pval_path, row.names = 1, check.names = FALSE))
  mode(observed_diff) <- "numeric"

  signif_labels <- matrix("", nrow = nrow(p_values), ncol = ncol(p_values))
  signif_labels[p_values < 0.05] <- "*"
  signif_labels[p_values < 0.01] <- "**"
  signif_labels[p_values < 0.001] <- "***"

  rownames(signif_labels) <- rownames(observed_diff)
  colnames(signif_labels) <- colnames(observed_diff)

  lim <- max(abs(observed_diff), na.rm = TRUE)
  palette <- colorRampPalette(brewer.pal(11, "RdBu"))(100)
  breaks <- seq(-lim, lim, length.out = length(palette) + 1)

  pheatmap(
    observed_diff,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = signif_labels,
    number_color = "black",
    color = palette,
    breaks = breaks,
    fontsize_row = 8,
    fontsize_col = 8,
    main = paste("Heatmap -", metric_name, "(", model_name, ")"),
    filename = output_path,
    width = 9,
    height = 7
  )
}
