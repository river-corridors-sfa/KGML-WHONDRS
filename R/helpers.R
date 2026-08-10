library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(stringr)
library(ggplot2)

root_file <- function(...) {
  file.path("..", ...)
}

raw_file <- function(...) {
  file.path("raw_data", ...)
}

derived_file <- function(...) {
  file.path("derived", ...)
}

figure_file <- function(...) {
  file.path("figures", ...)
}

unify_site_id <- function(x) {
  x <- as.character(x)
  id <- stringr::str_extract(x, "[0-9]{4}")
  ifelse(is.na(id), x, id)
}

position_from_id <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    stringr::str_detect(x, "(_ICR\\.|-|_)U\\b|SED_INC-U|FCS-U|BULK-U") ~ "U",
    stringr::str_detect(x, "(_ICR\\.|-|_)M\\b|SED_INC-M|FCS-M|BULK-M") ~ "M",
    stringr::str_detect(x, "(_ICR\\.|-|_)D\\b|SED_INC-D|FCS-D|BULK-D") ~ "D",
    TRUE ~ NA_character_
  )
}

rescale01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r)) || diff(r) == 0) {
    return(ifelse(is.na(x), NA_real_, 0))
  }
  (x - r[1]) / diff(r)
}

model_stats <- function(obs, pred) {
  idx <- complete.cases(obs, pred)
  obs <- obs[idx]
  pred <- pred[idx]

  if (length(obs) < 3) {
    return(tibble(
      n = length(obs), RMSE = NA_real_, Pearson_r = NA_real_,
      CCC = NA_real_, Bias = NA_real_, R2 = NA_real_
    ))
  }

  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  r <- suppressWarnings(cor(obs, pred))
  vx <- var(obs)
  vy <- var(pred)
  ccc <- (2 * r * sqrt(vx) * sqrt(vy)) / (vx + vy + (mean(obs) - mean(pred))^2)

  tibble(
    n = length(obs),
    RMSE = sqrt(mean((obs - pred)^2)),
    Pearson_r = r,
    CCC = ccc,
    Bias = mean(pred - obs),
    R2 = ifelse(ss_tot > 0, 1 - ss_res / ss_tot, NA_real_)
  )
}

calc_d50 <- function(clay, silt, fine, med, coarse) {
  d_min <- c(clay = 0.0005, silt = 0.002, fine = 0.053, med = 0.25, coarse = 0.5)
  d_max <- c(clay = 0.002, silt = 0.053, fine = 0.25, med = 0.5, coarse = 2.0)
  pct_vec <- c(clay = clay, silt = silt, fine = fine, med = med, coarse = coarse)

  if (!all(is.finite(pct_vec)) || sum(pct_vec, na.rm = TRUE) <= 0) {
    return(NA_real_)
  }

  log_d_mid <- (log10(d_min) + log10(d_max)) / 2
  d_mid <- 10^log_d_mid
  frac <- pct_vec / sum(pct_vec, na.rm = TRUE)
  F <- cumsum(frac)
  i <- which(F >= 0.5)[1]

  if (length(i) == 0 || is.na(i)) return(NA_real_)

  if (i == 1) {
    F0 <- 0
    log_d0 <- log10(d_min[1])
  } else {
    F0 <- F[i - 1]
    log_d0 <- log10(d_mid[i - 1])
  }

  F1 <- F[i]
  log_d1 <- log10(d_mid[i])
  if (!is.finite(F1 - F0) || F1 == F0) return(NA_real_)

  10^(log_d0 + ((0.5 - F0) / (F1 - F0)) * (log_d1 - log_d0))
}

safe_log10 <- function(x, eps = 1e-8) {
  log10(pmax(x, eps))
}

puor <- list(
  orange = "#e08214",
  purple = "#542788",
  light_orange = "#fdb863",
  light_purple = "#8073ac",
  dark_orange = "#b35806",
  dark_purple = "#2d004b",
  pale_orange = "#fee0b6",
  pale_purple = "#b2abd2"
)

puor_seq <- function(n) {
  if (n <= 0) return(character())

  strong_puor <- c(
    puor$dark_purple,
    puor$purple,
    puor$light_purple,
    puor$light_orange,
    puor$orange,
    puor$dark_orange
  )

  if (n <= length(strong_puor)) {
    return(strong_puor[seq_len(n)])
  }

  grDevices::colorRampPalette(strong_puor)(n)
}

model_colors <- c(
  "Mechanistic" = puor$purple,
  "Mechanistic fixed Vh0" = puor$purple,
  "Multiple linear regression" = puor$light_purple,
  "Pure RF" = puor$orange,
  "Pure XGBoost" = puor$dark_orange,
  "KGML RF NPOC" = puor$pale_purple,
  "KGML RF expanded" = puor$light_purple,
  "KGML RF FTICR" = puor$light_orange,
  "KGML RF no accessibility" = puor$pale_purple,
  "KGML RF full" = puor$light_orange,
  "KGML RF interaction" = puor$dark_purple,
  "KGML XGBoost full" = puor$dark_orange,
  "KGML RF no accessibility fixed Vh0" = puor$pale_purple,
  "KGML RF full fixed Vh0" = puor$light_orange,
  "KGML RF interaction fixed Vh0" = puor$dark_purple,
  "KGML XGBoost full fixed Vh0" = puor$dark_orange,
  "Pure RF" = puor$orange,
  "KGML RF interaction fixed Vh0 = 0.10" = puor$dark_purple
)

puor_values <- function(x) {
  x <- as.character(x)
  vals <- model_colors[x]
  missing <- is.na(vals)
  if (any(missing)) {
    vals[missing] <- puor_seq(sum(missing))
  }
  unname(vals)
}

scale_color_puor <- function(values, ...) {
  ggplot2::scale_color_manual(values = setNames(puor_values(values), values), ...)
}

scale_fill_puor <- function(values, ...) {
  ggplot2::scale_fill_manual(values = setNames(puor_values(values), values), ...)
}

theme_review <- function() {
  base_theme <- if (requireNamespace("ggpubr", quietly = TRUE)) {
    ggpubr::theme_pubr(base_size = 11, border = TRUE)
  } else {
    ggplot2::theme_classic(base_size = 11)
  }

  base_theme +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.5),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.35),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.35),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.4),
      legend.key = ggplot2::element_rect(fill = "white", color = NA)
    )
}
