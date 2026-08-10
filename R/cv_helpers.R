library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(rsample)

make_site_folds <- function(dat, mode = c("group_vfold", "leave_one_site_out"),
                            v = 5, repeats = 10, seed = 42) {
  mode <- match.arg(mode)
  set.seed(seed)

  if (mode == "group_vfold") {
    return(group_vfold_cv(dat, group = site, v = v, repeats = repeats))
  }

  sites <- sort(unique(dat$site))
  splits <- lapply(sites, function(s) {
    assessment_idx <- which(dat$site == s)
    analysis_idx <- which(dat$site != s)
    make_splits(list(analysis = analysis_idx, assessment = assessment_idx), dat)
  })

  tibble(
    splits = splits,
    id = paste0("LOSO_", sites)
  )
}

make_block_folds <- function(dat, block_col) {
  blocks <- sort(unique(dat[[block_col]]))
  blocks <- blocks[!is.na(blocks)]

  splits <- lapply(blocks, function(b) {
    assessment_idx <- which(dat[[block_col]] == b)
    analysis_idx <- which(dat[[block_col]] != b)
    make_splits(list(analysis = analysis_idx, assessment = assessment_idx), dat)
  })

  tibble(
    splits = splits,
    id = paste0("Holdout_", block_col, "_", blocks)
  )
}

add_spatial_blocks <- function(dat, k = 5, seed = 42) {
  site_coords <- dat %>%
    distinct(site, Latitude, Longitude) %>%
    filter(is.finite(Latitude), is.finite(Longitude))

  set.seed(seed)
  km <- kmeans(scale(site_coords[, c("Latitude", "Longitude")]), centers = k, nstart = 50)

  site_coords <- site_coords %>%
    mutate(
      spatial_cluster = paste0("cluster_", km$cluster),
      lon_band = ntile(Longitude, 3),
      lon_band = recode(as.character(lon_band),
                        "1" = "west",
                        "2" = "central",
                        "3" = "east")
    )

  dat %>%
    left_join(site_coords %>% select(site, spatial_cluster, lon_band), by = "site")
}

summarise_cv_predictions <- function(predictions, obs_col = "logR_obs",
                                     pred_col = "logR_pred",
                                     model_col = "Model") {
  predictions %>%
    group_by(.data[[model_col]]) %>%
    group_modify(~ model_stats(.x[[obs_col]], .x[[pred_col]])) %>%
    ungroup() %>%
    rename(Model = all_of(model_col)) %>%
    select(Model, everything())
}

summarise_cv_by_fold <- function(predictions, obs_col = "logR_obs",
                                 pred_col = "logR_pred",
                                 model_col = "Model",
                                 fold_col = "fold_id") {
  predictions %>%
    group_by(.data[[model_col]], .data[[fold_col]]) %>%
    group_modify(~ model_stats(.x[[obs_col]], .x[[pred_col]])) %>%
    ungroup() %>%
    rename(Model = all_of(model_col), fold_id = all_of(fold_col)) %>%
    group_by(Model) %>%
    summarise(
      n_folds = n(),
      across(c(RMSE, Pearson_r, CCC, Bias, R2),
             list(mean = ~ mean(.x, na.rm = TRUE),
                  sd = ~ sd(.x, na.rm = TRUE)),
             .names = "{.col}_{.fn}"),
      .groups = "drop"
    )
}

mean_impute_from_train <- function(train_df, test_df, predictors) {
  train_out <- train_df
  test_out <- test_df

  for (nm in predictors) {
    fill <- mean(train_out[[nm]], na.rm = TRUE)
    if (!is.finite(fill)) fill <- 0
    train_out[[nm]][!is.finite(train_out[[nm]]) | is.na(train_out[[nm]])] <- fill
    test_out[[nm]][!is.finite(test_out[[nm]]) | is.na(test_out[[nm]])] <- fill
  }

  list(train = train_out, test = test_out)
}

scale_from_train <- function(train_df, test_df, predictors) {
  train_out <- train_df
  test_out <- test_df

  for (nm in predictors) {
    center <- mean(train_out[[nm]], na.rm = TRUE)
    spread <- sd(train_out[[nm]], na.rm = TRUE)
    if (!is.finite(center)) center <- 0
    if (!is.finite(spread) || spread == 0) spread <- 1
    train_out[[nm]] <- (train_out[[nm]] - center) / spread
    test_out[[nm]] <- (test_out[[nm]] - center) / spread
  }

  list(train = train_out, test = test_out)
}

prune_correlated_train <- function(train_df, candidates, cutoff = 0.80) {
  candidates <- candidates[candidates %in% names(train_df)]
  candidates <- candidates[vapply(train_df[candidates], is.numeric, logical(1))]

  if (length(candidates) < 2) return(candidates)

  complete_counts <- vapply(train_df[candidates], function(x) sum(is.finite(x)), numeric(1))
  candidates <- candidates[complete_counts >= 5]
  if (length(candidates) < 2) return(candidates)

  corr_mat <- suppressWarnings(cor(train_df[, candidates, drop = FALSE],
                                   use = "pairwise.complete.obs"))
  corr_mat[!is.finite(corr_mat)] <- 0
  high_corr <- caret::findCorrelation(corr_mat, cutoff = cutoff)
  if (length(high_corr) > 0) candidates[-high_corr] else candidates
}

add_fold_mechanistic_terms <- function(train_df, test_df, biomass_col = "biomass_log10",
                                      fixed_vh0 = NULL) {
  if (!is.null(fixed_vh0)) {
    Vh0 <- fixed_vh0
  } else {
    Vh0 <- median(train_df$Vh, na.rm = TRUE)
    if (!is.finite(Vh0) || Vh0 <= 0) {
      Vh0 <- median(train_df$Vh0_global, na.rm = TRUE)
    }
  }

  add_terms <- function(dat) {
    dat %>%
      mutate(
        Vh0_fold = Vh0,
        biomass_for_model = .data[[biomass_col]],
        log_npoc = log10(npoc),
        inv_npoc = 1 / npoc,
        # Intermediate mechanistic exponent, retained for R_mech only.
        thermo_accessibility = abs(stoichMet_donor_mean) / (Vh0_fold * npoc),
        R_mech_fold = exp(-thermo_accessibility) * biomass_for_model,
        logR_obs = safe_log10(respC),
        logR_mech = safe_log10(R_mech_fold),
        dR = logR_obs - logR_mech,
        mech_x_npoc = logR_mech * npoc
      )
  }

  list(train = add_terms(train_df), test = add_terms(test_df), Vh0 = Vh0)
}
