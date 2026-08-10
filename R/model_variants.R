library(dplyr)
library(tidyr)
library(ranger)
library(xgboost)

core_predictors <- c(
  "npoc", "inv_npoc", "log_npoc",
  "biomass_log10", "delGcat_mean", "lambda_mean", "ne_mean",
  "stoichMet_donor_mean", "clay", "silt", "fine", "med", "coarse",
  "ssa_tot", "d50_mm"
)

no_access_predictors <- c(
  "biomass_log10", "delGcat_mean", "lambda_mean", "ne_mean",
  "stoichMet_donor_mean", "clay", "silt", "fine", "med", "coarse",
  "ssa_tot", "d50_mm"
)

kgml_access_predictors <- c(
  "logR_mech", "npoc", "inv_npoc", "log_npoc",
  "biomass_log10", "delGcat_mean", "lambda_mean", "ne_mean",
  "stoichMet_donor_mean", "clay", "silt", "fine", "med", "coarse",
  "ssa_tot", "d50_mm"
)

kgml_interaction_predictors <- c(
  kgml_access_predictors,
  "mech_x_npoc"
)

fticr_predictors <- c(
  "n_peaks", "shannon",
  "NOSC_med", "NOSC_iqr", "NOSC_mean", "NOSC_p10", "NOSC_p90", "NOSC_range90",
  "DBE_med", "DBE_mean", "DBE_iqr",
  "AI_med", "AI_iqr",
  "delG_med", "delG_iqr",
  "lambda_med", "lambda_iqr",
  "Cnum_med", "Cnum_iqr",
  "frac_NOSC_pos", "frac_NOSC_low", "frac_AI_high", "DBE_C_ratio_mean",
  "frac_ConHC", "frac_Lipid", "frac_AminoSugar", "frac_Protein", "frac_Tannin",
  "frac_UnsatHC", "frac_Lignin", "frac_Carb", "frac_Other",
  "ratio_Protein_Lignin", "ratio_Lipid_Carb", "ratio_UnsatHC_ConHC"
)

fit_predict_lm <- function(train_df, test_df, response, predictors) {
  predictors <- predictors[predictors %in% names(train_df)]
  fixed <- mean_impute_from_train(train_df, test_df, predictors)
  formula <- reformulate(predictors, response = response)
  fit <- lm(formula, data = fixed$train)
  as.numeric(predict(fit, newdata = fixed$test))
}

fit_predict_rf <- function(train_df, test_df, response, predictors,
                           num.trees = 800, seed = 7) {
  predictors <- predictors[predictors %in% names(train_df)]
  fixed <- mean_impute_from_train(train_df, test_df, predictors)
  formula <- reformulate(predictors, response = response)
  fit <- ranger(
    formula = formula,
    data = fixed$train,
    num.trees = num.trees,
    min.node.size = 5,
    importance = "permutation",
    seed = seed
  )
  predict(fit, data = fixed$test)$predictions
}

fit_predict_xgb <- function(train_df, test_df, response, predictors,
                            nrounds = 120, seed = 7) {
  predictors <- predictors[predictors %in% names(train_df)]
  fixed <- mean_impute_from_train(train_df, test_df, predictors)
  scaled <- scale_from_train(fixed$train, fixed$test, predictors)

  dtrain <- xgb.DMatrix(
    data = as.matrix(scaled$train[, predictors, drop = FALSE]),
    label = scaled$train[[response]]
  )
  dtest <- xgb.DMatrix(data = as.matrix(scaled$test[, predictors, drop = FALSE]))

  set.seed(seed)
  fit <- xgb.train(
    params = list(
      objective = "reg:squarederror",
      eta = 0.04,
      max_depth = 2,
      subsample = 0.8,
      colsample_bytree = 0.8,
      min_child_weight = 3
    ),
    data = dtrain,
    nrounds = nrounds,
    verbose = 0
  )

  as.numeric(predict(fit, dtest))
}

run_revision_cv <- function(dat, folds, include_fticr = TRUE,
                            cor_cutoff = 0.80, fixed_vh0 = NULL) {
  results <- vector("list", nrow(folds))

  for (i in seq_len(nrow(folds))) {
    train_raw <- analysis(folds$splits[[i]])
    test_raw <- assessment(folds$splits[[i]])
    fold_id <- if ("id" %in% names(folds)) as.character(folds$id[[i]]) else paste0("Fold", i)

    terms <- add_fold_mechanistic_terms(train_raw, test_raw, fixed_vh0 = fixed_vh0)
    train_df <- terms$train
    test_df <- terms$test

    fticr_keep <- character(0)
    if (include_fticr) {
      fticr_keep <- prune_correlated_train(train_df, fticr_predictors, cutoff = cor_cutoff)
    }

    pure_predictors <- unique(c(core_predictors, fticr_keep))
    no_access <- unique(c(no_access_predictors, fticr_keep))
    full_kgml <- unique(c(kgml_access_predictors, fticr_keep))
    interaction_kgml <- unique(c(kgml_interaction_predictors, fticr_keep))

    base_cols <- test_df %>%
      transmute(
        site, Position, fold_id = fold_id, Vh0_fold = terms$Vh0,
        respC, npoc, logR_obs, logR_mech
      )

    fold_predictions <- bind_rows(
      base_cols %>%
        mutate(Model = "Mechanistic", logR_pred = logR_mech),

      base_cols %>%
        mutate(Model = "Multiple linear regression",
               logR_pred = fit_predict_lm(train_df, test_df, "logR_obs", pure_predictors)),

      base_cols %>%
        mutate(Model = "Pure RF",
               logR_pred = fit_predict_rf(train_df, test_df, "logR_obs", pure_predictors)),

      base_cols %>%
        mutate(Model = "Pure XGBoost",
               logR_pred = fit_predict_xgb(train_df, test_df, "logR_obs", pure_predictors)),

      base_cols %>%
        mutate(Model = "KGML RF no accessibility",
               logR_pred = logR_mech + fit_predict_rf(train_df, test_df, "dR", no_access)),

      base_cols %>%
        mutate(Model = "KGML RF full",
               logR_pred = logR_mech + fit_predict_rf(train_df, test_df, "dR", full_kgml)),

      base_cols %>%
        mutate(Model = "KGML RF interaction",
               logR_pred = logR_mech + fit_predict_rf(train_df, test_df, "dR", interaction_kgml)),

      base_cols %>%
        mutate(Model = "KGML XGBoost full",
               logR_pred = logR_mech + fit_predict_xgb(train_df, test_df, "dR", full_kgml))
    ) %>%
      mutate(
        R_obs = 10^logR_obs,
        R_pred = 10^logR_pred,
        residual_log = logR_obs - logR_pred
      )

    results[[i]] <- fold_predictions
  }

  bind_rows(results)
}

fit_final_rf <- function(dat, response, predictors, cor_cutoff = 0.80) {
  predictors <- predictors[predictors %in% names(dat)]
  fticr_in_model <- intersect(predictors, fticr_predictors)
  non_fticr <- setdiff(predictors, fticr_predictors)
  fticr_keep <- prune_correlated_train(dat, fticr_in_model, cutoff = cor_cutoff)
  predictors <- unique(c(non_fticr, fticr_keep))
  fixed <- mean_impute_from_train(dat, dat, predictors)$train

  fit <- ranger(
    formula = reformulate(predictors, response = response),
    data = fixed,
    num.trees = 1000,
    min.node.size = 5,
    importance = "permutation",
    seed = 7
  )

  list(model = fit, data = fixed, predictors = predictors)
}
