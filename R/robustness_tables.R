suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

normalize_robustness_models <- function(models_df){
  if (is.null(models_df)) return(tibble::tibble())
  models_df <- tibble::as_tibble(models_df, .name_repair = "minimal")
  if (!nrow(models_df)) return(models_df)

  needed_num <- c("k","I2","tau2","QM","QMp","pval","beta")
  for (nm in needed_num) {
    if (!nm %in% names(models_df)) models_df[[nm]] <- NA_real_
  }
  if (!"model" %in% names(models_df)) models_df$model <- NA_character_
  if (!"df_spline" %in% names(models_df)) models_df$df_spline <- NA_real_
  if (!"term" %in% names(models_df)) models_df$term <- "dose_diff"
  if (!"estimate" %in% names(models_df)) {
    models_df$estimate <- if ("beta" %in% names(models_df)) models_df$beta else NA_real_
  }

  models_df %>%
    mutate(
      model = as.character(.data$model),
      term = as.character(.data$term),
      model = case_when(
        !is.na(.data$df_spline) & str_detect(.data$model, "^spline$") ~ paste0("spline_df", as.integer(.data$df_spline)),
        TRUE ~ .data$model
      ),
      term = coalesce(.data$term, "dose_diff")
    )
}

sig_stars <- function(p){
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

# ---- molecule-level --------------------------------------------------------

robustness_molecule_linear_vs_spline <- function(models_df, outfile_csv, alpha = 0.05){
  stopifnot(is.data.frame(models_df))
  dir.create(dirname(outfile_csv), recursive = TRUE, showWarnings = FALSE)
  if (is.null(models_df) || !nrow(models_df)) {
    out <- tibble::tibble()
    write_csv(out, outfile_csv)
    return(invisible(out))
  }
  models_df <- normalize_robustness_models(models_df)
  if (!"QM" %in% names(models_df)) models_df$QM <- NA_real_
  if (!"QMp" %in% names(models_df)) models_df$QMp <- NA_real_
  if (!"I2" %in% names(models_df)) models_df$I2 <- NA_real_
  if (!"tau2" %in% names(models_df)) models_df$tau2 <- NA_real_
  
  # linear dose term
  lin <- models_df %>%
    filter(.data$model == "linear", .data$term == "dose_diff") %>%
    mutate(
      I2 = if ("I2" %in% names(.)) suppressWarnings(as.numeric(I2)) else NA_real_,
      tau2 = if ("tau2" %in% names(.)) suppressWarnings(as.numeric(tau2)) else NA_real_
    ) %>%
    transmute(
      molecule,
      k_linear = k,
      I2_linear = I2,
      tau2_linear = tau2,
      beta_linear = estimate,
      p_linear = pval,
      stars_linear = sig_stars(p_linear),
      dir_linear = case_when(
        is.na(beta_linear) ~ NA_character_,
        beta_linear > 0 ~ "positive",
        beta_linear < 0 ~ "negative",
        TRUE ~ "zero"
      )
    )
  
  # spline omnibus (take one row per molecule per spline model)
  spl <- models_df %>%
    filter(str_detect(.data$model, "^spline")) %>%
    mutate(
      I2 = if ("I2" %in% names(.)) suppressWarnings(as.numeric(I2)) else NA_real_,
      tau2 = if ("tau2" %in% names(.)) suppressWarnings(as.numeric(tau2)) else NA_real_,
      QM = if ("QM" %in% names(.)) suppressWarnings(as.numeric(QM)) else NA_real_,
      QMp = if ("QMp" %in% names(.)) suppressWarnings(as.numeric(QMp)) else NA_real_
    ) %>%
    group_by(molecule, model) %>%
    summarise(
      k_spline = first(k),
      I2_spline = first(I2),
      tau2_spline = first(tau2),
      QM_spline = first(QM),
      p_spline = first(QMp),
      stars_spline = sig_stars(p_spline),
      .groups = "drop"
    ) %>%
    mutate(df_spline = str_extract(model, "df\\d+") %>% str_remove("df") %>% suppressWarnings(as.integer(.))) %>%
    arrange(molecule, df_spline) %>%
    group_by(molecule) %>%
    slice_head(n = 1) %>%   # keep the first spline spec (usually df3). change if you want "best"
    ungroup() %>%
    select(molecule, spline_model = model, df_spline, k_spline, I2_spline, tau2_spline, QM_spline, p_spline, stars_spline)
  
  out <- full_join(lin, spl, by = "molecule") %>%
    mutate(
      sig_linear = !is.na(p_linear) & p_linear < alpha,
      sig_spline = !is.na(p_spline) & p_spline < alpha,
      robustness = case_when(
        sig_linear & sig_spline ~ "both_significant",
        sig_linear & !sig_spline ~ "linear_only",
        !sig_linear & sig_spline ~ "spline_only",
        TRUE ~ "neither_significant"
      )
    ) %>%
    arrange(molecule)
  
  write_csv(out, outfile_csv)
  invisible(out)
}

# ---- AE × molecule ---------------------------------------------------------

robustness_ae_molecule_linear_vs_spline <- function(models_df, outfile_csv, alpha = 0.05){
  stopifnot(is.data.frame(models_df))
  dir.create(dirname(outfile_csv), recursive = TRUE, showWarnings = FALSE)
  if (is.null(models_df) || !nrow(models_df)) {
    out <- tibble::tibble()
    write_csv(out, outfile_csv)
    return(invisible(out))
  }
  models_df <- normalize_robustness_models(models_df)
  if (!"QM" %in% names(models_df)) models_df$QM <- NA_real_
  if (!"QMp" %in% names(models_df)) models_df$QMp <- NA_real_
  
  lin <- models_df %>%
    filter(.data$model == "linear", .data$term == "dose_diff") %>%
    transmute(
      ae_term, molecule,
      k_linear = k,
      beta_linear = estimate,
      p_linear = pval,
      stars_linear = sig_stars(p_linear),
      dir_linear = case_when(
        is.na(beta_linear) ~ NA_character_,
        beta_linear > 0 ~ "positive",
        beta_linear < 0 ~ "negative",
        TRUE ~ "zero"
      )
    )
  
  spl <- models_df %>%
    filter(str_detect(.data$model, "^spline")) %>%
    mutate(
      QM = if ("QM" %in% names(.)) suppressWarnings(as.numeric(QM)) else NA_real_,
      QMp = if ("QMp" %in% names(.)) suppressWarnings(as.numeric(QMp)) else NA_real_
    ) %>%
    group_by(ae_term, molecule, model) %>%
    summarise(
      k_spline = first(k),
      QM_spline = first(QM),
      p_spline = first(QMp),
      stars_spline = sig_stars(p_spline),
      .groups = "drop"
    ) %>%
    mutate(df_spline = str_extract(model, "df\\d+") %>% str_remove("df") %>% suppressWarnings(as.integer(.))) %>%
    arrange(ae_term, molecule, df_spline) %>%
    group_by(ae_term, molecule) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(ae_term, molecule, spline_model = model, df_spline, k_spline, QM_spline, p_spline, stars_spline)
  
  out <- full_join(lin, spl, by = c("ae_term", "molecule")) %>%
    mutate(
      sig_linear = !is.na(p_linear) & p_linear < alpha,
      sig_spline = !is.na(p_spline) & p_spline < alpha,
      robustness = case_when(
        sig_linear & sig_spline ~ "both_significant",
        sig_linear & !sig_spline ~ "linear_only",
        !sig_linear & sig_spline ~ "spline_only",
        TRUE ~ "neither_significant"
      )
    ) %>%
    arrange(molecule, ae_term)
  
  write_csv(out, outfile_csv)
  invisible(out)
}
