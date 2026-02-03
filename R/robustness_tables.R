suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

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
  
  # linear dose term
  lin <- models_df %>%
    filter(model == "linear", term == "dose_diff") %>%
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
    filter(str_detect(model, "^spline")) %>%
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
  
  lin <- models_df %>%
    filter(model == "linear", term == "dose_diff") %>%
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
    filter(str_detect(model, "^spline")) %>%
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
