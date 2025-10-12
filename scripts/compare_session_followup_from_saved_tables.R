#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

fmt_p <- function(p){
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "<0.001",
                formatC(p, format = "f", digits = 3)))
}

sig_stars <- function(p){
  case_when(
    is.na(p)        ~ "",
    p < 0.001       ~ "***",
    p < 0.01        ~ "**",
    p < 0.05        ~ "*",
    TRUE            ~ ""
  )
}

read_csv_safe <- function(path){
  if (!file.exists(path)) {
    warning("Missing file: ", path)
    return(NULL)
  }
  tryCatch(readr::read_csv(path, show_col_types = FALSE),
           error = function(e) { warning("Read error on ", path, ": ", e$message); NULL })
}

keep_or_all <- function(df){
  if (is.null(df)) return(NULL)
  if ("model" %in% names(df) && any(df$model == "spline_df3")) {
    df %>% filter(model == "spline_df3")
  } else df
}

collapse_byMol <- function(df){
  if (is.null(df)) return(NULL)
  df %>%
    mutate(term = ifelse(is.na(term), "", term)) %>%
    filter(!grepl("intrcpt", term, ignore.case = TRUE)) %>%
    group_by(molecule) %>%
    summarise(
      k_total = suppressWarnings(max(k, na.rm = TRUE)),
      I2      = suppressWarnings(max(I2, na.rm = TRUE)),
      tau2    = suppressWarnings(max(tau2, na.rm = TRUE)),
      p_overall = suppressWarnings(
        if ("QMp" %in% names(.)) min(QMp, na.rm = TRUE) else min(pval, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    mutate(stars = sig_stars(p_overall))
}

collapse_byAE <- function(df){
  if (is.null(df)) return(NULL)
  df %>%
    mutate(term = ifelse(is.na(term), "", term)) %>%
    filter(!grepl("intrcpt", term, ignore.case = TRUE)) %>%
    group_by(ae_term, molecule) %>%
    summarise(
      k_total = suppressWarnings(max(k, na.rm = TRUE)),
      p_overall = suppressWarnings(
        if ("QMp" %in% names(.)) min(QMp, na.rm = TRUE) else min(pval, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    mutate(stars = sig_stars(p_overall))
}

count_sig_ae <- function(df_coll){
  if (is.null(df_coll)) return(tibble(molecule = character(), n_sig = integer()))
  df_coll %>%
    filter(!is.na(p_overall)) %>%
    mutate(sig = p_overall < 0.05) %>%
    group_by(molecule) %>%
    summarise(n_sig = sum(sig, na.rm = TRUE), .groups = "drop")
}

compare_session_followup_from_saved_tables <- function(
    dir_session = here::here("results", "session", "tables"),
    dir_followup = here::here("results", "follow_up", "tables"),
    out_dir = here::here("results", "compare", "tables")
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  f_s_byMol   <- file.path(dir_session,  "significance_by_molecule_models.csv")
  f_f_byMol   <- file.path(dir_followup, "significance_by_molecule_models.csv")
  f_s_byAE    <- file.path(dir_session,  "significance_by_ae_models.csv")
  f_f_byAE    <- file.path(dir_followup, "significance_by_ae_models.csv")
  f_s_aggMol  <- file.path(dir_session,  "significance_agg_by_molecule.csv")
  f_f_aggMol  <- file.path(dir_followup, "significance_agg_by_molecule.csv")
  f_s_aggAEM  <- file.path(dir_session,  "significance_agg_by_ae_molecule.csv")
  f_f_aggAEM  <- file.path(dir_followup, "significance_agg_by_ae_molecule.csv")

  byMol_s  <- keep_or_all(read_csv_safe(f_s_byMol))
  byMol_f  <- keep_or_all(read_csv_safe(f_f_byMol))
  byAE_s   <- keep_or_all(read_csv_safe(f_s_byAE))
  byAE_f   <- keep_or_all(read_csv_safe(f_f_byAE))
  aggMol_s <- keep_or_all(read_csv_safe(f_s_aggMol))
  aggMol_f <- keep_or_all(read_csv_safe(f_f_aggMol))
  aggAEM_s <- keep_or_all(read_csv_safe(f_s_aggAEM))
  aggAEM_f <- keep_or_all(read_csv_safe(f_f_aggAEM))

  byMol_s_coll <- collapse_byMol(byMol_s) %>% mutate(window = "session")
  byMol_f_coll <- collapse_byMol(byMol_f) %>% mutate(window = "follow_up")

  cmp_byMol <- full_join(
    byMol_s_coll %>% rename(k_session = k_total, I2_session = I2, tau2_session = tau2,
                            p_overall_session = p_overall, stars_session = stars),
    byMol_f_coll %>% rename(k_follow = k_total, I2_follow = I2, tau2_follow = tau2,
                            p_overall_follow = p_overall, stars_follow = stars),
    by = "molecule"
  )
  write_csv(cmp_byMol, file.path(out_dir, "compare_by_molecule_overall.csv"))

  byAE_s_coll <- collapse_byAE(byAE_s) %>% mutate(window = "session")
  byAE_f_coll <- collapse_byAE(byAE_f) %>% mutate(window = "follow_up")

  cmp_byAE <- full_join(
    byAE_s_coll %>% rename(k_session = k_total, p_session = p_overall, stars_session = stars),
    byAE_f_coll %>% rename(k_follow  = k_total, p_follow  = p_overall, stars_follow  = stars),
    by = c("ae_term","molecule")
  )
  write_csv(cmp_byAE, file.path(out_dir, "compare_by_ae_molecule.csv"))

  cmp_agg_byMol <- full_join(
    (aggMol_s %>% rename(k_session = k_total, QM_session = QM, p_session = QMp, stars_session = stars)),
    (aggMol_f %>% rename(k_follow  = k_total, QM_follow  = QM, p_follow  = QMp, stars_follow  = stars)),
    by = "molecule"
  )
  write_csv(cmp_agg_byMol, file.path(out_dir, "compare_agg_by_molecule.csv"))

  cmp_agg_byAEM <- full_join(
    (aggAEM_s %>% rename(k_session = k_total, QM_session = QM, p_session = QMp, stars_session = stars)),
    (aggAEM_f %>% rename(k_follow  = k_total, QM_follow  = QM, p_follow  = QMp, stars_follow  = stars)),
    by = c("ae_term","molecule")
  )
  write_csv(cmp_agg_byAEM, file.path(out_dir, "compare_agg_by_ae_molecule.csv"))

  sig_counts_session <- count_sig_ae(byAE_s_coll)
  sig_counts_follow  <- count_sig_ae(byAE_f_coll)

  topline <- cmp_agg_byMol %>%
    select(molecule, p_session, p_follow, stars_session, stars_follow) %>%
    left_join(sig_counts_session %>% rename(n_sig_session = n_sig), by = "molecule") %>%
    left_join(sig_counts_follow  %>% rename(n_sig_follow = n_sig),  by = "molecule") %>%
    mutate(
      p_session_fmt = fmt_p(p_session),
      p_follow_fmt  = fmt_p(p_follow)
    ) %>%
    select(
      molecule,
      p_session_fmt,
      stars_session,
      n_sig_session,
      p_follow_fmt,
      stars_follow,
      n_sig_follow
    )

  write_csv(topline, file.path(out_dir, "compare_topline_molecule.csv"))

  invisible(list(
    cmp_by_molecule = cmp_byMol,
    cmp_by_ae_molecule = cmp_byAE,
    cmp_agg_by_molecule = cmp_agg_byMol,
    cmp_agg_by_ae_molecule = cmp_agg_byAEM,
    topline = topline
  ))
}

if (!interactive() && sys.nframe() == 0) {
  compare_session_followup_from_saved_tables()
}
