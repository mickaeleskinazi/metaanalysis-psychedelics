suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(metafor); library(tidyr); library(readr); library(stringr)
})

# stars helper
.sig_stars <- function(p){
  dplyr::case_when(
    is.na(p)            ~ "",
    p < 0.001           ~ "***",
    p < 0.01            ~ "**",
    p < 0.05            ~ "*",
    TRUE                ~ ""
  )
}

# 1) Per-window slope within each molecule: yi ~ dose_mg (REML)
dr_fit_per_window <- function(es, min_k_per_window = 2){
  stopifnot(all(c("molecule","dose_mg","time_window","yi","vi") %in% names(es)))
  es %>%
    filter(!is.na(dose_mg), is.finite(yi), is.finite(vi)) %>%
    group_by(molecule, time_window) %>%
    group_modify(~{
      dat <- .x
      if (nrow(dat) < min_k_per_window) return(NULL)
      m <- tryCatch(rma(yi ~ dose_mg, vi = vi, data = dat, method = "REML"), error = function(e) NULL)
      if (is.null(m)) return(NULL)
      tibble(
        molecule       = dat$molecule[[1]],
        time_window    = dat$time_window[[1]],
        k              = m$k,
        beta_dose      = as.numeric(coef(m)["dose_mg"]),
        se_dose        = as.numeric(m$se["dose_mg"]),
        z_dose         = as.numeric(m$zval["dose_mg"]),
        p_dose         = as.numeric(m$pval["dose_mg"]),
        ci_lb          = as.numeric(m$ci.lb[["dose_mg"]]),
        ci_ub          = as.numeric(m$ci.ub[["dose_mg"]]),
        I2             = tryCatch(m$I2, error=function(e) NA_real_),
        tau2           = tryCatch(m$tau2, error=function(e) NA_real_)
      )
    }) %>% ungroup() %>%
    mutate(stars = .sig_stars(p_dose))
}

# 2) Interaction test: yi ~ dose_mg * time_window within each molecule
dr_test_session_vs_followup <- function(es, min_k_total = 4){
  stopifnot(all(c("molecule","dose_mg","time_window","yi","vi") %in% names(es)))
  es %>%
    filter(!is.na(dose_mg), is.finite(yi), is.finite(vi)) %>%
    group_by(molecule) %>%
    group_modify(~{
      dat <- .x
      if (length(unique(dat$time_window)) < 2 || nrow(dat) < min_k_total) return(NULL)
      m <- tryCatch(rma(yi ~ dose_mg * time_window, vi = vi, data = dat, method = "REML"), error=function(e) NULL)
      if (is.null(m)) return(NULL)
      tibble(
        molecule          = dat$molecule[[1]],
        beta_dose_main    = as.numeric(coef(m)["dose_mg"]),
        beta_interaction  = as.numeric(coef(m)[grep("^dose_mg:time_window", names(coef(m)))]),
        p_interaction     = as.numeric(m$pval[grep("^dose_mg:time_window", names(m$pval))])
      )
    }) %>% ungroup() %>%
    mutate(stars_interaction = .sig_stars(p_interaction))
}

# 3) Convenience wrapper: compute & export all tables
make_dr_window_tables <- function(
    es,
    out_dir = "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/results_compare",
    min_k_per_window = 2,
    min_k_total      = 4
){
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  per_window <- dr_fit_per_window(es, min_k_per_window = min_k_per_window)
  write_csv(per_window, file.path(out_dir, "dr_per_window_by_molecule.csv"))
  
  # Wide table: one row per molecule, session vs follow_up side-by-side
  per_window_wide <- per_window %>%
    select(molecule, time_window, beta_dose, ci_lb, ci_ub, p_dose, stars, k, I2, tau2) %>%
    mutate(time_window = case_when(
      time_window %in% c("session","Session") ~ "session",
      time_window %in% c("follow_up","follow-up","Follow_up") ~ "follow_up",
      TRUE ~ time_window
    )) %>%
    pivot_wider(
      id_cols = molecule,
      names_from = time_window,
      values_from = c(beta_dose, ci_lb, ci_ub, p_dose, stars, k, I2, tau2),
      names_sep = "."
    )
  write_csv(per_window_wide, file.path(out_dir, "dr_per_window_by_molecule_wide.csv"))
  
  interact <- dr_test_session_vs_followup(es, min_k_total = min_k_total)
  write_csv(interact, file.path(out_dir, "dr_session_vs_followup_interaction.csv"))
  
  # Merge a compact publication table
  pub_tab <- per_window_wide %>%
    left_join(interact %>% select(molecule, beta_interaction, p_interaction, stars_interaction), by="molecule") %>%
    mutate(
      session_effect   = ifelse(!is.na(beta_dose.session),
                                sprintf("%.3f [%.3f, %.3f]%s", beta_dose.session, ci_lb.session, ci_ub.session,
                                        ifelse(is.na(stars.session), "", paste0(" ", stars.session))), NA),
      followup_effect  = ifelse(!is.na(beta_dose.follow_up),
                                sprintf("%.3f [%.3f, %.3f]%s", beta_dose.follow_up, ci_lb.follow_up, ci_ub.follow_up,
                                        ifelse(is.na(stars.follow_up), "", paste0(" ", stars.follow_up))), NA),
      interaction_str  = ifelse(!is.na(p_interaction),
                                sprintf("β×window=%.3f; p=%.3g%s",
                                        beta_interaction, p_interaction,
                                        ifelse(is.na(stars_interaction), "", paste0(" ", stars_interaction))), NA)
    ) %>%
    transmute(
      Molecule = molecule,
      `Session slope (β, 95% CI)`  = session_effect,
      `Follow-up slope (β, 95% CI)`= followup_effect,
      `Window interaction`         = interaction_str
    )
  
  write_csv(pub_tab, file.path(out_dir, "dr_session_followup_publication_table.csv"))
  
  invisible(list(per_window = per_window,
                 per_window_wide = per_window_wide,
                 interaction = interact,
                 publication_table = pub_tab))
}