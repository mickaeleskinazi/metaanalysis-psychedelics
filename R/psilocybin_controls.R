# =========================
# FILE: R/psilocybin_controls.R
# =========================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readxl)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(metafor)
  library(here)
})

# ---- small helpers -----------------------------------------------------------

.extract_year <- function(author_year){
  as.integer(str_extract(as.character(author_year), "(19|20)\\d{2}"))
}

.safe_dir <- function(path){
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

.p_to_stars <- function(p){
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

# ---- core control runner -----------------------------------------------------
# Requires your existing functions from the pipeline:
# - load_data()          [data_ingest.R]
# - build_pairwise_2x2() [data_ingest.R]
# - build_escalc()       [dose_response_models.R]
# - norm_window_label()  [run_main_analysis.R OR redefined here if needed]

run_psilocybin_controls <- function(
    data_xlsx,
    sheet = "Feuil1",
    out_dir,
    window = "session",
    era_cutoff = 2018,
    min_k = 2
){
  .safe_dir(out_dir)
  
  # --- load ------------------------------------------------------------------
  raw <- load_data(data_xlsx, sheet = sheet)
  
  if ("time_window" %in% names(raw)) {
    # if you already have norm_window_label() in global env use it
    if (exists("norm_window_label", mode = "function")) {
      raw <- raw %>% mutate(time_window = norm_window_label(time_window))
    } else {
      raw <- raw %>%
        mutate(
          time_window = tolower(trimws(as.character(time_window))),
          time_window = dplyr::case_when(
            time_window %in% c("session","acute","in-session","in session") ~ "session",
            time_window %in% c("followup","follow-up","follow up","follow_up") ~ "follow_up",
            TRUE ~ time_window
          )
        )
    }
  }
  
  if (!"author_year" %in% names(raw)) {
    stop("run_psilocybin_controls(): column 'author_year' not found in raw. Add it to the Excel or ingest.")
  }
  
  raw <- raw %>% mutate(study_year = .extract_year(author_year))
  
  # --- subset (psilocybin + window) -------------------------------------------
  raw_w <- raw %>%
    filter(time_window == window) %>%
    mutate(molecule = toupper(as.character(molecule))) %>%
    filter(molecule == "PSILOCYBIN")
  
  if (!nrow(raw_w)) stop("No PSILOCYBIN rows for window = ", window)
  
  # --- build ES (same logic as main) ------------------------------------------
  if (!exists("default_ref_policies", mode = "list") && !exists("default_ref_policies", mode = "environment")) {
    stop("run_psilocybin_controls(): 'default_ref_policies' not found (it is defined in run_main_analysis.R).")
  }
  
  contrasts <- build_pairwise_2x2(raw_w, ref_policies = default_ref_policies)
  if (is.null(contrasts) || !nrow(contrasts)) stop("No contrasts could be built for psilocybin.")
  
  es <- build_escalc(contrasts) %>%
    left_join(raw_w %>% distinct(study_id, author_year, study_year), by = "study_id") %>%
    filter(is.finite(yi), is.finite(vi), is.finite(dose_mg))
  
  # quick sanity
  readr::write_csv(es, file.path(out_dir, "psilocybin_session_es.csv"))
  
  # ----------------------------------------------------------------------------
  # 1) TOPLINE: global linear dose slope
  # ----------------------------------------------------------------------------
  fit_global <- tryCatch(
    metafor::rma(yi, vi, mods = ~ dose_mg, data = es, method = "REML"),
    error = function(e) NULL
  )
  
  topline <- if (is.null(fit_global)) {
    tibble(
      window = window, molecule = "PSILOCYBIN",
      beta = NA_real_, p_slope = NA_real_, stars = "",
      k = nrow(es),
      year_min = min(es$study_year, na.rm = TRUE),
      year_max = max(es$study_year, na.rm = TRUE),
      dose_min = min(es$dose_mg, na.rm = TRUE),
      dose_max = max(es$dose_mg, na.rm = TRUE)
    )
  } else {
    p <- as.numeric(fit_global$QMp)
    tibble(
      window = window, molecule = "PSILOCYBIN",
      beta = as.numeric(coef(fit_global)[2]),
      p_slope = p,
      stars = .p_to_stars(p),
      k = fit_global$k,
      I2 = fit_global$I2,
      tau2 = fit_global$tau2,
      year_min = min(es$study_year, na.rm = TRUE),
      year_max = max(es$study_year, na.rm = TRUE),
      dose_min = min(es$dose_mg, na.rm = TRUE),
      dose_max = max(es$dose_mg, na.rm = TRUE)
    )
  }
  readr::write_csv(topline, file.path(out_dir, "psilocybin_session_topline.csv"))
  
  # ----------------------------------------------------------------------------
  # 2) META-REG TERMS: add study_year + dose × era split
  # ----------------------------------------------------------------------------
  es <- es %>%
    mutate(
      era = dplyr::case_when(
        is.na(study_year) ~ NA_character_,
        study_year <= era_cutoff ~ paste0("≤", era_cutoff),
        TRUE ~ paste0(">", era_cutoff)
      )
    )
  
  # model A: dose only (same as topline but keeps term table format)
  mod_A <- tryCatch(metafor::rma(yi, vi, mods = ~ dose_mg, data = es, method = "REML"), error = function(e) NULL)
  
  # model B: year-adjusted (dose + study_year)
  mod_B <- tryCatch(
    metafor::rma(yi, vi, mods = ~ dose_mg + study_year, data = es %>% filter(!is.na(study_year)), method = "REML"),
    error = function(e) NULL
  )
  
  # model C: era interaction (dose * era), requires both eras present
  mod_C <- NULL
  if (length(na.omit(unique(es$era))) >= 2) {
    mod_C <- tryCatch(
      metafor::rma(yi, vi, mods = ~ dose_mg * era, data = es %>% filter(!is.na(era)), method = "REML"),
      error = function(e) NULL
    )
  }
  
  pack_terms <- function(m, model_name){
    if (is.null(m)) return(tibble())
    co <- as.data.frame(m$beta)
    se <- as.data.frame(m$se)
    z  <- as.data.frame(m$zval)
    p  <- as.data.frame(m$pval)
    tibble(
      model = model_name,
      term  = rownames(co),
      estimate = as.numeric(co[,1]),
      se = as.numeric(se[,1]),
      z = as.numeric(z[,1]),
      p = as.numeric(p[,1]),
      stars = .p_to_stars(as.numeric(p[,1])),
      k = m$k,
      QM = m$QM,
      QMp = m$QMp
    )
  }
  
  terms <- bind_rows(
    pack_terms(mod_A, "dose_only"),
    pack_terms(mod_B, "dose_plus_year"),
    pack_terms(mod_C, "dose_x_era")
  )
  
  readr::write_csv(terms, file.path(out_dir, "psilocybin_session_meta_reg_terms.csv"))
  
  # ----------------------------------------------------------------------------
  # 3) LEAVE-ONE-OUT by study_id (global slope robustness)
  # ----------------------------------------------------------------------------
  loo <- es %>%
    distinct(study_id) %>%
    pull(study_id) %>%
    purrr::map_dfr(function(sid){
      dat <- es %>% filter(study_id != sid)
      m <- tryCatch(metafor::rma(yi, vi, mods = ~ dose_mg, data = dat, method = "REML"), error = function(e) NULL)
      if (is.null(m)) return(tibble())
      p <- as.numeric(m$QMp)
      tibble(
        left_out_study = sid,
        beta = as.numeric(coef(m)[2]),
        p_slope = p,
        stars = .p_to_stars(p),
        k = m$k
      )
    }) %>%
    left_join(es %>% distinct(study_id, author_year, study_year), by = c("left_out_study" = "study_id")) %>%
    arrange(beta)
  
  readr::write_csv(loo, file.path(out_dir, "psilocybin_session_leave_one_out.csv"))
  
  # ----------------------------------------------------------------------------
  # FIGURES
  # ----------------------------------------------------------------------------
  
  # A) bubble scatter: yi vs dose, size ~ inverse SE, color ~ year
  figA <- ggplot(
    es %>% mutate(w = 1 / sqrt(vi)),
    aes(x = dose_mg, y = yi, size = w, color = study_year)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey50") +
    geom_point(alpha = 0.75) +
    scale_size_continuous(range = c(1.2, 6.0), breaks = pretty_breaks(4)) +
    scale_color_viridis_c(option = "C", end = 0.95, na.value = "grey70") +
    labs(
      title = "Psilocybin (session): effect sizes vs dose",
      subtitle = "Point size ~ precision (1/SE); color ~ study year",
      x = "Dose (mg)",
      y = "Log odds ratio (yi)",
      color = "Year",
      size = "1/SE"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
  
  ggsave(file.path(out_dir, "FigA_psilocybin_scatter_yi_vs_dose.pdf"), figA, width = 8.6, height = 5.2, units = "in", dpi = 300)
  
  # B) fitted lines by era (dose-only within era)
  pred_by_era <- es %>%
    filter(!is.na(era)) %>%
    group_by(era) %>%
    group_modify(function(dat, key){
      if (nrow(dat) < 4 || n_distinct(dat$dose_mg) < 2) return(tibble())
      m <- tryCatch(metafor::rma(yi, vi, mods = ~ dose_mg, data = dat, method = "REML"), error = function(e) NULL)
      if (is.null(m)) return(tibble())
      grid <- seq(min(dat$dose_mg), max(dat$dose_mg), length.out = 80)
      pr <- predict(m, newmods = grid)
      tibble(
        era = key$era,
        dose_mg = grid,
        fit = pr$pred,
        lwr = pr$ci.lb,
        upr = pr$ci.ub,
        beta = as.numeric(coef(m)[2]),
        p_slope = as.numeric(m$QMp)
      )
    }) %>%
    ungroup() %>%
    mutate(stars = .p_to_stars(p_slope))
  
  if (nrow(pred_by_era)) {
    figB <- ggplot(pred_by_era, aes(x = dose_mg, y = fit, color = era, fill = era)) +
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, linewidth = 0) +
      geom_line(linewidth = 1.1) +
      labs(
        title = "Psilocybin (session): dose–response split by era",
        subtitle = "Separate linear fits within era; ribbon = 95% CI",
        x = "Dose (mg)",
        y = "Predicted log OR (AEs)"
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom") +
      guides(fill = "none")
    
    # add stars to legend labels (keeps it readable)
    era_lab <- pred_by_era %>%
      group_by(era) %>%
      summarise(st = first(stars), .groups = "drop") %>%
      mutate(lab = paste0(era, " ", st))
    figB <- figB + scale_color_discrete(labels = setNames(era_lab$lab, era_lab$era))
    
    ggsave(file.path(out_dir, "FigB_psilocybin_dr_by_era.pdf"), figB, width = 8.6, height = 4.9, units = "in", dpi = 300)
  }
  
  # C) leave-one-out slope plot
  if (nrow(loo)) {
    figC <- ggplot(loo, aes(x = reorder(left_out_study, beta), y = beta)) +
      geom_hline(yintercept = if (!is.null(fit_global)) as.numeric(coef(fit_global)[2]) else 0,
                 linetype = "dashed", linewidth = 0.5, color = "grey40") +
      geom_point(aes(color = p_slope), size = 2.2) +
      coord_flip() +
      scale_color_viridis_c(option = "D", end = 0.95, direction = -1, na.value = "grey70") +
      labs(
        title = "Psilocybin (session): leave-one-out dose slope",
        subtitle = "Each point = slope after removing one study; dashed line = full-data slope",
        x = NULL,
        y = "Dose slope (beta)",
        color = "p-value"
      ) +
      theme_minimal(base_size = 11)
    
    ggsave(file.path(out_dir, "FigC_psilocybin_leave_one_out_beta.pdf"), figC, width = 8.6, height = 6.8, units = "in", dpi = 300)
  }
  
  invisible(list(
    es = es,
    topline = topline,
    terms = terms,
    leave_one_out = loo
  ))
}

# =============================================================================
# Compatibility wrapper for run_main_analysis.R
# Ensures make_paper_figures_for_window() exists and works with dr_ae objects
# =============================================================================

make_paper_figures_for_window <- function(
    window_value,
    es,
    dr_mol,
    dr_ae,
    out_dir,
    only_session = TRUE
){
  if (isTRUE(only_session) && window_value != "session") return(invisible(NULL))
  
  fig_dir <- file.path(out_dir, "paper_figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  
  # dr_mol may be list(preds, models) OR directly preds
  preds_mol <- if (is.list(dr_mol) && "preds" %in% names(dr_mol)) dr_mol$preds else dr_mol
  
  # dr_ae may be list(preds, models) OR directly preds
  preds_ae  <- if (is.list(dr_ae)  && "preds" %in% names(dr_ae))  dr_ae$preds  else dr_ae
  
  es <- es %>%
    mutate(year_bin = dplyr::case_when(
      is.na(study_year) ~ NA_character_,
      study_year <= 2019 ~ "≤2019",
      study_year <= 2022 ~ "2020–2022",
      TRUE ~ "≥2023"
    ))
  
  pred_by_bin <- es %>%
    filter(!is.na(year_bin)) %>%
    group_by(year_bin) %>%
    group_modify(function(dat, key){
      if (nrow(dat) < 4 || n_distinct(dat$dose_mg) < 2) return(tibble())
      m <- tryCatch(
        metafor::rma(yi, vi, mods = ~ dose_mg, data = dat, method = "REML"),
        error = function(e) NULL
      )
      if (is.null(m)) return(tibble())
      grid <- seq(min(dat$dose_mg), max(dat$dose_mg), length.out = 80)
      pr <- predict(m, newmods = grid)
      tibble(
        dose_mg = grid,
        fit = pr$pred,
        lwr = pr$ci.lb,
        upr = pr$ci.ub,
        p = as.numeric(m$QMp)
      )
    }) %>%
    ungroup()
  
  if (nrow(pred_by_bin)) {
    figY <- ggplot(pred_by_bin, aes(dose_mg, fit, color = year_bin, fill = year_bin)) +
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, linewidth = 0) +
      geom_line(linewidth = 1.1) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom") +
      guides(fill = "none") +
      labs(
        title = "Psilocybin (session): dose–response by study-year bin",
        subtitle = "Linear fits within bins; ribbon = 95% CI",
        x = "Dose (mg)", y = "Predicted log OR (AEs)"
      )
    ggsave(file.path(out_dir, "FigY_psilocybin_dr_by_year_bin.pdf"),
           figY, width = 8.6, height = 4.9, units = "in", dpi = 300)
  }
  
  # ---- Call your existing plotting functions (names must match your file) ----
  # Fig1 (session DR by molecule)
  if (exists("plot_fig1_dr_session_three_molecules", mode = "function")) {
    plot_fig1_dr_session_three_molecules(
      preds_mol  = preds_mol,
      es_session = es,
      outfile    = file.path(fig_dir, "Fig1_DR_session_3molecules.pdf")
    )
  }
  
  # Fig2 (forest main AEs)
  if (exists("plot_fig2_forest_main_aes_by_molecule", mode = "function")) {
    plot_fig2_forest_main_aes_by_molecule(
      es      = es,
      outfile = file.path(fig_dir, "Fig2_Forest_mainAEs_4molecules.pdf")
    )
  }
  
  # Fig3 (shared AEs DR normalized)
  if (exists("plot_fig3_dr_shared_aes_normalized", mode = "function")) {
    plot_fig3_dr_shared_aes_normalized(
      preds_ae = preds_ae,
      es       = es,
      outfile  = file.path(fig_dir, "Fig3_DR_sharedAEs_normalized.pdf")
    )
  }
  
  invisible(TRUE)
}

# --- compatibility alias (do not remove) --------------------------------------
if (!exists("make_paper_figures_for_window") && exists("make_paper_figures_for_session")) {
  make_paper_figures_for_window <- function(window_value, es, dr_mol, dr_ae, out_dir, only_session = TRUE){
    # window_value ignored: old wrapper expected "session" only anyway
    make_paper_figures_for_session(
      es_session            = es,
      dr_mol_session        = dr_mol,
      dr_ae_session_for_plot = dr_ae,
      out_dir_session       = out_dir
    )
  }
}