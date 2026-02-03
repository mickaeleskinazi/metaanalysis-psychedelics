# =========================
# FILE: run_main_analysis.R
# =========================
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(readxl)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(readr)
})

# =============================================================================
# SOURCES
# =============================================================================
source(here::here("R", "compat_map_groups.R"))
source(here::here("R", "data_ingest.R"))            # load_data(), build_pairwise_2x2()
source(here::here("R", "dose_response_models.R"))   # build_escalc(), run_dr_*
source(here::here("R", "forest_plots.R"))
source(here::here("R", "dose_response_plots.R"))
source(here::here("R", "master_plots.R"))
source(here::here("R", "session_followup_dr_plots.R"))
source(here::here("R", "session_followup_forest_plots.R"))
source(here::here("R", "bubble_weight_plots.R"))
source(here::here("R", "psilocybin_controls.R"))
source(here::here("R", "publication_figures.R"))

# =============================================================================
# UTILITIES
# =============================================================================
safe_dir <- function(path){
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

norm_window_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    x %in% c("session","acute","in-session","in session") ~ "session",
    x %in% c("followup","follow-up","follow up","follow_up") ~ "follow_up",
    TRUE ~ x
  )
}

extract_study_year <- function(author_year){
  as.integer(stringr::str_extract(as.character(author_year), "(19|20)\\d{2}"))
}

build_es_for_window <- function(raw, window_value, ref_policies) {
  raw_w <- dplyr::filter(raw, time_window == window_value)
  if (!nrow(raw_w)) return(tibble())
  
  contrasts <- build_pairwise_2x2(raw_w, ref_policies = ref_policies)
  if (is.null(contrasts) || !nrow(contrasts)) return(tibble())
  
  es <- build_escalc(contrasts)
  
  if ("study_id" %in% names(es) && all(c("study_id","study_year") %in% names(raw_w))) {
    es <- es %>% left_join(raw_w %>% distinct(study_id, study_year), by = "study_id")
  }
  
  es
}

# =============================================================================
# DEFAULT REFERENCE POLICIES
# =============================================================================
default_ref_policies <- list(
  MDMA       = c("inactive_placebo","active_non_psy_placebo","active_placebo"),
  LSD        = c("inactive_placebo","active_non_psy_placebo","active_placebo"),
  PSILOCYBIN = c("inactive_placebo","active_non_psy_placebo","active_placebo"),
  AYAHUASCA  = c("inactive_placebo","active_non_psy_placebo","active_placebo"),
  .default   = c("inactive_placebo","active_non_psy_placebo","active_placebo")
)

# =============================================================================
# MAIN
# =============================================================================
run_main_analysis <- function(
    data_xlsx = here::here("data", "Adverse-events-dose-v5.xlsx"),
    sheet     = "Feuil1",
    out_dir   = here::here("results", "main"),
    min_k     = 2,
    df_spline = 3
) {
  if (!file.exists(data_xlsx)) stop("Data file not found: ", data_xlsx)
  safe_dir(out_dir)
  
  message("1) Load & harmonize data …")
  raw <- load_data(data_xlsx, sheet = sheet) %>%
    mutate(time_window = norm_window_label(time_window))
  
  if ("author_year" %in% names(raw)) {
    raw <- raw %>% mutate(study_year = extract_study_year(author_year))
  } else {
    raw <- raw %>% mutate(study_year = NA_integer_)
  }
  
  if (exists("run_psilocybin_controls", mode = "function")) {
    message("→ Controls: psilocybin (session) …")
    try(run_psilocybin_controls(
      data_xlsx   = data_xlsx,
      sheet       = sheet,
      out_dir     = here::here("results", "controls_psilocybin"),
      window      = "session",
      era_cutoff  = 2018
    ), silent = TRUE)
  }
  
  windows <- intersect(c("session","follow_up"), unique(raw$time_window))
  if (!length(windows)) stop("No valid time windows in data.")
  
  models_to_run <- list(
    linear = list(model = "linear", df_spline = NA_integer_),
    spline = list(model = "spline", df_spline = df_spline)
  )
  
  run_window <- function(w){
    
    message("→ Window '", w, "': build effect sizes")
    es <- build_es_for_window(raw, w, default_ref_policies)
    if (!nrow(es)) return(NULL)
    
    out_w <- file.path(out_dir, w)
    safe_dir(out_w)
    
    # Bubble weights (optional)
    try({
      bubble_dir <- file.path(out_w, "paper_figures")
      safe_dir(bubble_dir)
      plot_bubble_weights_by_molecule(
        es       = es,
        outfile  = file.path(bubble_dir, paste0("FigX_Bubble_weights_", w, ".pdf")),
        top_n_ae = 16,
        min_k_ae = 3
      )
    }, silent = TRUE)
    
    # Forest plots (classic) (optional)
    try({
      make_forest_plots(es, file.path(out_w, "forest_plots"))
      make_forest_plots_per_molecule_pdf(es, file.path(out_w, "forest_plots_by_molecule"))
      make_forest_summary_per_molecule(es, file.path(out_w, "forest_plots_summary"))
    }, silent = TRUE)
    
    dr <- purrr::imap(models_to_run, function(spec, tag){
      
      message("   → DR ", tag, " …")
      
      out_tag <- file.path(out_w, tag)
      safe_dir(out_tag)
      
     dr_mol_obs <- run_dr_by_molecule(es, min_k=min_k, grid="observed",
                                 model=spec$model, df_spline=spec$df_spline)

dr_mol_plot <- run_dr_by_molecule(es, min_k=min_k, grid="continuous", n_grid=200,
                                  model=spec$model, df_spline=spec$df_spline)

# AE always linear
dr_ae_obs  <- run_dr_by_ae(es, min_k=min_k, grid="observed",
                           model="linear", df_spline=NA_integer_)
dr_ae_plot <- run_dr_by_ae(es, min_k=min_k, grid="continuous", n_grid=120,
                           model="linear", df_spline=NA_integer_)
      
      # --- tables
      tab_dir <- file.path(out_tag, "tables")
      safe_dir(tab_dir)
      
      write_csv(dr_mol_obs$models, file.path(tab_dir, paste0("dr_models_by_molecule_", w, "_", tag, ".csv")))
      write_csv(dr_ae_obs$models,  file.path(tab_dir, paste0("dr_models_by_ae_", w, "_", tag, ".csv")))
      
      dr_ae_sig <- dr_ae_obs$models %>%
        mutate(
          p_adj_bh_global = p.adjust(pval, method = "BH")
        ) %>%
        group_by(molecule) %>%
        mutate(
          p_adj_bh_by_mol = p.adjust(pval, method = "BH")
        ) %>%
        ungroup() %>%
        mutate(
          sig_p05          = !is.na(pval) & pval < 0.05,
          sig_fdr_global05 = !is.na(p_adj_bh_global) & p_adj_bh_global < 0.05,
          sig_fdr_mol05    = !is.na(p_adj_bh_by_mol) & p_adj_bh_by_mol < 0.05
        ) %>%
        arrange(pval)
      
      write_csv(dr_ae_sig, file.path(tab_dir, paste0("dr_models_by_ae_", w, "_", tag, "_with_sig_flags.csv")))
      write_csv(filter(dr_ae_sig, sig_p05),
                file.path(tab_dir, paste0("dr_models_by_ae_", w, "_", tag, "_significant_p05.csv")))
      write_csv(filter(dr_ae_sig, sig_fdr_global05),
                file.path(tab_dir, paste0("dr_models_by_ae_", w, "_", tag, "_significant_fdr05_global.csv")))
      write_csv(filter(dr_ae_sig, sig_fdr_mol05),
                file.path(tab_dir, paste0("dr_models_by_ae_", w, "_", tag, "_significant_fdr05_by_molecule.csv")))
      
      # --- standard DR plots (optional)
      try({
        dr_dir <- file.path(out_tag, "dose_response")
        safe_dir(dr_dir)
        plot_dr_by_molecule_split(dr_mol_plot$preds, outdir = file.path(dr_dir, "by_molecule"))
        plot_dr_per_molecule_across_ae_facets(dr_ae_plot$preds, outdir = file.path(dr_dir, "by_ae_facets"))
        plot_dr_per_ae_normalized_dose(dr_ae_plot$preds, outdir = file.path(dr_dir, "by_ae_normalized"))
      }, silent = TRUE)
      
      # --- master plots (optional)
      try({
        master_dir <- file.path(out_tag, "master")
        safe_dir(master_dir)
        plot_master_dr_by_molecule(dr_mol_plot$preds, outfile = file.path(master_dir, paste0("master_dr_by_molecule_", tag, ".pdf")))
        plot_master_dr_by_ae(preds = dr_ae_plot$preds, outfile = file.path(master_dir, paste0("master_dr_by_ae_", tag, ".pdf")), significant_only = FALSE)
      }, silent = TRUE)
      
      # --- paper figs (session only, per model)
      if (w == "session") {
        make_paper_figures_for_window(
          window_value = "session",
          es = es,
          dr_mol_plot = dr_mol_plot,
          dr_ae_plot  = dr_ae_plot,
          out_dir = out_tag,
          model = spec$model,
          df_spline = df_spline,
          only_session = TRUE
        )
      }
      
      list(mol_obs = dr_mol_obs, ae_obs = dr_ae_obs, mol_plot = dr_mol_plot, ae_plot = dr_ae_plot)
    })
    
    list(es = es, dr = dr)
  }
  
  results <- purrr::map(windows, run_window)
  names(results) <- windows
  results <- purrr::compact(results)
  
  if (all(c("session","follow_up") %in% names(results))) {
    message("→ Paper comparison figures (session vs follow-up): linear + spline")
    make_paper_figures_comparison(results, out_dir = out_dir, df_spline = df_spline)
  } else {
    message("⚠️ Follow-up window missing: comparison not produced.")
  }
  
  invisible(results)
}

if (identical(environment(), globalenv())) {
  run_main_analysis()
}