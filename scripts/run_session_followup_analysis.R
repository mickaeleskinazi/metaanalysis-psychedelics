#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(purrr)
})

source(here::here("R", "compat_map_groups.R"))
source(here::here("R", "data_ingest.R"))
source(here::here("R", "dose_response_models.R"))
source(here::here("R", "forest_plots.R"))
source(here::here("R", "dose_response_plots.R"))
source(here::here("R", "significance_tables.R"))
source(here::here("R", "session_followup_tables.R"))
source(here::here("R", "session_followup_dr_plots.R"))
source(here::here("R", "session_followup_forest_plots.R"))

default_ref_policies <- list(
  MDMA       = c("inactive_placebo", "active_non_psy_placebo", "active_placebo"),
  LSD        = c("inactive_placebo", "active_placebo", "active_non_psy_placebo"),
  PSILOCYBIN = c("inactive_placebo", "active_non_psy_placebo", "active_placebo"),
  AYAHUASCA  = c("inactive_placebo", "active_placebo", "active_non_psy_placebo"),
  .default   = c("inactive_placebo", "active_placebo", "active_non_psy_placebo")
)

run_session_followup_analysis <- function(
    data_xlsx = here::here("data", "Adverse-events-dose-v5.xlsx"),
    sheet = "Feuil1",
    out_dir_session = here::here("results", "session"),
    out_dir_followup = here::here("results", "follow_up"),
    out_dir_compare = here::here("results", "compare"),
    min_k = 2,
    fit_spline = TRUE
) {
  if (!file.exists(data_xlsx)) {
    stop("Data file not found: ", data_xlsx)
  }

  raw_all <- load_data(data_xlsx, sheet = sheet)

  build_and_save_window <- function(window_value, out_dir) {
    if (!nrow(raw_all)) return(NULL)
    es <- build_es_for_window(raw_all, window_value, ref_policies = default_ref_policies)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    dr_mol <- run_dr_by_molecule(es, min_k = min_k, fit_spline = fit_spline)
    dr_ae  <- run_dr_by_ae(es, min_k = min_k, fit_spline = fit_spline)

    tables_dir <- file.path(out_dir, "tables")
    dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
    save_significance_table(dr_mol$models, file.path(tables_dir, "significance_by_molecule_models.csv"))
    save_significance_table(dr_ae$models,  file.path(tables_dir, "significance_by_ae_models.csv"))
    save_agg_significance_by_molecule(dr_mol$models, file.path(tables_dir, "significance_agg_by_molecule.csv"))
    save_agg_significance_by_ae_molecule(dr_ae$models,  file.path(tables_dir, "significance_agg_by_ae_molecule.csv"))

    make_forest_summary_per_molecule(es, file.path(out_dir, "forest_plots_summary"))
    plot_dr_by_molecule_split(dr_mol$preds, dr_mol$models, file.path(out_dir, "dose_response/by_molecule_split"))

    list(es = es, dr_mol = dr_mol, dr_ae = dr_ae)
  }

  res_session  <- build_and_save_window("session", out_dir_session)
  res_followup <- build_and_save_window("follow_up", out_dir_followup)

  if (is.null(res_session) || is.null(res_followup)) {
    warning("Both session and follow-up windows are required for comparison outputs.")
    return(invisible(list(session = res_session, follow_up = res_followup)))
  }

  es_all <- bind_rows(
    res_session$es %>% mutate(time_window = "session"),
    res_followup$es %>% mutate(time_window = "follow_up")
  )
  dir.create(out_dir_compare, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir_compare, "dose_response"), recursive = TRUE, showWarnings = FALSE)

  message("Comparative forest plots …")
  forest_compare_all_molecules(
    es_all,
    outdir = file.path(out_dir_compare, "forest_by_ae_side_by_side")
  )
  forest_compare_all_molecules_combined(
    es_all,
    outfile = file.path(out_dir_compare, "forest_combined_all_molecules.pdf")
  )

  message("Comparative dose–response overlays …")
  dr_compare_all_molecules(
    preds_session  = res_session$dr_mol$preds,
    preds_followup = res_followup$dr_mol$preds,
    outfile        = file.path(out_dir_compare, "dose_response", "dr_session_vs_followup.pdf")
  )

  message("Tables comparing session vs follow-up …")
  make_dr_window_tables(
    es = es_all,
    out_dir = file.path(out_dir_compare, "tables"),
    min_k_per_window = min_k,
    min_k_total = max(4, 2 * min_k)
  )

  invisible(list(
    session = res_session,
    follow_up = res_followup,
    combined_es = es_all
  ))
}

if (identical(environment(), globalenv())) {
  run_session_followup_analysis()
}
