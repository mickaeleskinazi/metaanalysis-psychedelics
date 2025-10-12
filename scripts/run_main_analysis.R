#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(readxl)
  library(dplyr)
  library(purrr)
  library(stringr)
})

source(here::here("R", "data_ingest.R"))
source(here::here("R", "dose_response_models.R"))
source(here::here("R", "forest_plots.R"))
source(here::here("R", "dose_response_plots.R"))
source(here::here("R", "significance_tables.R"))
source(here::here("R", "master_plots.R"))
source(here::here("R", "publication_tables.R"))
source(here::here("R", "session_followup_tables.R"))
source(here::here("R", "session_followup_dr_plots.R"))
source(here::here("R", "session_followup_forest_plots.R"))

# Default reference arm preferences -------------------------------------------------
default_ref_policies <- list(
  MDMA       = c("inactive_placebo", "active_non_psy_placebo", "active_placebo"),
  LSD        = c("inactive_placebo", "active_placebo", "active_non_psy_placebo"),
  PSILOCYBIN = c("inactive_placebo", "active_non_psy_placebo", "active_placebo"),
  AYAHUASCA  = c("inactive_placebo", "active_placebo", "active_non_psy_placebo"),
  .default   = c("inactive_placebo", "active_placebo", "active_non_psy_placebo")
)

run_main_analysis <- function(
    data_xlsx = here::here("data", "Adverse-events-dose-v5.xlsx"),
    sheet = "Feuil1",
    out_dir = here::here("results", "main"),
    min_k = 2,
    fit_spline = TRUE,
    make_paper_tables = TRUE,
    paper_dir = here::here("results", "paper_tables"),
    compare_dir = file.path(out_dir, "compare")) {
  if (!file.exists(data_xlsx)) {
    stop("Data file not found: ", data_xlsx)
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  message("1) Load & harmonize …")
  raw <- load_data(data_xlsx, sheet = sheet)
  if (!"time_window" %in% names(raw)) {
    stop("Loaded data does not contain a 'time_window' column.")
  }
  raw <- raw %>% mutate(time_window = norm_window_label(time_window))

  available_windows <- unique(raw$time_window)
  if (!length(available_windows)) {
    stop("No time windows found in the dataset.")
  }

  window_order <- c("session", "follow_up")
  ordered_windows <- unique(c(intersect(window_order, available_windows),
                              setdiff(available_windows, window_order)))

  run_window <- function(window_value) {
    out_dir_window <- file.path(out_dir, window_value)
    paper_dir_window <- file.path(paper_dir, window_value)
    dir.create(out_dir_window, recursive = TRUE, showWarnings = FALSE)

    message(sprintf("→ Window '%s': build contrasts …", window_value))
    es <- build_es_for_window(raw, window_value, ref_policies = default_ref_policies)
    if (is.null(es) || !nrow(es)) {
      warning("No effect sizes constructed for window '", window_value, "'.")
      return(NULL)
    }

    message(sprintf("→ Window '%s': dose–response models …", window_value))
    dr_mol <- run_dr_by_molecule(es, min_k = min_k, fit_spline = fit_spline)
    dr_ae  <- run_dr_by_ae(es, min_k = min_k, fit_spline = fit_spline)

    tables_dir <- file.path(out_dir_window, "tables")
    dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
    sig_mol <- save_significance_table(dr_mol$models, file.path(tables_dir, "significance_by_molecule_models.csv"))
    sig_ae  <- save_significance_table(dr_ae$models,  file.path(tables_dir, "significance_by_ae_models.csv"))
    agg_mol <- save_agg_significance_by_molecule(dr_mol$models, file.path(tables_dir, "significance_agg_by_molecule.csv"))
    agg_ae  <- save_agg_significance_by_ae_molecule(dr_ae$models,  file.path(tables_dir, "significance_agg_by_ae_molecule.csv"))

    message(sprintf("→ Window '%s': forest plots …", window_value))
    make_forest_plots(es, file.path(out_dir_window, "forest_plots"))
    make_forest_plots_per_molecule_pdf(es, file.path(out_dir_window, "forest_plots_by_molecule"))
    make_forest_summary_per_molecule(es, file.path(out_dir_window, "forest_plots_summary"))

    message(sprintf("→ Window '%s': dose–response plots …", window_value))
    plot_dr_by_molecule_split(dr_mol$preds, dr_mol$models, file.path(out_dir_window, "dose_response/by_molecule_split"))
    plot_dr_per_molecule_across_ae_facets(dr_ae$preds, dr_ae$models, file.path(out_dir_window, "dose_response/by_molecule_facets"))
    plot_dr_per_ae_normalized_dose(dr_ae$preds, file.path(out_dir_window, "dose_response/by_ae_normalized"))

    message(sprintf("→ Window '%s': master plots …", window_value))
    master_dir <- file.path(out_dir_window, "master")
    dir.create(master_dir, recursive = TRUE, showWarnings = FALSE)
    plot_master_dr_by_molecule(dr_mol$preds, file.path(master_dir, "master_dr_by_molecule.pdf"))
    plot_master_forest_by_molecule(es, file.path(master_dir, "master_forest_by_molecule.pdf"), min_k = 3)
    plot_master_dr_by_ae(
      preds               = dr_ae$preds,
      models              = dr_ae$models,
      outfile             = file.path(master_dir, "master_dr_by_ae.pdf"),
      max_ae_per_molecule = 20,
      significant_only    = TRUE
    )

    if (isTRUE(make_paper_tables)) {
      message(sprintf("→ Window '%s': publication tables …", window_value))
      dir.create(paper_dir_window, recursive = TRUE, showWarnings = FALSE)
      make_all_paper_tables(
        path_mol_agg    = file.path(tables_dir, "significance_agg_by_molecule.csv"),
        path_ae_mol_agg = file.path(tables_dir, "significance_agg_by_ae_molecule.csv"),
        path_models_ae  = file.path(tables_dir, "significance_by_ae_models.csv"),
        output_dir      = paper_dir_window
      )
    }

    list(
      es = es,
      dr_molecule = dr_mol,
      dr_ae = dr_ae,
      tables = list(
        significance_by_molecule = sig_mol,
        significance_by_ae = sig_ae,
        agg_by_molecule = agg_mol,
        agg_by_ae_molecule = agg_ae
      ),
      out_dir = out_dir_window,
      paper_dir = paper_dir_window
    )
  }

  window_results <- purrr::map(ordered_windows, run_window)
  names(window_results) <- ordered_windows
  window_results <- purrr::compact(window_results)

  comparison_outputs <- NULL

  if (all(c("session", "follow_up") %in% names(window_results))) {
    message("→ Comparative outputs (session vs follow-up) …")
    dir.create(compare_dir, recursive = TRUE, showWarnings = FALSE)
    compare_tables_dir <- file.path(compare_dir, "tables")
    dir.create(compare_tables_dir, recursive = TRUE, showWarnings = FALSE)

    es_all <- dplyr::bind_rows(
      window_results$session$es %>% mutate(time_window = "session"),
      window_results$follow_up$es %>% mutate(time_window = "follow_up")
    )

    forest_compare_all_molecules(
      es_all,
      outdir = file.path(compare_dir, "forest_by_ae_side_by_side")
    )
    forest_compare_all_molecules_combined(
      es_all,
      outfile = file.path(compare_dir, "forest_combined_all_molecules.pdf")
    )

    dr_compare_all_molecules(
      preds_session  = window_results$session$dr_molecule$preds,
      preds_followup = window_results$follow_up$dr_molecule$preds,
      outfile        = file.path(compare_dir, "dose_response", "dr_session_vs_followup.pdf")
    )

    make_dr_window_tables(
      es = es_all,
      out_dir = compare_tables_dir,
      min_k_per_window = min_k,
      min_k_total = max(4, 2 * min_k)
    )

    session_agg_ae <- window_results$session$tables$agg_by_ae_molecule %>%
      mutate(
        p_session = p_overall,
        stars_session = stars,
        k_session = k_total,
        QM_session = QM,
        QMp_session = QMp
      ) %>%
      select(ae_term, molecule, p_session, stars_session, k_session, QM_session, QMp_session)

    follow_agg_ae <- window_results$follow_up$tables$agg_by_ae_molecule %>%
      mutate(
        p_follow = p_overall,
        stars_follow = stars,
        k_follow = k_total,
        QM_follow = QM,
        QMp_follow = QMp
      ) %>%
      select(ae_term, molecule, p_follow, stars_follow, k_follow, QM_follow, QMp_follow)

    ae_compare <- dplyr::full_join(session_agg_ae, follow_agg_ae, by = c("ae_term", "molecule")) %>%
      mutate(
        sig_session = !is.na(p_session) & p_session < 0.05,
        sig_follow  = !is.na(p_follow)  & p_follow  < 0.05,
        window_presence = dplyr::case_when(
          sig_session & sig_follow ~ "session_and_follow_up",
          sig_session & !sig_follow ~ "session_only",
          !sig_session & sig_follow ~ "follow_up_only",
          TRUE ~ "not_significant"
        )
      )
    readr::write_csv(ae_compare, file.path(compare_tables_dir, "ae_significance_by_window.csv"))

    ae_counts <- ae_compare %>%
      group_by(molecule) %>%
      summarise(
        n_session_only = sum(window_presence == "session_only", na.rm = TRUE),
        n_follow_only  = sum(window_presence == "follow_up_only", na.rm = TRUE),
        n_both         = sum(window_presence == "session_and_follow_up", na.rm = TRUE),
        n_not_sig      = sum(window_presence == "not_significant", na.rm = TRUE),
        .groups = "drop"
      )
    readr::write_csv(ae_counts, file.path(compare_tables_dir, "ae_significance_counts_by_molecule.csv"))

    session_agg_mol <- window_results$session$tables$agg_by_molecule %>%
      transmute(
        molecule,
        p_session = p_overall,
        stars_session = stars,
        k_session = k_total
      )
    follow_agg_mol <- window_results$follow_up$tables$agg_by_molecule %>%
      transmute(
        molecule,
        p_follow = p_overall,
        stars_follow = stars,
        k_follow = k_total
      )

    topline <- full_join(session_agg_mol, follow_agg_mol, by = "molecule") %>%
      left_join(ae_counts, by = "molecule") %>%
      mutate(
        p_session_fmt = ifelse(is.na(p_session), "", ifelse(p_session < 0.001, "<0.001", formatC(p_session, format = "f", digits = 3))),
        p_follow_fmt  = ifelse(is.na(p_follow), "", ifelse(p_follow < 0.001, "<0.001", formatC(p_follow, format = "f", digits = 3)))
      )
    readr::write_csv(topline, file.path(compare_tables_dir, "topline_session_followup_summary.csv"))

    comparison_outputs <- list(
      es_all = es_all,
      ae_compare = ae_compare,
      ae_counts = ae_counts,
      topline = topline
    )
  } else {
    message("⚠️ Skipping comparative outputs: both 'session' and 'follow_up' windows are required.")
  }

  invisible(list(
    raw = raw,
    windows = window_results,
    comparison = comparison_outputs
  ))
}

if (identical(environment(), globalenv())) {
  run_main_analysis()
}
