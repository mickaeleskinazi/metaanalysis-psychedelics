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

# Default reference arm preferences -------------------------------------------------
default_ref_policies <- list(
  MDMA       = c("inactive_placebo", "active_non_psy_placebo", "active_placebo"),
  LSD        = c("inactive_placebo", "active_placebo", "active_non_psy_placebo"),
  PSILOCYBIN = c("inactive_placebo", "active_non_psy_placebo", "active_placebo"),
  AYAHUASCA  = c("inactive_placebo", "active_placebo", "active_non_psy_placebo"),
  .default   = c("inactive_placebo", "active_placebo", "active_non_psy_placebo")
)

run_main_analysis <- function(
    data_xlsx = here::here("data", "Adverse-events-dose-v5_follow_up.xlsx"),
    sheet = "Sheet1",
    out_dir = here::here("results", "main"),
    min_k = 2,
    fit_spline = TRUE,
    make_paper_tables = TRUE,
    paper_dir = here::here("results", "paper_tables")
) {
  if (!file.exists(data_xlsx)) {
    stop("Data file not found: ", data_xlsx)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  message("1) Load & harmonize …")
  raw <- load_data(data_xlsx, sheet = sheet)

  message("2) Build pairwise contrasts (reference-aware) …")
  contr <- build_pairwise_2x2(raw, ref_policies = default_ref_policies)

  message("3) Escalc …")
  es <- build_escalc(contr)

  message("4) Dose–response by molecule …")
  dr_mol <- run_dr_by_molecule(es, min_k = min_k, fit_spline = fit_spline)

  message("5) Dose–response by adverse event …")
  dr_ae  <- run_dr_by_ae(es, min_k = min_k, fit_spline = fit_spline)

  tables_dir <- file.path(out_dir, "tables")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

  message("6) Save significance tables …")
  save_significance_table(dr_mol$models, file.path(tables_dir, "significance_by_molecule_models.csv"))
  save_significance_table(dr_ae$models,  file.path(tables_dir, "significance_by_ae_models.csv"))
  save_agg_significance_by_molecule(dr_mol$models, file.path(tables_dir, "significance_agg_by_molecule.csv"))
  save_agg_significance_by_ae_molecule(dr_ae$models,  file.path(tables_dir, "significance_agg_by_ae_molecule.csv"))

  message("7) Forest plots …")
  make_forest_plots(es, file.path(out_dir, "forest_plots"))
  make_forest_plots_per_molecule_pdf(es, file.path(out_dir, "forest_plots_by_molecule"))
  make_forest_summary_per_molecule(es, file.path(out_dir, "forest_plots_summary"))

  message("8) Dose–response plots …")
  plot_dr_by_molecule_split(dr_mol$preds, dr_mol$models, file.path(out_dir, "dose_response/by_molecule_split"))
  plot_dr_per_molecule_across_ae_facets(dr_ae$preds, dr_ae$models, file.path(out_dir, "dose_response/by_molecule_facets"))
  plot_dr_per_ae_normalized_dose(dr_ae$preds, file.path(out_dir, "dose_response/by_ae_normalized"))

  message("9) Master plots …")
  master_dir <- file.path(out_dir, "master")
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
    message("10) Publication tables …")
    make_all_paper_tables(
      path_mol_agg    = file.path(tables_dir, "significance_agg_by_molecule.csv"),
      path_ae_mol_agg = file.path(tables_dir, "significance_agg_by_ae_molecule.csv"),
      path_models_ae  = file.path(tables_dir, "significance_by_ae_models.csv"),
      output_dir      = paper_dir
    )
  }

  invisible(list(
    raw = raw,
    contrasts = contr,
    escalc = es,
    dr_molecule = dr_mol,
    dr_ae = dr_ae
  ))
}

if (identical(environment(), globalenv())) {
  run_main_analysis()
}
