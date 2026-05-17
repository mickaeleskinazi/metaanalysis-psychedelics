#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(purrr)
  library(stringr)
})

source(here::here("R", "compat_map_groups.R"))
source(here::here("R", "data_ingest.R"))
source(here::here("R", "dose_response_models.R"))
source(here::here("R", "significance_tables.R"))
source(here::here("R", "publication_tables.R"))
source(here::here("R", "session_followup_tables.R"))
source(here::here("R", "publication_figures.R"))
source(here::here("R", "absolute_rate_tables.R"))


# Default reference arm preferences -------------------------------------------------
default_ref_policies <- list(
  MDMA       = c("inactive_placebo", "active_non_psy_placebo", "active_placebo", "placebo_actif", "placebo_inactif", "placebo_actif_non_psychedelic"),
  LSD        = c("inactive_placebo", "active_non_psy_placebo", "active_placebo", "placebo_actif", "placebo_inactif", "placebo_actif_non_psychedelic"),
  PSILOCYBIN = c("inactive_placebo", "active_non_psy_placebo", "active_placebo", "placebo_actif", "placebo_inactif", "placebo_actif_non_psychedelic"),
  AYAHUASCA  = c("inactive_placebo", "active_non_psy_placebo", "active_placebo", "placebo_actif", "placebo_inactif", "placebo_actif_non_psychedelic"),
  .default   = c("inactive_placebo", "active_non_psy_placebo", "active_placebo", "placebo_actif", "placebo_inactif", "placebo_actif_non_psychedelic")
)

select_model_preds <- function(dr_result, preferred = c("spline", "linear")) {
  preferred <- match.arg(preferred)
  if (is.null(dr_result) || is.null(dr_result$preds) || !nrow(dr_result$preds)) {
    return(tibble::tibble())
  }

  models_available <- unique(as.character(dr_result$preds$model))
  model_to_use <- if (preferred %in% models_available) preferred else models_available[[1]]
  dr_result$preds %>% filter(.data$model == model_to_use)
}

select_model_preds_by_molecule <- function(preds,
                                           molecules,
                                           preferred = c("spline", "linear")) {
  if (is.null(preds) || !nrow(preds)) return(tibble::tibble())

  preds <- preds %>%
    mutate(molecule = toupper(as.character(.data$molecule)))

  purrr::map_dfr(molecules, function(molecule_value) {
    dat <- preds %>% filter(.data$molecule == molecule_value)
    if (!nrow(dat)) return(tibble::tibble())

    model_to_use <- preferred[preferred %in% unique(as.character(dat$model))][[1]]
    dat %>% filter(.data$model == model_to_use)
  })
}

make_publication_figures <- function(window_results,
                                     out_dir,
                                     min_k = 2) {
  if (!all(c("session", "follow_up") %in% names(window_results))) {
    message("Skipping publication figures: both 'session' and 'follow_up' windows are required.")
    return(invisible(NULL))
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  preds_session <- select_model_preds(window_results$session$dr_molecule, preferred = "spline")
  preds_follow_up_all <- window_results$follow_up$dr_molecule$preds
  preds_follow_up <- tibble::tibble()

  if (!is.null(preds_follow_up_all) && nrow(preds_follow_up_all)) {
    preds_follow_up <- select_model_preds_by_molecule(
      preds = preds_follow_up_all,
      molecules = c("LSD", "MDMA"),
      preferred = c("spline", "linear")
    )
  }

  message("-> Publication figure 1 ...")
  plot_fig1_global_dr_session_followup(
    preds_session = preds_session,
    preds_followup = preds_follow_up,
    outfile = file.path(out_dir, "Fig1_Global_DR_session_followup.pdf")
  )

  message("-> Publication figure 2 ...")
  plot_fig2_dr_top_aes_spline_facets(
    preds_ae = window_results$session$dr_ae$preds,
    models_ae = window_results$session$dr_ae$models,
    outfile = file.path(out_dir, "Fig2_DR_sharedAEs_spline_df3.pdf"),
    min_k = max(3, min_k),
    top_n_ae = Inf,
    significant_only = FALSE
  )

  invisible(list(
    fig1 = file.path(out_dir, "Fig1_Global_DR_session_followup.pdf"),
    fig2 = file.path(out_dir, "Fig2_DR_sharedAEs_spline_df3.pdf")
  ))
}

run_main_analysis <- function(
    data_xlsx = here::here("data", "Adverse-events-dose-v5.xlsx"),
    sheet = "Feuil1",
    out_dir = here::here("results", "main"),
    min_k = 2,
    fit_spline = TRUE,
    make_paper_tables = TRUE,
    paper_dir = here::here("results", "paper_tables"),
    figures_dir = here::here("results", "figures_paper"),
    include_non_primary_windows = FALSE,
    clean_legacy_main_outputs = TRUE) {
  if (!file.exists(data_xlsx)) {
    stop("Data file not found: ", data_xlsx)
  }

  if (isTRUE(clean_legacy_main_outputs)) {
    legacy_dirs <- file.path(out_dir, c("compare", "session", "follow_up"))
    existing_legacy_dirs <- legacy_dirs[dir.exists(legacy_dirs)]
    if (length(existing_legacy_dirs)) {
      message("-> Cleaning legacy main-analysis outputs ...")
      unlink(existing_legacy_dirs, recursive = TRUE, force = TRUE)
    }
  }

  dir.create(paper_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  unlink(file.path(figures_dir, "*.pdf"), force = TRUE)

  message("1) Load & harmonize ...")
  raw <- load_data(data_xlsx, sheet = sheet)
  if (!"time_window" %in% names(raw)) {
    stop("Loaded data does not contain a 'time_window' column.")
  }

  raw <- raw %>%
    mutate(
      time_window = norm_window_label(time_window),
      time_window = tolower(trimws(time_window)),
      time_window = dplyr::case_when(
        time_window %in% c("followup", "follow-up", "follow up", "follow_up") ~ "follow_up",
        time_window %in% c("session", "acute", "in-session", "in session") ~ "session",
        TRUE ~ time_window
      )
    )

  available_windows <- unique(raw$time_window)
  if (any(available_windows %in% c("followup", "follow-up", "follow up"))) {
    stop("Non-normalized time_window values remain: ", paste(available_windows, collapse = ", "))
  }
  if (!length(available_windows)) {
    stop("No time windows found in the dataset.")
  }

  window_order <- c("session", "follow_up")
  if (isTRUE(include_non_primary_windows)) {
    ordered_windows <- unique(c(intersect(window_order, available_windows), setdiff(available_windows, window_order)))
  } else {
    ordered_windows <- intersect(window_order, available_windows)
  }

  run_window <- function(window_value) {
    paper_dir_window <- file.path(paper_dir, window_value)
    model_dir_window <- file.path(paper_dir_window, "model_inputs")
    dir.create(model_dir_window, recursive = TRUE, showWarnings = FALSE)

    message(sprintf("-> Window '%s': build contrasts ...", window_value))
    es <- build_es_for_window(raw, window_value, ref_policies = default_ref_policies)
    if (is.null(es) || !nrow(es)) {
      warning("No effect sizes constructed for window '", window_value, "'.")
      return(NULL)
    }

    message(sprintf("-> Window '%s': dose-response models ...", window_value))
    dr_mol_linear <- run_dr_by_molecule(es, min_k = min_k, model = "linear", grid = "continuous", n_grid = 160)
    dr_ae_linear  <- run_dr_by_ae(es, min_k = min_k, model = "linear", grid = "continuous", n_grid = 160)

    if (isTRUE(fit_spline)) {
      dr_mol_spline <- run_dr_by_molecule(es, min_k = min_k, model = "spline", df_spline = 3, grid = "continuous", n_grid = 160)
      dr_ae_spline  <- run_dr_by_ae(es, min_k = min_k, model = "spline", df_spline = 3, grid = "continuous", n_grid = 160)

      dr_mol <- list(
        preds = bind_rows(dr_mol_linear$preds, dr_mol_spline$preds),
        models = bind_rows(dr_mol_linear$models, dr_mol_spline$models)
      )
      dr_ae <- list(
        preds = bind_rows(dr_ae_linear$preds, dr_ae_spline$preds),
        models = bind_rows(dr_ae_linear$models, dr_ae_spline$models)
      )
    } else {
      dr_mol <- dr_mol_linear
      dr_ae <- dr_ae_linear
    }

    message(sprintf("-> Window '%s': publication tables ...", window_value))
    sig_mol <- save_significance_table(
      dr_mol$models,
      file.path(model_dir_window, "significance_by_molecule_models.csv")
    )
    sig_ae <- save_significance_table(
      dr_ae$models,
      file.path(model_dir_window, "significance_by_ae_models.csv")
    )
    agg_mol <- save_agg_significance_by_molecule(
      dr_mol$models,
      file.path(model_dir_window, "significance_agg_by_molecule.csv")
    )
    agg_ae <- save_agg_significance_by_ae_molecule(
      dr_ae$models,
      file.path(model_dir_window, "significance_agg_by_ae_molecule.csv")
    )

    if (isTRUE(make_paper_tables) && window_value %in% c("session", "follow_up")) {
      make_all_paper_tables(
        path_mol_agg = file.path(model_dir_window, "significance_agg_by_molecule.csv"),
        path_ae_mol_agg = file.path(model_dir_window, "significance_agg_by_ae_molecule.csv"),
        path_models_ae = file.path(model_dir_window, "significance_by_ae_models.csv"),
        output_dir = paper_dir_window
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
      paper_dir = paper_dir_window
    )
  }

  message("-> Supplementary absolute-rate tables ...")
  absolute_rate_outputs <- make_absolute_rate_tables(
    raw = raw,
    out_dir = file.path(paper_dir, "absolute_rates"),
    min_total_n = 1,
    top_n_ae_per_group = 10,
    write_outputs = TRUE
  )

  window_results <- purrr::map(ordered_windows, run_window)
  names(window_results) <- ordered_windows
  window_results <- purrr::compact(window_results)

  figure_outputs <- make_publication_figures(
    window_results = window_results,
    out_dir = figures_dir,
    min_k = min_k
  )

  invisible(list(
    raw = raw,
    windows = window_results,
    figures = figure_outputs,
    absolute_rates = absolute_rate_outputs
  ))
}

if (identical(environment(), globalenv())) {
  run_main_analysis()
}
