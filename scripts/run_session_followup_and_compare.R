# /Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/scripts/run_session_followup_and_compare.R

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readxl)
  library(purrr)
})

# ---- Source project scripts (absolute paths) ----
source("/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/scripts/utils_data.R")
source("/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/scripts/analysis_dose_response.R")
source("/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/scripts/analysis_forest_plots.R")
source("/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/scripts/analysis_plots_dr.R")
source("/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/scripts/analysis_tables.R")
# comparison helpers
source("/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/scripts/forest_compare_session_followup.R")
source("/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/scripts/dr_compare_session_followup.R")

# ---- Reference policies ----
default_ref_policies <- list(
  MDMA       = c("inactive_placebo","active_non_psy_placebo","active_placebo"),
  LSD        = c("inactive_placebo","active_placebo","active_non_psy_placebo"),
  PSILOCYBIN = c("inactive_placebo","active_non_psy_placebo","active_placebo"),
  AYAHUASCA  = c("inactive_placebo","active_placebo","active_non_psy_placebo"),
  .default   = c("inactive_placebo","active_placebo","active_non_psy_placebo")
)

# ---- Small helpers ----
diagnose_df <- function(df, label="df"){
  cat("----", label, "----\n")
  cat("Names:", paste(names(df), collapse=", "), "\n")
  cat("Rows:", nrow(df), "\n")
  if ("time_window" %in% names(df)) print(table(df$time_window, useNA="ifany"))
  if (all(c("n","events") %in% names(df))) {
    bad <- sum(df$n < df$events, na.rm=TRUE)
    cat("Rows with n<events:", bad, "\n")
  }
  cat("Preview:\n"); print(utils::head(df, 4))
  cat("--------------\n")
}

# Map your input (v5 or raw) to the *legacy schema required by build_pairwise_2x2*:
# REQUIRED (legacy): study_id, molecule, ae_term, time_window, group, arm_type, dose_mg, events, n
map_to_legacy_schema <- function(df){
  df <- as.data.frame(df)
  
  # 1) Build/rename counts to legacy names: n, events
  if (!"n" %in% names(df)) {
    if ("n_total" %in% names(df)) {
      df <- dplyr::rename(df, n = n_total)
    } else if ("n_participants_arm" %in% names(df)) {
      df <- dplyr::rename(df, n = n_participants_arm)
    }
  }
  if (!"events" %in% names(df) && "n_events" %in% names(df)) {
    df <- dplyr::rename(df, events = n_events)
  }
  
  # 2) Create 'group' from arm_id (preferred) or arm_type if present
  if (!"group" %in% names(df)) {
    if ("arm_id" %in% names(df)) {
      df$group <- as.character(df$arm_id)
    } else if ("arm_type" %in% names(df)) {
      df$group <- as.character(df$arm_type)
    } else {
      df$group <- NA_character_
    }
  }
  # Ensure 'arm_type' exists (some downstream code inspects it)
  if (!"arm_type" %in% names(df)) {
    df$arm_type <- df$group
  }
  
  # 3) Normalize, coerce types (but DO NOT coerce events to integer yet)
  df <- df %>%
    dplyr::mutate(
      study_id    = as.character(study_id),
      molecule    = as.character(molecule),
      ae_term     = as.character(ae_term),
      time_window = tolower(trimws(gsub("\u00A0"," ", as.character(time_window), fixed = TRUE))),
      time_window = stringr::str_replace_all(time_window, "[-\\s]+", "_"),
      time_window = stringr::str_replace_all(time_window, "_+", "_"),
      group       = as.character(group),
      arm_type    = as.character(arm_type),
      dose_mg     = suppressWarnings(as.numeric(dose_mg)),
      n           = suppressWarnings(as.numeric(n)),
      events      = suppressWarnings(as.numeric(events))
    )
  
  # 4) Detect if 'events' looks like a proportion (0–1) and convert to counts
  #    Rule: if max(events, na.rm=T) <= 1 AND some events have decimals, treat as proportion
  looks_like_prop <- is.finite(suppressWarnings(max(df$events, na.rm = TRUE))) &&
    suppressWarnings(max(df$events, na.rm = TRUE)) <= 1 &&
    any(abs(df$events - round(df$events)) > .Machine$double.eps^0.5, na.rm = TRUE)
  
  if (looks_like_prop) {
    df <- df %>%
      dplyr::mutate(
        events = round(pmax(0, pmin(1, events)) * n)
      )
  }
  
  # 5) Final coercion to integers and basic validity filters
  df <- df %>%
    dplyr::mutate(
      n      = as.integer(round(n)),
      events = as.integer(round(events))
    ) %>%
    dplyr::filter(
      !is.na(study_id), !is.na(molecule), !is.na(ae_term),
      !is.na(time_window), !is.na(group),
      is.finite(dose_mg),
      !is.na(n), !is.na(events),
      n >= events, n > 0
    )
  
  # 6) Final column presence check for build_pairwise_2x2
  required <- c("study_id","molecule","ae_term","time_window","group","arm_type","dose_mg","events","n")
  miss <- setdiff(required, names(df))
  if (length(miss)) {
    stop("map_to_legacy_schema(): missing columns after mapping: ", paste(miss, collapse=", "))
  }
  df
}
# ---- Main runner (compare session vs follow_up without touching build_pairwise_2x2) ----
run_session_followup_and_compare <- function(
    base_xlsx        = "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/data/Adverse-events-dose-v5.xlsx",
    sheet            = "Feuil1",
    out_root_session = "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/results_session",
    out_root_followup= "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/results_followup",
    out_root_compare = "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/results_compare",
    min_k            = 2,
    fit_spline       = TRUE,
    ref_policies     = default_ref_policies
){
  message("== Reading Excel …")
  stopifnot(file.exists(base_xlsx))
  d_raw <- readxl::read_excel(base_xlsx, sheet = sheet)
  message("Columns in Excel: ", paste(names(d_raw), collapse=", "))
  
  message("== Mapping to legacy schema required by build_pairwise_2x2() …")
  d_legacy <- map_to_legacy_schema(d_raw)
  diagnose_df(d_legacy, "d_legacy (post-map)")
  
  # split windows
  df_session  <- d_legacy %>% filter(time_window == "session")
  df_followup <- d_legacy %>% filter(time_window == "follow_up")
  
  message("Rows session:  ", nrow(df_session))
  message("Rows follow_up:", nrow(df_followup))
  
  # helper to run full pipeline for one window
  .run_one_window <- function(df_in, out_root){
    dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
    message("→ Pairwise contrasts …")
    contr <- build_pairwise_2x2(df_in, ref_policies = ref_policies)
    
    message("→ Escalc …")
    es <- build_escalc(contr)
    
    message("→ Dose–response by molecule …")
    dr_mol <- run_dr_by_molecule(es, min_k = min_k, fit_spline = fit_spline)
    
    message("→ Dose–response by AE …")
    dr_ae  <- run_dr_by_ae(es, min_k = min_k, fit_spline = fit_spline)
    
    # tables + a couple of figures
    dir.create(file.path(out_root,"tables"), recursive = TRUE, showWarnings = FALSE)
    save_significance_table(dr_mol$models, file.path(out_root,"tables/significance_by_molecule_models.csv"))
    save_significance_table(dr_ae$models,  file.path(out_root,"tables/significance_by_ae_models.csv"))
    save_agg_significance_by_molecule(dr_mol$models, file.path(out_root,"tables/significance_agg_by_molecule.csv"))
    save_agg_significance_by_ae_molecule(dr_ae$models,  file.path(out_root,"tables/significance_agg_by_ae_molecule.csv"))
    
    make_forest_summary_per_molecule(es, file.path(out_root,"forest_plots_summary"))
    plot_dr_by_molecule_split(dr_mol$preds, dr_mol$models, file.path(out_root,"dose_response/by_molecule_split"))
    
    list(es=es, dr_mol=dr_mol, dr_ae=dr_ae)
  }
  
  res_s <- NULL; res_f <- NULL
  if (nrow(df_session))  { message("== Running SESSION ==");   res_s <- .run_one_window(df_session,  out_root_session) } else message("⚠️ No session rows.")
  if (nrow(df_followup)) { message("== Running FOLLOW-UP =="); res_f <- .run_one_window(df_followup, out_root_followup) } else message("⚠️ No follow_up rows.")
  
  if (is.null(res_s) || is.null(res_f)) {
    message("⚠️ Comparison skipped (need both windows).")
    return(invisible(list(session=res_s, follow_up=res_f)))
  }
  
  # comparative outputs
  es_all <- bind_rows(
    res_s$es %>% mutate(time_window = "session"),
    res_f$es %>% mutate(time_window = "follow_up")
  )
  dir.create(out_root_compare, recursive = TRUE, showWarnings = FALSE)
  
  message("== Comparative FOREST (side-by-side + combined) ==")
  forest_compare_all_molecules(
    es_all,
    outdir = file.path(out_root_compare, "forest_by_ae_side_by_side"),
    min_k_per_window = max(2, min_k),
    order_by = "pval",
    star_offset = 0.06
  )
  forest_compare_all_molecules_combined(
    es_all,
    outfile = file.path(out_root_compare, "forest_combined_all_molecules.pdf"),
    min_k_per_window = max(2, min_k),
    order_by = "pval",
    star_offset = 0.06,
    ncol = 2
  )
  
  message("== Comparative DOSE–RESPONSE (molecule facets) ==")
  dr_compare_all_molecules(
    preds_session  = res_s$dr_mol$preds,
    preds_followup = res_f$dr_mol$preds,
    outfile = file.path(out_root_compare, "dr_by_molecule_session_vs_followup.pdf")
  )
  
  message("✅ Done. Outputs:")
  message(" - ", normalizePath(out_root_session,  mustWork = FALSE))
  message(" - ", normalizePath(out_root_followup, mustWork = FALSE))
  message(" - ", normalizePath(out_root_compare,  mustWork = FALSE))
  
  invisible(list(session=res_s, follow_up=res_f, es_all=es_all))
}