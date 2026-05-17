suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

# Build supplementary absolute-rate tables to contextualize OR-based analyses.
# Required in `raw`: molecule, ae_term, time_window, study_id, n
# Optional in `raw`: absolute_events (preferred), events (fallback)
make_absolute_rate_tables <- function(raw,
                                      out_dir,
                                      min_total_n = 1,
                                      top_n_ae_per_group = 10,
                                      arm_types_keep = "active",
                                      write_outputs = TRUE) {
  needed <- c("molecule", "ae_term", "time_window", "study_id", "n")
  missing_cols <- setdiff(needed, names(raw))
  if (length(missing_cols)) {
    stop(
      "make_absolute_rate_tables() missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  dat <- raw %>%
    mutate(
      molecule = toupper(as.character(molecule)),
      ae_term = str_squish(str_to_lower(as.character(ae_term))),
      time_window = as.character(time_window),
      study_id = as.character(study_id),
      arm_id = if ("arm_id" %in% names(.)) as.character(arm_id) else NA_character_,
      group = if ("group" %in% names(.)) as.character(group) else NA_character_,
      arm_type = if ("arm_type" %in% names(.)) as.character(arm_type) else NA_character_,
      dose_mg = if ("dose_mg" %in% names(.)) suppressWarnings(as.numeric(dose_mg)) else NA_real_,
      events = if ("events" %in% names(.)) suppressWarnings(as.numeric(events)) else NA_real_,
      absolute_events = if ("absolute_events" %in% names(.)) suppressWarnings(as.numeric(absolute_events)) else NA_real_,
      n = suppressWarnings(as.numeric(n)),
      events_for_absolute = dplyr::coalesce(absolute_events, events),
      events_source_row = dplyr::case_when(
        !is.na(absolute_events) ~ "absolute_events",
        !is.na(events) ~ "events",
        TRUE ~ "missing"
      )
    ) %>%
    filter(is.null(arm_types_keep) | is.na(arm_type) | arm_type %in% arm_types_keep) %>%
    filter(
      !is.na(molecule),
      !is.na(ae_term),
      !is.na(time_window),
      !is.na(study_id),
      is.finite(events_for_absolute),
      is.finite(n),
      n > 0
    )

  arm_cols <- c("time_window", "molecule", "study_id", "arm_id", "group", "arm_type", "dose_mg")
  arm_denominators <- dat %>%
    distinct(across(all_of(arm_cols)), n) %>%
    group_by(time_window, molecule) %>%
    summarise(
      n_studies = n_distinct(study_id),
      n_arms = n(),
      n_total = sum(n, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_total >= min_total_n)
  
  # -------------------------------------------------------------------------
  # Base table at arm granularity: study × arm × molecule × window × AE.
  # `n` is deduplicated per arm instead of summed over repeated AE rows.
  # -------------------------------------------------------------------------
  by_ae_arm <- dat %>%
    group_by(time_window, molecule, study_id, arm_id, group, arm_type, dose_mg, ae_term) %>%
    summarise(
      events_total = sum(events_for_absolute, na.rm = TRUE),
      n_total = max(n, na.rm = TRUE),
      rate = ifelse(n_total > 0, events_total / n_total, NA_real_),
      rate_pct = 100 * rate,
      .groups = "drop"
    )

  by_ae_study <- by_ae_arm %>%
    group_by(time_window, molecule, study_id, ae_term) %>%
    summarise(
      n_arms = n(),
      events_total = sum(events_total, na.rm = TRUE),
      n_total = sum(n_total, na.rm = TRUE),
      rate = ifelse(n_total > 0, events_total / n_total, NA_real_),
      rate_pct = 100 * rate,
      .groups = "drop"
    )
  
  # -------------------------------------------------------------------------
  # AE-level pooled table: molecule × window × AE
  # -------------------------------------------------------------------------
  by_ae <- by_ae_study %>%
    group_by(time_window, molecule, ae_term) %>%
    summarise(
      n_studies = n_distinct(study_id),
      events_total = sum(events_total, na.rm = TRUE),
      n_total = sum(n_total, na.rm = TRUE),
      rate = ifelse(n_total > 0, events_total / n_total, NA_real_),
      rate_pct = 100 * rate,
      .groups = "drop"
    ) %>%
    filter(n_total >= min_total_n) %>%
    mutate(events_source = "row-wise preference: absolute_events else events") %>%
    arrange(time_window, molecule, desc(rate_pct), ae_term)
  
  # -------------------------------------------------------------------------
  # Global "any AE" table: molecule × window
  # Priority: explicit any-AE rows if available
  # Fallback: sum across AE rows (flagged as caution)
  # -------------------------------------------------------------------------
  any_ae_labels <- c(
    "any adverse event", "any ae", "overall adverse events", "all adverse events"
  )
  
  dat_any <- by_ae_study %>%
    filter(tolower(trimws(ae_term)) %in% any_ae_labels)
  
  if (nrow(dat_any) > 0) {
    global <- dat_any %>%
      group_by(time_window, molecule) %>%
      summarise(
        n_studies = n_distinct(study_id),
        events_total = sum(events_total, na.rm = TRUE),
        n_total = sum(n_total, na.rm = TRUE),
        rate = ifelse(n_total > 0, events_total / n_total, NA_real_),
        rate_pct = 100 * rate,
        .groups = "drop"
      ) %>%
      filter(n_total >= min_total_n) %>%
      mutate(
        global_definition = "ae_term == any adverse event",
        events_source = "row-wise preference: absolute_events else events"
      ) %>%
      arrange(time_window, molecule)
  } else {
    global <- arm_denominators %>%
      mutate(
        events_total = NA_real_,
        rate = NA_real_,
        rate_pct = NA_real_,
        global_definition = "denominator only: no explicit any-AE rows; AE rows cannot be summed without double-counting participants",
        events_source = "row-wise preference: absolute_events else events"
      ) %>%
      select(time_window, molecule, n_studies, events_total, n_total, rate, rate_pct, global_definition, events_source) %>%
      arrange(time_window, molecule)
  }
  
  # -------------------------------------------------------------------------
  # Top-N AE per molecule × window
  # -------------------------------------------------------------------------
  top_ae <- by_ae %>%
    group_by(time_window, molecule) %>%
    slice_max(order_by = rate_pct, n = top_n_ae_per_group, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(time_window, molecule, desc(rate_pct), ae_term)

  fmt_rate_cell <- function(events, denominator, pct) {
    dplyr::case_when(
      is.na(events) | is.na(denominator) | is.na(pct) ~ "",
      TRUE ~ sprintf(
        "%s/%s (%.1f%%)",
        format(round(events, 1), trim = TRUE, scientific = FALSE),
        format(round(denominator, 0), trim = TRUE, scientific = FALSE),
        pct
      )
    )
  }

  # -------------------------------------------------------------------------
  # Clinician-facing supplementary table.
  # Keeps the molecule-level participant denominator visible, while each AE
  # percentage uses the AE-specific denominator to avoid underestimating events
  # that were not reported by every study.
  # -------------------------------------------------------------------------
  clinical_long <- by_ae %>%
    left_join(
      arm_denominators %>%
        select(
          time_window,
          molecule,
          n_studies_molecule = n_studies,
          n_total_molecule = n_total
        ),
      by = c("time_window", "molecule")
    ) %>%
    transmute(
      time_window,
      molecule,
      n_studies_molecule,
      n_total_molecule,
      ae_term,
      n_studies_ae = n_studies,
      n_total_ae = n_total,
      events_total,
      event_pct = rate_pct,
      events_per_1000 = 1000 * rate,
      display = fmt_rate_cell(events_total, n_total, rate_pct)
    ) %>%
    arrange(time_window, molecule, desc(event_pct), ae_term)

  clinical_wide <- clinical_long %>%
    mutate(
      ae_column = str_replace_all(str_to_lower(ae_term), "[^a-z0-9]+", "_"),
      ae_column = str_replace_all(ae_column, "^_|_$", ""),
      ae_column = paste0("ae_", ae_column)
    ) %>%
    select(
      time_window,
      molecule,
      n_studies_molecule,
      n_total_molecule,
      ae_column,
      display
    ) %>%
    group_by(time_window, molecule, n_studies_molecule, n_total_molecule, ae_column) %>%
    summarise(
      display = paste(unique(display[display != ""]), collapse = "; "),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = ae_column,
      values_from = display,
      values_fill = list(display = "")
    ) %>%
    arrange(time_window, molecule)
  
  if (isTRUE(write_outputs)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    readr::write_csv(global, file.path(out_dir, "Sx_absolute_rates_global.csv"))
    readr::write_csv(by_ae, file.path(out_dir, "Sx_absolute_rates_by_ae.csv"))
    readr::write_csv(top_ae, file.path(out_dir, "Sx_absolute_rates_topAE.csv"))
    readr::write_csv(clinical_long, file.path(out_dir, "Sx_absolute_rates_clinician_summary.csv"))
    readr::write_csv(clinical_wide, file.path(out_dir, "Sx_absolute_rates_clinician_wide.csv"))
    readr::write_csv(by_ae_study, file.path(out_dir, "Sx_absolute_rates_by_ae_study.csv"))
    readr::write_csv(by_ae_arm, file.path(out_dir, "Sx_absolute_rates_by_ae_arm.csv"))
    readr::write_csv(arm_denominators, file.path(out_dir, "Sx_absolute_rates_denominators.csv"))
    
    note <- c(
      "Absolute-rate supplementary tables",
      paste0("- included arm_type values: ", ifelse(is.null(arm_types_keep), "all", paste(arm_types_keep, collapse = ", "))),
      "- events_total uses absolute_events when available, otherwise events",
      "- AE-level table aggregates study × arm × molecule × window × AE first, deduplicating N per arm, then pools",
      "- clinician_summary keeps molecule-level total participants visible and uses AE-specific denominators for each AE percentage",
      "- clinician_wide reports each AE as events/N (%) for quick clinical reading",
      "- global table uses 'any adverse event' rows when available",
      "- if no any-AE rows are available, the global table reports denominator only; summing AE rows would double-count participants"
    )
    writeLines(note, con = file.path(out_dir, "README_absolute_rates.txt"))
  }
  
  list(
    global = global,
    by_ae = by_ae,
    top_ae = top_ae,
    clinical_summary = clinical_long,
    clinical_wide = clinical_wide,
    by_ae_study = by_ae_study,
    by_ae_arm = by_ae_arm,
    denominators = arm_denominators
  )
}
