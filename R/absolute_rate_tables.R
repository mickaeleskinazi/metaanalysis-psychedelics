suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# Build supplementary absolute-rate tables to contextualize OR-based analyses.
# Required columns in `raw`: molecule, ae_term, time_window, study_id, n
# Optional columns: absolute_events (preferred), events (fallback)
make_absolute_rate_tables <- function(raw,
                                      out_dir,
                                      min_total_n = 1,
                                      top_n_ae_per_group = 10,
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
      ae_term = as.character(ae_term),
      time_window = as.character(time_window),
      study_id = as.character(study_id),
      events = if ("events" %in% names(.)) suppressWarnings(as.numeric(events)) else NA_real_,
      absolute_events = if ("absolute_events" %in% names(.)) suppressWarnings(as.numeric(absolute_events)) else NA_real_,
      n = suppressWarnings(as.numeric(n)),
      events_for_absolute = dplyr::coalesce(absolute_events, events),
      events_source = dplyr::case_when(
        !is.na(absolute_events) ~ "absolute_events",
        !is.na(events) ~ "events",
        TRUE ~ "missing"
      )
    ) %>%
    filter(
      !is.na(molecule),
      !is.na(ae_term),
      !is.na(time_window),
      !is.na(study_id),
      is.finite(events_for_absolute),
      is.finite(n),
      n > 0
    )
  
  # ---------------------------------------------------------------------------
  # 1) AE-level absolute rates (base analytique)
  #    On agrège d'abord à l'échelle étude × molécule × fenêtre × AE
  #    pour éviter les duplications de lignes.
  # ---------------------------------------------------------------------------
  by_ae_study <- dat %>%
    group_by(time_window, molecule, study_id, ae_term) %>%
    summarise(
      events_total = sum(events_for_absolute, na.rm = TRUE),
      n_total = sum(n, na.rm = TRUE),
      rate = ifelse(n_total > 0, events_total / n_total, NA_real_),
      rate_pct = 100 * rate,
      .groups = "drop"
    )
  
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
  
  # ---------------------------------------------------------------------------
  # 2) Global "any AE" par molécule × fenêtre
  #    Priorité: utiliser explicitement ae_term == "any adverse event"
  #    pour éviter de sommer des n à travers tous les AE (double-comptage).
  # ---------------------------------------------------------------------------
  any_ae_labels <- c(
    "any adverse event", "any ae", "overall adverse events", "all adverse events"
  )
  
  dat_any <- dat %>%
    filter(tolower(trimws(ae_term)) %in% any_ae_labels)
  
  if (nrow(dat_any) > 0) {
    global <- dat_any %>%
      group_by(time_window, molecule, study_id) %>%
      summarise(
        events_total = sum(events_for_absolute, na.rm = TRUE),
        n_total = sum(n, na.rm = TRUE),
        .groups = "drop"
      ) %>%
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
    # fallback explicite si "any adverse event" absent
    global <- by_ae %>%
      group_by(time_window, molecule) %>%
      summarise(
        n_studies = max(n_studies, na.rm = TRUE),
        events_total = sum(events_total, na.rm = TRUE),
        n_total = sum(n_total, na.rm = TRUE),
        rate = ifelse(n_total > 0, events_total / n_total, NA_real_),
        rate_pct = 100 * rate,
        .groups = "drop"
      ) %>%
      filter(n_total >= min_total_n) %>%
      mutate(
        global_definition = "fallback: sum across AE rows (interpret with caution)",
        events_source = "row-wise preference: absolute_events else events"
      ) %>%
      arrange(time_window, molecule)
  }
  
  # ---------------------------------------------------------------------------
  # 3) Top AE par molécule × fenêtre
  # ---------------------------------------------------------------------------
  top_ae <- by_ae %>%
    group_by(time_window, molecule) %>%
    slice_max(order_by = rate_pct, n = top_n_ae_per_group, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(time_window, molecule, desc(rate_pct), ae_term)
  
  if (isTRUE(write_outputs)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    readr::write_csv(global, file.path(out_dir, "Sx_absolute_rates_global.csv"))
    readr::write_csv(by_ae, file.path(out_dir, "Sx_absolute_rates_by_ae.csv"))
    readr::write_csv(top_ae, file.path(out_dir, "Sx_absolute_rates_topAE.csv"))
    readr::write_csv(by_ae_study, file.path(out_dir, "Sx_absolute_rates_by_ae_study.csv"))
    
    note <- c(
      "Absolute-rate supplementary tables",
      "- events_total uses absolute_events when available, otherwise events",
      "- AE-level tables aggregate at study × molecule × window × AE before pooling",
      "- Global table uses ae_term='any adverse event' when available",
      "- If absent, fallback sums AE rows (possible denominator duplication; interpret cautiously)"
    )
    writeLines(note, con = file.path(out_dir, "README_absolute_rates.txt"))
  }
  
  list(
    global = global,
    by_ae = by_ae,
    top_ae = top_ae,
    by_ae_study = by_ae_study
  )
}