suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# Build supplementary absolute-rate tables to contextualize OR-based analyses.
# raw must contain at least: molecule, ae_term, time_window, events, n
# optional: absolute_events (preferred for clinician-facing absolute tables)
make_absolute_rate_tables <- function(raw,
                                      out_dir,
                                      min_total_n = 1,
                                      top_n_ae_per_group = 10,
                                      write_outputs = TRUE) {
  needed <- c("molecule", "ae_term", "time_window", "n")
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
    filter(!is.na(molecule), !is.na(ae_term), !is.na(time_window), is.finite(events_for_absolute), is.finite(n), n > 0)

  # Global absolute rates: "any AE" aggregate over all AE rows.
  global <- dat %>%
    group_by(time_window, molecule) %>%
    summarise(
      events_total = sum(events_for_absolute, na.rm = TRUE),
      n_total = sum(n, na.rm = TRUE),
      rate = events_total / n_total,
      rate_pct = 100 * rate,
      .groups = "drop"
    ) %>%
    filter(n_total >= min_total_n) %>%
    mutate(events_source = "row-wise preference: absolute_events else events") %>%
    arrange(time_window, molecule)

  # AE-level absolute rates by molecule and window.
  by_ae <- dat %>%
    group_by(time_window, molecule, ae_term) %>%
    summarise(
      events_total = sum(events_for_absolute, na.rm = TRUE),
      n_total = sum(n, na.rm = TRUE),
      rate = events_total / n_total,
      rate_pct = 100 * rate,
      .groups = "drop"
    ) %>%
    filter(n_total >= min_total_n) %>%
    mutate(events_source = "row-wise preference: absolute_events else events") %>%
    arrange(time_window, molecule, desc(rate_pct), ae_term)

  # Compact supplementary table: top N AEs per molecule×window.
  top_ae <- by_ae %>%
    group_by(time_window, molecule) %>%
    slice_max(order_by = rate_pct, n = top_n_ae_per_group, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(events_source = "row-wise preference: absolute_events else events") %>%
    arrange(time_window, molecule, desc(rate_pct), ae_term)

  if (isTRUE(write_outputs)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    readr::write_csv(global, file.path(out_dir, "Sx_absolute_rates_global.csv"))
    readr::write_csv(by_ae, file.path(out_dir, "Sx_absolute_rates_by_ae.csv"))
    readr::write_csv(top_ae, file.path(out_dir, "Sx_absolute_rates_topAE.csv"))

    note <- c(
      "Absolute-rate supplementary tables",
      "- rates are descriptive: events_total / n_total (events_total uses absolute_events when available, otherwise events)",
      "- values are aggregated at arm level from extracted trial data",
      "- these tables complement OR-based meta-analytic outputs"
    )
    writeLines(note, con = file.path(out_dir, "README_absolute_rates.txt"))
  }

  list(global = global, by_ae = by_ae, top_ae = top_ae)
}
