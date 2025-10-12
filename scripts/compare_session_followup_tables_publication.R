#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(gt)
})

source(here::here("scripts", "compare_session_followup_from_saved_tables.R"))

make_publication_tables <- function(
    dir_session = here::here("results_session", "tables"),
    dir_followup = here::here("results_followup", "tables"),
    out_dir = here::here("results_compare", "tables")) {
  res <- compare_session_followup_from_saved_tables(dir_session, dir_followup, out_dir)
  topline <- res$topline

  if (is.null(topline) || !nrow(topline)) {
    warning("No topline results to format.")
    return(invisible(NULL))
  }

  formatted <- topline %>%
    mutate(
      `Session p` = p_session_fmt,
      `Session sig` = stars_session,
      `# significant AE (session)` = n_sig_session,
      `Follow-up p` = p_follow_fmt,
      `Follow-up sig` = stars_follow,
      `# significant AE (follow-up)` = n_sig_follow
    ) %>%
    select(Molecule = molecule, `Session p`, `Session sig`, `# significant AE (session)`,
           `Follow-up p`, `Follow-up sig`, `# significant AE (follow-up)`)

  gt_tbl <- formatted %>%
    gt() %>%
    tab_header(title = md("**Session vs follow-up dose–response summary**")) %>%
    cols_align("center", columns = -Molecule) %>%
    tab_source_note(md("Stars indicate significance thresholds: * p<0.05, ** p<0.01, *** p<0.001."))

  html_path <- file.path(out_dir, "compare_topline_molecule.html")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  gtsave(gt_tbl, html_path)

  invisible(formatted)
}

if (!interactive() && sys.nframe() == 0) {
  make_publication_tables()
}
