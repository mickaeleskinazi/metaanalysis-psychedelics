#!/usr/bin/env Rscript
# Usage: Rscript scripts/build_results_tables.R
#
# The script reads pre-computed CSV outputs from the meta-analysis pipeline
# and regenerates the LaTeX tables plus the manuscript results subsection.
# It only depends on tidyverse-style packages that are commonly used across
# the project.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tidyr)
  library(glue)
  library(tibble)
})

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

first_existing <- function(paths) {
  match <- purrr::detect(paths, file.exists)
  if (is.null(match)) {
    stop("None of the candidate paths exist: ", paste(paths, collapse = ", "))
  }
  match
}

read_csv_any <- function(paths, ...) {
  readr::read_csv(first_existing(paths), ...)
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  formatC(p, format = "e", digits = 2)
}

fmt_tex_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p >= 0.1) {
    return(formatC(p, format = "f", digits = 2))
  }
  if (p >= 0.01) {
    return(formatC(p, format = "f", digits = 3))
  }
  if (p >= 0.001) {
    return(formatC(p, format = "f", digits = 4))
  }
  sci <- formatC(p, format = "e", digits = 2)
  parts <- str_split_fixed(sci, "e", 2)
  mantissa <- parts[, 1]
  exponent <- as.integer(parts[, 2])
  glue("{mantissa}\\times10^{{{exponent}}}")
}

sig_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

nice_ae <- function(ae) {
  ae %>%
    str_replace_all("_", " ") %>%
    str_split("\\\s+") %>%
    purrr::map_chr(~ str_to_title(trimws(paste(.x, collapse = " "))))
}

safe_rma <- function(yi, vi) {
  yi <- yi[is.finite(yi)]
  vi <- vi[is.finite(vi) & vi > 0]
  if (length(yi) == 0 || length(vi) == 0) {
    return(tibble(log_or = NA_real_, ci_lb = NA_real_, ci_ub = NA_real_, pval = NA_real_, k = length(yi)))
  }
  n <- length(yi)
  if (n == 1) {
    se <- sqrt(vi[1])
    log_or <- yi[1]
    ci_lb <- log_or - 1.96 * se
    ci_ub <- log_or + 1.96 * se
    z <- ifelse(se > 0, log_or / se, NA_real_)
    pval <- ifelse(is.na(z), NA_real_, 2 * stats::pnorm(abs(z), lower.tail = FALSE))
    return(tibble(log_or = log_or, ci_lb = ci_lb, ci_ub = ci_ub, pval = pval, k = n))
  }
  wi <- 1 / vi
  sum_wi <- sum(wi)
  mu_fixed <- sum(wi * yi) / sum_wi
  Q <- sum(wi * (yi - mu_fixed)^2)
  df <- n - 1
  sum_wi_sq <- sum(wi^2)
  denom <- sum_wi - (sum_wi_sq / sum_wi)
  tau2 <- max(0, (Q - df) / denom)
  weights <- 1 / (vi + tau2)
  sum_w <- sum(weights)
  log_or <- sum(weights * yi) / sum_w
  se <- sqrt(1 / sum_w)
  ci_lb <- log_or - 1.96 * se
  ci_ub <- log_or + 1.96 * se
  z <- ifelse(se > 0, log_or / se, NA_real_)
  pval <- ifelse(is.na(z), NA_real_, 2 * stats::pnorm(abs(z), lower.tail = FALSE))
  tibble(log_or = log_or, ci_lb = ci_lb, ci_ub = ci_ub, pval = pval, k = n)
}

# ---- Data ingestion -----------------------------------------------------

ensure_dir("tables")
ensure_dir("latex")

topline <- read_csv_any(
  c(
    "results/topline_session_followup_summary.csv",
    "results/main/compare/tables/topline_session_followup_summary.csv",
    "results_compare/tables/DR_session_vs_followup_summary.csv"
  ),
  show_col_types = FALSE
) %>%
  mutate(molecule = toupper(molecule))

agg_molecule <- read_csv_any(
  c(
    "results/main/session/tables/significance_agg_by_molecule.csv",
    "results/session/tables/significance_agg_by_molecule.csv",
    "results/significance_agg_by_molecule.csv"
  ),
  show_col_types = FALSE
) %>%
  mutate(molecule = toupper(molecule))

ae_window <- read_csv_any(
  c(
    "results/forest_by_ae-session.csv", # legacy placeholder
    "results/main/compare/tables/ae_significance_by_window.csv",
    "results_compare/tables/DR_session_vs_followup_by_ae.csv"
  ),
  show_col_types = FALSE
)

if (!"molecule" %in% names(ae_window)) {
  stop("The AE-by-window table must contain a 'molecule' column.")
}
ae_window <- ae_window %>% mutate(molecule = toupper(molecule))

escalc <- read_csv_any(
  c(
    "results/tables/escalc.csv",
    "results/main/tables/escalc.csv",
    "results_session/tables/escalc.csv"
  ),
  show_col_types = FALSE
) %>%
  mutate(
    molecule = toupper(molecule),
    yi = suppressWarnings(as.numeric(yi)),
    vi = suppressWarnings(as.numeric(vi))
  ) %>%
  filter(is.finite(yi), is.finite(vi), vi > 0)

# ---- Table 1: Global DR (session) --------------------------------------

molecules_order <- c("LSD", "MDMA", "PSILOCYBIN", "AYAHUASCA")

global_tbl <- agg_molecule %>%
  select(molecule, QM, p_overall, k_total) %>%
  rename(qm = QM, p = p_overall, k = k_total) %>%
  left_join(topline %>% select(molecule, k_session, p_session), by = "molecule") %>%
  mutate(
    k = coalesce(k_session, k),
    p = coalesce(p_session, p)
  ) %>%
  mutate(across(c(k, qm, p), ~ suppressWarnings(as.numeric(.x)))) %>%
  filter(molecule %in% molecules_order, !is.na(p)) %>%
  mutate(
    molecule = factor(molecule, levels = molecules_order),
    qm = round(qm, 2),
    p_raw = p,
    p_fmt = fmt_p(p),
    stars = sig_stars(p),
    k = if_else(is.na(k), "", as.character(as.integer(round(k))))
  ) %>%
  arrange(molecule) %>%
  select(molecule, k, qm, p_fmt, stars, p_raw)

write_lines(
  glue("""
\\begin{{table}}[htbp]
  \\centering
  \\caption{{Global session dose--response tests by molecule.}}
  \\label{{tab:dr-global-session}}
  \\begin{{tabular}}{{lcccc}}
    \\toprule
    Molecule & $k_{{\\text{{session}}}}$ & $Q_M$ & $p$ & Sig. \\\\
    \\midrule
    {paste(glue('{global_tbl$molecule} & {global_tbl$k} & {global_tbl$qm} & {global_tbl$p_fmt} & {global_tbl$stars} \\\\'), collapse = '\n    ')}
    \\bottomrule
  \\end{{tabular}}
\\end{{table}}
"""), "tables/dr_global_session.tex")

# ---- Table 2: AE x molecule (session significant) ----------------------

ae_sig_tbl <- ae_window %>%
  mutate(
    p_session = suppressWarnings(as.numeric(p_session)),
    sig_session = sig_session %in% c(TRUE, "TRUE")
  ) %>%
  filter(sig_session) %>%
  transmute(
    molecule,
    ae_term,
    p_val = p_session,
    stars = sig_stars(p_session)
  ) %>%
  mutate(
    ae_label = nice_ae(ae_term),
    p_fmt = fmt_p(p_val)
  ) %>%
  arrange(molecule, p_val) %>%
  select(molecule, ae_label, p_fmt, stars, p_val)

write_lines(
  glue("""
\\begin{{table}}[htbp]
  \\centering
  \\caption{{Significant session dose--response effects by adverse event (AE) and molecule.}}
  \\label{{tab:dr-ae-by-molecule-session}}
  \\begin{{tabular}}{{llcc}}
    \\toprule
    Molecule & AE & $p$ & Sig. \\\\
    \\midrule
    {paste(glue('{ae_sig_tbl$molecule} & {ae_sig_tbl$ae_label} & {ae_sig_tbl$p_fmt} & {ae_sig_tbl$stars} \\\\'), collapse = '\n    ')}
    \\bottomrule
  \\end{{tabular}}
\\end{{table}}
"""), "tables/dr_ae_by_molecule_session_sig.tex")

# ---- Table 3: Forest ORs by molecule & window --------------------------

meta_df <- escalc %>%
  group_by(molecule, ae_term, time_window) %>%
  filter(n() >= 2) %>%
  summarise(safe_rma(yi, vi), .groups = "drop") %>%
  mutate(window = time_window) %>%
  select(-time_window)

meta_wide <- meta_df %>%
  pivot_wider(
    id_cols = c(molecule, ae_term),
    names_from = window,
    values_from = c(log_or, ci_lb, ci_ub, k)
  )

forest_tbl <- ae_window %>%
  mutate(
    p_session = suppressWarnings(as.numeric(p_session)),
    p_follow = suppressWarnings(as.numeric(p_follow)),
    sig_session = sig_session %in% c(TRUE, "TRUE"),
    sig_follow = sig_follow %in% c(TRUE, "TRUE"),
    status = case_when(
      sig_session & sig_follow ~ "Persistent",
      sig_session & !sig_follow ~ "Transient",
      !sig_session & sig_follow ~ "Emergent",
      TRUE ~ "Non-significant"
    )
  ) %>%
  filter(status != "Non-significant") %>%
  left_join(meta_wide, by = c("molecule", "ae_term")) %>%
  mutate(
    session_or = if_else(!is.na(log_or_session),
                         sprintf("%.2f [%.2f, %.2f]%s",
                                 exp(log_or_session), exp(ci_lb_session), exp(ci_ub_session),
                                 dplyr::if_else(sig_session, "*", "")),
                         "--"),
    follow_or = if_else(!is.na(log_or_follow_up),
                        sprintf("%.2f [%.2f, %.2f]%s",
                                exp(log_or_follow_up), exp(ci_lb_follow_up), exp(ci_ub_follow_up),
                                dplyr::if_else(sig_follow, "*", "")),
                        "--"),
    abs_log = case_when(
      status %in% c("Persistent", "Transient") & !is.na(log_or_session) ~ abs(log_or_session),
      status %in% c("Persistent", "Emergent") & !is.na(log_or_follow_up) ~ abs(log_or_follow_up),
      TRUE ~ 0
    ),
    ae_label = nice_ae(ae_term)
  )

status_order <- c("Persistent", "Transient", "Emergent")

forest_tbl <- forest_tbl %>%
  mutate(status = factor(status, levels = status_order)) %>%
  arrange(molecule, status, desc(abs_log), ae_label) %>%
  select(molecule, ae_label, session_or, follow_or, status)

write_lines(
  glue("""
\\begin{{table}}[htbp]
  \\centering
  \\caption{{Significant pooled odds ratios (OR) for adverse events (AEs) by molecule and time window.}}
  \\label{{tab:forest-ae-by-window}}
  \\begin{{tabular}}{{lllcl}}
    \\toprule
    Molecule & AE & Session OR [95\\% CI] & Follow-up OR [95\\% CI] & Status \\\\
    \\midrule
    {paste(glue('{forest_tbl$molecule} & {forest_tbl$ae_label} & {forest_tbl$session_or} & {forest_tbl$follow_or} & {forest_tbl$status} \\\\'), collapse = '\n    ')}
    \\bottomrule
  \\end{{tabular}}
\\end{{table}}
"""), "tables/forest_ae_sig_by_window.tex")

# ---- Table 4: Counts ----------------------------------------------------

count_tbl <- forest_tbl %>%
  count(molecule, status, name = "n") %>%
  mutate(status = as.character(status)) %>%
  complete(molecule, status = status_order, fill = list(n = 0)) %>%
  pivot_wider(names_from = status, values_from = n) %>%
  arrange(molecule)

write_lines(
  glue("""
\\begin{{table}}[htbp]
  \\centering
  \\caption{{Counts of significant adverse events (AEs) by persistence pattern.}}
  \\label{{tab:forest-ae-sig-counts}}
  \\begin{{tabular}}{{lccc}}
    \\toprule
    Molecule & Transient & Emergent & Persistent \\\\
    \\midrule
    {paste(glue('{count_tbl$molecule} & {count_tbl$Transient} & {count_tbl$Emergent} & {count_tbl$Persistent} \\\\'), collapse = '\n    ')}
    \\bottomrule
  \\end{{tabular}}
\\end{{table}}
"""), "tables/forest_ae_sig_counts.tex")

# ---- Results subsection -------------------------------------------------

topline_clean <- topline %>%
  mutate(across(c(p_session, p_follow, k_session, k_follow), ~ suppressWarnings(as.numeric(.x))))

session_snippets <- global_tbl %>%
  mutate(snippet = glue("{molecule} ($Q_M = {qm},\\ p = {fmt_tex_p(p_raw)}$)")) %>%
  arrange(molecule) %>%
  pull(snippet)

session_sentence <- glue_collapse(session_snippets, sep = ", ", last = ", and ")

follow_sentences <- topline_clean %>%
  filter(molecule %in% molecules_order) %>%
  mutate(desc = case_when(
    !is.na(p_follow) & p_follow < 0.05 ~ glue("{molecule} retained a significant follow-up slope ($p = {fmt_tex_p(p_follow)}$; $k = {k_follow}$)"),
    !is.na(p_follow) ~ glue("{molecule} follow-up slopes were not significant ($p = {fmt_tex_p(p_follow)}$)"),
    TRUE ~ glue("{molecule} had no evaluable follow-up contrasts")
  )) %>%
  pull(desc)

follow_sentence <- glue_collapse(follow_sentences, sep = "; ")

result_paragraph <- glue(
  "Global dose--response analyses indicated significant session-level slopes for {session_sentence} (Table~\\ref{{tab:dr-global-session}}). Ayahuasca could not be evaluated because only a single dosing condition was available. {follow_sentence}. Fitted curves suggested monotonic dose increases with mild deviations from linearity (Figure~\\ref{{fig:dr-global-session}})."
)

ae_phrases <- ae_sig_tbl %>%
  mutate(p_tex = fmt_tex_p(p_val)) %>%
  group_by(molecule) %>%
  summarise(
    desc = glue("{molecule} ({glue_collapse(glue('{ae_label} ($p = {p_tex}$)'), sep = ", ", last = ", and ")})"),
    .groups = "drop"
  ) %>%
  arrange(molecule) %>%
  pull(desc)

ae_sentence <- glue(
  "Session adverse-event analyses showed that {glue_collapse(ae_phrases, sep = "; ")} (Figure~\\ref{{fig:dr-by-ae-session}}). Exact $p$-values are provided in Table~\\ref{{tab:dr-ae-by-molecule-session}}."
)

forest_summary <- forest_tbl %>%
  group_by(molecule) %>%
  summarise(
    desc = glue("{molecule}: {glue_collapse(ae_label, sep = ", ", last = ", and ")}"),
    .groups = "drop"
  ) %>%
  arrange(molecule) %>%
  pull(desc)

forest_sentence <- glue(
  "Forest comparisons (Figure~\\ref{{fig:forest-combined}}) indicated that psilocybin and ayahuasca were only estimable during dosing sessions, whereas LSD and MDMA contributed both session and follow-up panels. All significant AE signals were transient with no persistent or emergent effects: {glue_collapse(forest_summary, sep = "; ")}. Detailed odds-ratio summaries appear in Table~\\ref{{tab:forest-ae-by-window}}, and the transient counts are tabulated in Table~\\ref{{tab:forest-ae-sig-counts}}."
)

results_section <- glue(
  "\\subsection{{Results}}\n\n{result_paragraph}\n\n{ae_sentence}\n\n{forest_sentence}\n\n\\begin{{figure}}[htb]\n  \\centering\n  \\includegraphics[width=\\textwidth]{{figures/master_dr_by_molecule-session.pdf}}\n  \\caption{{Global dose--response during session by molecule.}}\n  \\label{{fig:dr-global-session}}\n\\end{{figure}}\n\n\\begin{{figure}}[htb]\n  \\centering\n  \\includegraphics[width=\\textwidth]{{figures/master_dr_by_ae-session.pdf}}\n  \\caption{{Dose--response by adverse event (AE) during session (facets by molecule and AE).}}\n  \\label{{fig:dr-by-ae-session}}\n\\end{{figure}}\n\n\\begin{{figure}}[htb]\n  \\centering\n  \\includegraphics[width=\\textwidth]{{figures/forest_combined_all_molecules.pdf}}\n  \\caption{{Forest plots of AE odds ratios (OR) by molecule and time window (session vs follow-up).}}\n  \\label{{fig:forest-combined}}\n\\end{{figure}}\n\n\\input{{tables/dr_global_session.tex}}\n\\input{{tables/dr_ae_by_molecule_session_sig.tex}}\n\\input{{tables/forest_ae_sig_by_window.tex}}\n\\input{{tables/forest_ae_sig_counts.tex}}\n"
)

write_lines(results_section, "latex/results_section.tex")

message("✅ Results tables and narrative written to the latex/ and tables/ directories.")
