suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(gt)
})

fmt_p <- function(p){
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "<0.001", formatC(p, format="f", digits=3)))
}

sig_flag <- function(p, alpha = 0.05){
  !is.na(p) & p < alpha
}

load_robustness_4 <- function(root = here::here("results","main")){
  paths <- tibble::tribble(
    ~window,      ~type,        ~path,
    "session",    "molecule",    file.path(root, "session",   "robustness", "robustness_linear_vs_spline_by_molecule.csv"),
    "session",    "ae",          file.path(root, "session",   "robustness", "robustness_linear_vs_spline_by_ae_molecule.csv"),
    "follow_up",  "molecule",    file.path(root, "follow_up", "robustness", "robustness_linear_vs_spline_by_molecule.csv"),
    "follow_up",  "ae",          file.path(root, "follow_up", "robustness", "robustness_linear_vs_spline_by_ae_molecule.csv")
  )
  
  missing <- paths %>% filter(!file.exists(path))
  if (nrow(missing)) {
    stop("Missing robustness files:\n", paste(missing$path, collapse = "\n"))
  }
  
  mol <- paths %>% filter(type=="molecule") %>%
    mutate(dat = map(path, ~ read_csv(.x, show_col_types = FALSE))) %>%
    select(window, dat) %>% tidyr::unnest(dat)
  
  ae <- paths %>% filter(type=="ae") %>%
    mutate(dat = map(path, ~ read_csv(.x, show_col_types = FALSE))) %>%
    select(window, dat) %>% tidyr::unnest(dat)
  
  list(molecule = mol, ae = ae)
}

make_table_S1_molecule_robustness <- function(mol_df, outdir){
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  tab <- mol_df %>%
    transmute(
      Window = window,
      Molecule = molecule,
      `β (linear)` = round(beta_linear, 4),
      `p (linear)` = fmt_p(p_linear),
      `p (spline omnibus)` = fmt_p(p_spline),
      `k` = k_linear,
      Robustness = robustness
    ) %>%
    arrange(Window, Molecule)
  
  gt_tbl <- tab %>%
    gt() %>%
    tab_header(title = md("**Robustness of global dose–response (linear vs spline)**")) %>%
    cols_align("center", columns = everything())
  
  gtsave(gt_tbl, file.path(outdir, "S1_robustness_global_by_molecule.html"))
  writeLines(as.character(gt_tbl), con = file.path(outdir, "S1_robustness_global_by_molecule.txt"))
  readr::write_csv(tab, file.path(outdir, "S1_robustness_global_by_molecule.csv"))
  
  invisible(tab)
}

make_table_S2_counts_ae_robustness <- function(ae_df, outdir){
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  counts <- ae_df %>%
    group_by(window, molecule, robustness) %>%
    summarise(n = n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = robustness, values_from = n, values_fill = 0) %>%
    mutate(total = rowSums(across(where(is.numeric)))) %>%
    arrange(window, molecule)
  
  gt_tbl <- counts %>%
    rename(Window = window, Molecule = molecule) %>%
    gt() %>%
    tab_header(title = md("**AE-level robustness counts (linear vs spline)**"),
               subtitle = md("_Counts of AE×molecule tests by robustness category; full results provided as Supplementary Data (CSV)._")) %>%
    cols_align("center", columns = -c(Window, Molecule))
  
  gtsave(gt_tbl, file.path(outdir, "S2_robustness_counts_AE_by_molecule.html"))
  readr::write_csv(counts, file.path(outdir, "S2_robustness_counts_AE_by_molecule.csv"))
  
  invisible(counts)
}

make_topline_results_text <- function(mol_df, ae_df, alpha = 0.05){
  # global
  global <- mol_df %>%
    group_by(window) %>%
    summarise(
      n_molecules = n(),
      n_both_sig = sum(robustness == "both_significant", na.rm = TRUE),
      .groups = "drop"
    )
  
  # AE-level counts of any signal (linear or spline)
  ae_any <- ae_df %>%
    mutate(any_sig = sig_flag(p_linear, alpha) | sig_flag(p_spline, alpha)) %>%
    group_by(window, molecule) %>%
    summarise(
      n_tests = n(),
      n_any_sig = sum(any_sig, na.rm = TRUE),
      n_both_sig = sum(robustness == "both_significant", na.rm = TRUE),
      n_linear_only = sum(robustness == "linear_only", na.rm = TRUE),
      n_spline_only = sum(robustness == "spline_only", na.rm = TRUE),
      .groups = "drop"
    )
  
  list(global = global, ae_any = ae_any)
}

run_robustness_reporting <- function(
    root_results = here::here("results","main"),
    outdir = here::here("results","paper_tables","robustness_summary")
){
  dat <- load_robustness_4(root = root_results)
  mol <- dat$molecule
  ae  <- dat$ae
  
  S1 <- make_table_S1_molecule_robustness(mol, outdir)
  S2 <- make_table_S2_counts_ae_robustness(ae, outdir)
  
  txt <- make_topline_results_text(mol, ae)
  
  # write short “Results-ready” text snippets
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(
      "ROBUSTNESS RESULTS (AUTO-SUMMARY)",
      "",
      paste0("Global (by molecule):"),
      paste0(apply(txt$global, 1, function(r){
        sprintf("Window %s: %s molecules analysed; %s showed concordant significance across linear and spline models.",
                r[["window"]], r[["n_molecules"]], r[["n_both_sig"]])
      }), collapse = "\n"),
      "",
      "AE-level (any signal = p_linear<0.05 OR p_spline<0.05):",
      paste0(apply(txt$ae_any, 1, function(r){
        sprintf("Window %s, %s: %s/%s AE×molecule tests showed any signal (both=%s; linear-only=%s; spline-only=%s).",
                r[["window"]], r[["molecule"]], r[["n_any_sig"]], r[["n_tests"]],
                r[["n_both_sig"]], r[["n_linear_only"]], r[["n_spline_only"]])
      }), collapse = "\n")
    ),
    file.path(outdir, "robustness_results_autotext.txt")
  )
  
  invisible(list(S1=S1, S2=S2, txt=txt))
}
