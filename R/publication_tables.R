suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(gt)
  library(knitr)
  library(kableExtra)
  library(flextable)
  library(officer)
})

sig_stars_vec <- function(p){
  out <- rep("", length(p))
  out[!is.na(p) & p < 0.001] <- "***"
  out[!is.na(p) & p >= 0.001 & p < 0.01] <- "**"
  out[!is.na(p) & p >= 0.01  & p < 0.05] <- "*"
  out
}

fmt_p <- function(p){
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "<0.001",
                formatC(p, format = "f", digits = 3)))
}

fmt_num <- function(x, d = 2){
  ifelse(is.na(x), "", formatC(x, format = "f", digits = d))
}

add_or_cols <- function(df, est = "estimate", lo = "ci_low", hi = "ci_high"){
  if (all(c(est, lo, hi) %in% names(df))){
    df %>%
      mutate(OR   = exp(.data[[est]]),
             ORlo = exp(.data[[lo]]),
             ORhi = exp(.data[[hi]]))
  } else df
}

save_as_html <- function(gt_tbl, outfile_html){
  dir.create(dirname(outfile_html), recursive = TRUE, showWarnings = FALSE)
  gtsave(gt_tbl, outfile_html)
}

save_as_tex <- function(df, caption, outfile_tex, col_names = NULL, align = NULL){
  dir.create(dirname(outfile_tex), recursive = TRUE, showWarnings = FALSE)
  tex <- kable(df, format = "latex", booktabs = TRUE, escape = TRUE,
               caption = caption, col.names = col_names, align = align) %>%
    kable_styling(full_width = FALSE, position = "center", latex_options = c("hold_position"))
  writeLines(tex, outfile_tex)
}

save_as_docx <- function(df, caption, outfile_docx){
  dir.create(dirname(outfile_docx), recursive = TRUE, showWarnings = FALSE)
  ft <- flextable(df)
  ft <- set_caption(ft, caption)
  doc <- read_docx()
  doc <- body_add_flextable(doc, ft)
  print(doc, target = outfile_docx)
}

make_table_dr_global_by_molecule <- function(path, output_dir){
  stopifnot(file.exists(path))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  x <- read_csv(path, show_col_types = FALSE) %>%
    mutate(
      p = coalesce(p_overall, QMp),
      Stars = sig_stars_vec(p),
      `k (total)` = k_total,
      `Q (model)` = fmt_num(QM, 2),
      `p (overall)` = fmt_p(p)
    ) %>%
    select(Molecule = molecule, `k (total)`, `Q (model)`, `p (overall)`, Stars)

  gt_tbl <- x %>%
    gt() %>%
    tab_header(
      title = md("**Global dose–response per molecule (all AE pooled)**")
    ) %>%
    cols_align("center", columns = everything())

  save_as_html(gt_tbl, file.path(output_dir, "01_dr_global_by_molecule.html"))
  save_as_tex(x,
              "Global dose–response per molecule (all AE pooled)",
              file.path(output_dir, "01_dr_global_by_molecule.tex"))
  save_as_docx(x,
               "Global dose–response per molecule (all AE pooled)",
               file.path(output_dir, "01_dr_global_by_molecule.docx"))
  invisible(x)
}

make_table_matrix_ae_by_molecule <- function(path, output_dir, min_k = 2){
  stopifnot(file.exists(path))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  x <- read_csv(path, show_col_types = FALSE) %>%
    mutate(
      p = coalesce(p_overall, QMp),
      star = sig_stars_vec(p),
      k_total = suppressWarnings(as.integer(k_total))
    ) %>%
    filter(is.na(min_k) | k_total >= min_k)

  mat <- x %>%
    transmute(ae_term, molecule,
              value = ifelse(is.na(p), "", paste0(star, " (p=", fmt_p(p), ")"))) %>%
    distinct() %>%
    pivot_wider(names_from = molecule, values_from = value) %>%
    arrange(ae_term)

  gt_tbl <- mat %>%
    gt() %>%
    tab_header(
      title = md("**Adverse events × molecules (significance)**"),
      subtitle = md("_Stars show overall significance for each AE within molecule; blank = NS_")
    ) %>%
    cols_align("center", columns = -ae_term) %>%
    cols_label(ae_term = "Adverse event")

  save_as_html(gt_tbl, file.path(output_dir, "02_dr_by_ae_molecule_matrix.html"))
  save_as_tex(mat %>% rename(`Adverse event` = ae_term),
              "Adverse events × molecules (significance)",
              file.path(output_dir, "02_dr_by_ae_molecule_matrix.tex"))
  save_as_docx(mat %>% rename(`Adverse event` = ae_term),
               "Adverse events × molecules (significance)",
               file.path(output_dir, "02_dr_by_ae_molecule_matrix.docx"))
  invisible(mat)
}

make_table_top_ae_per_molecule <- function(path, output_dir, top_n = 10, min_k = 2){
  stopifnot(file.exists(path))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  x <- read_csv(path, show_col_types = FALSE) %>%
    mutate(
      p = coalesce(p_overall, QMp),
      star = sig_stars_vec(p),
      k_total = suppressWarnings(as.integer(k_total))
    ) %>%
    filter(is.na(min_k) | k_total >= min_k)

  top <- x %>%
    arrange(molecule, p) %>%
    group_by(molecule) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    transmute(
      Molecule = molecule,
      `Adverse event` = ae_term,
      `k` = k_total,
      `p (overall)` = fmt_p(p),
      Stars = star
    )

  gt_tbl <- top %>%
    gt() %>%
    tab_header(title = md("**Top adverse events per molecule (ordered by p)**")) %>%
    cols_align("center", columns = everything())

  save_as_html(gt_tbl, file.path(output_dir, "03_topAE_per_molecule.html"))
  save_as_tex(top,
              "Top adverse events per molecule (ordered by p)",
              file.path(output_dir, "03_topAE_per_molecule.tex"))
  save_as_docx(top,
               "Top adverse events per molecule (ordered by p)",
               file.path(output_dir, "03_topAE_per_molecule.docx"))
  invisible(top)
}

make_appendix_model_terms <- function(path_models, outfile_csv){
  stopifnot(file.exists(path_models))
  dir.create(dirname(outfile_csv), recursive = TRUE, showWarnings = FALSE)

  x <- read_csv(path_models, show_col_types = FALSE) %>%
    add_or_cols(est = "estimate", lo = "ci_low", hi = "ci_high") %>%
    mutate(
      `β (log-OR)` = fmt_num(estimate, 3),
      `95% CI (log-OR)` = ifelse(is.na(ci_low) | is.na(ci_high), "",
                                 paste0(fmt_num(ci_low,3), " to ", fmt_num(ci_high,3))),
      `OR [95% CI]` = ifelse(is.na(OR) | is.na(ORlo) | is.na(ORhi), "",
                             paste0(fmt_num(OR,2), " [", fmt_num(ORlo,2), "–", fmt_num(ORhi,2), "]")),
      `p` = fmt_p(pval),
      Stars = sig_stars_vec(pval),
      `I² (%)` = fmt_num(I2, 1),
      `τ²` = fmt_num(tau2, 3)
    ) %>%
    select(
      Molecule = molecule,
      `Adverse event` = ae_term,
      Model = model,
      Term = term,
      `β (log-OR)`, `95% CI (log-OR)`,
      `OR [95% CI]`,
      `z`, `p`, Stars,
      `k`, `I² (%)`, `τ²`, QE, QM, QMp
    )

  write_csv(x, outfile_csv)
  invisible(x)
}

make_all_paper_tables <- function(
    path_mol_agg,
    path_ae_mol_agg,
    path_models_ae,
    output_dir
) {
  message("Table 1: global dose–response by molecule …")
  t1 <- make_table_dr_global_by_molecule(path_mol_agg, output_dir)

  message("Table 2: AE × molecule matrix …")
  t2 <- make_table_matrix_ae_by_molecule(path_ae_mol_agg, output_dir, min_k = 2)

  message("Table 3: top AE per molecule …")
  t3 <- make_table_top_ae_per_molecule(path_ae_mol_agg, output_dir, top_n = 10, min_k = 2)

  message("Appendix: model terms …")
  app <- make_appendix_model_terms(path_models_ae, file.path(output_dir, "04_model_terms_appendix.csv"))

  invisible(list(t1 = t1, t2 = t2, t3 = t3, appendix = app))
}
