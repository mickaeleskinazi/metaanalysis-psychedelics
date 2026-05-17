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

fmt_star_p <- function(star, p){
  p_fmt <- fmt_p(p)
  ifelse(is.na(p) | p_fmt == "", "", ifelse(star == "", paste0("p=", p_fmt), paste(star, paste0("p=", p_fmt))))
}

paper_table_palette <- function(){
  list(
    header_bg = "#F4F1FA",
    stripe_bg = "#FAFAFC",
    border = "#D8D3E6",
    text = "#24212B",
    muted = "#5E5A66"
  )
}

style_gt_paper <- function(gt_tbl, source_note = NULL){
  pal <- paper_table_palette()
  gt_tbl <- gt_tbl %>%
    opt_table_font(font = list(google_font("Source Sans Pro"), default_fonts())) %>%
    opt_row_striping(row_striping = TRUE) %>%
    tab_options(
      table.background.color = "white",
      table.border.top.color = pal$border,
      table.border.top.width = px(1),
      table.border.bottom.color = pal$border,
      table.border.bottom.width = px(1),
      heading.border.bottom.color = pal$border,
      heading.title.font.size = px(16),
      heading.subtitle.font.size = px(12),
      column_labels.background.color = pal$header_bg,
      column_labels.font.weight = "bold",
      column_labels.border.top.color = pal$border,
      column_labels.border.bottom.color = pal$border,
      row.striping.background_color = pal$stripe_bg,
      table.font.color = pal$text,
      data_row.padding = px(5),
      source_notes.font.size = px(10)
    )
  if (!is.null(source_note)) {
    gt_tbl <- gt_tbl %>% tab_source_note(md(source_note))
  }
  gt_tbl
}

style_flextable_paper <- function(ft, sig_cols = character(), merge_first_col = FALSE){
  pal <- paper_table_palette()
  ft <- ft %>%
    theme_booktabs() %>%
    bg(part = "header", bg = pal$header_bg) %>%
    color(part = "all", color = pal$text) %>%
    bold(part = "header") %>%
    fontsize(part = "all", size = 8.5) %>%
    fontsize(part = "header", size = 9) %>%
    align(part = "all", align = "center") %>%
    valign(part = "all", valign = "center") %>%
    padding(part = "all", padding = 3)

  first_col <- names(ft$body$dataset)[[1]]
  ft <- ft %>% align(j = first_col, align = "left", part = "all")

  if (isTRUE(merge_first_col)) {
    ft <- ft %>%
      merge_v(j = first_col) %>%
      bold(j = first_col, part = "body")
  }

  for (col in intersect(sig_cols, names(ft$body$dataset))) {
    sig_rows <- grepl("\\*", ft$body$dataset[[col]])
    if (any(sig_rows, na.rm = TRUE)) {
      ft <- ft %>% bold(i = sig_rows, j = col, bold = TRUE, part = "body")
    }
  }

  ft %>% autofit()
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
    kable_styling(
      full_width = FALSE,
      position = "center",
      font_size = 8,
      latex_options = c("hold_position", "striped", "scale_down")
    )
  writeLines(tex, outfile_tex)
}

save_as_docx <- function(df, caption, outfile_docx, sig_cols = character(), merge_first_col = FALSE){
  dir.create(dirname(outfile_docx), recursive = TRUE, showWarnings = FALSE)
  ft <- flextable(df) %>%
    set_caption(caption) %>%
    style_flextable_paper(sig_cols = sig_cols, merge_first_col = merge_first_col)
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
      title = md("**Global dose-response by molecule**"),
      subtitle = md("_All adverse events pooled within each molecule_")
    ) %>%
    cols_align("center", columns = everything()) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(columns = c(`p (overall)`, Stars), rows = Stars != "")
    ) %>%
    style_gt_paper(source_note = "Stars denote model-level evidence for dose-response: * p < 0.05; ** p < 0.01; *** p < 0.001.")

  save_as_html(gt_tbl, file.path(output_dir, "01_dr_global_by_molecule.html"))
  save_as_tex(x,
              "Global dose-response by molecule",
              file.path(output_dir, "01_dr_global_by_molecule.tex"))
  save_as_docx(x,
               "Global dose-response by molecule",
               file.path(output_dir, "01_dr_global_by_molecule.docx"),
               sig_cols = "Stars")
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
              value = fmt_star_p(star, p)) %>%
    distinct() %>%
    pivot_wider(names_from = molecule, values_from = value) %>%
    mutate(across(-ae_term, ~ replace_na(.x, ""))) %>%
    arrange(ae_term)

  gt_tbl <- mat %>%
    gt() %>%
    tab_header(
      title = md("**Adverse events by molecule**"),
      subtitle = md("_Dose-response evidence for each adverse event-molecule pair_")
    ) %>%
    cols_align("center", columns = -ae_term) %>%
    cols_align("left", columns = ae_term) %>%
    cols_label(ae_term = "Adverse event") %>%
    style_gt_paper(source_note = "Blank cells indicate that no eligible model was available. Stars denote p-value thresholds.")

  for (col in setdiff(names(mat), "ae_term")) {
    gt_tbl <- gt_tbl %>%
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_body(columns = all_of(col), rows = grepl("\\*", .data[[col]]))
      )
  }

  save_as_html(gt_tbl, file.path(output_dir, "02_dr_by_ae_molecule_matrix.html"))
  save_as_tex(mat %>% rename(`Adverse event` = ae_term),
              "Adverse events by molecule",
              file.path(output_dir, "02_dr_by_ae_molecule_matrix.tex"))
  save_as_docx(mat %>% rename(`Adverse event` = ae_term),
               "Adverse events by molecule",
               file.path(output_dir, "02_dr_by_ae_molecule_matrix.docx"),
               sig_cols = setdiff(names(mat), "ae_term"))
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
      Sig. = star
    )

  gt_tbl <- top %>%
    gt() %>%
    tab_header(
      title = md("**Top adverse events by molecule**"),
      subtitle = md("_Ordered by overall dose-response p-value within molecule_")
    ) %>%
    tab_row_group(label = "AYAHUASCA", rows = Molecule == "AYAHUASCA") %>%
    tab_row_group(label = "LSD", rows = Molecule == "LSD") %>%
    tab_row_group(label = "MDMA", rows = Molecule == "MDMA") %>%
    tab_row_group(label = "PSILOCYBIN", rows = Molecule == "PSILOCYBIN") %>%
    cols_hide(columns = Molecule) %>%
    cols_align("center", columns = everything()) %>%
    cols_align("left", columns = `Adverse event`) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(columns = c(`p (overall)`, Sig.), rows = Sig. != "")
    ) %>%
    style_gt_paper(source_note = "Stars denote p-value thresholds for the selected model.")

  save_as_html(gt_tbl, file.path(output_dir, "03_topAE_per_molecule.html"))
  save_as_tex(top,
              "Top adverse events by molecule",
              file.path(output_dir, "03_topAE_per_molecule.tex"))
  save_as_docx(top,
               "Top adverse events by molecule",
               file.path(output_dir, "03_topAE_per_molecule.docx"),
               sig_cols = "Sig.",
               merge_first_col = TRUE)
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
