suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(readr)
})

# =============================================================================
# Helpers
# =============================================================================

.ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

.safe_filename <- function(x) {
  x %>%
    tolower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "_") %>%
    stringr::str_replace_all("^_+|_+$", "")
}

.require_cols <- function(df, cols, where = "data") {
  miss <- setdiff(cols, names(df))
  if (length(miss)) stop("Missing columns in ", where, ": ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

# =============================================================================
# 1) Plot: DR curves split by molecule (one PDF per molecule)
# Signature kept compatible with your pipeline:
#   plot_dr_by_molecule_split(preds, models, outdir)
# models is ignored (preds-only) but accepted for backward compatibility.
# =============================================================================
plot_dr_by_molecule_split <- function(preds, models = NULL, outdir) {
  .require_cols(preds, c("molecule", "dose_mg", "fit", "lwr", "upr"), where = "preds")
  .ensure_dir(outdir)
  
  mols <- sort(unique(preds$molecule))
  if (!length(mols)) return(invisible(NULL))
  
  for (mol in mols) {
    df <- preds %>% filter(molecule == mol) %>% arrange(dose_mg)
    if (!nrow(df)) next
    
    p <- ggplot(df, aes(x = dose_mg, y = fit)) +
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.20) +
      geom_line(linewidth = 0.9) +
      labs(
        title = paste0("Dose–response: ", mol),
        x = "Dose (mg)",
        y = "Effect size (log OR)"
      ) +
      theme_minimal(base_size = 12)
    
    ggsave(
      filename = file.path(outdir, paste0("dr_", .safe_filename(mol), ".pdf")),
      plot = p, width = 6.5, height = 4.5
    )
  }
  
  invisible(NULL)
}

# =============================================================================
# 2) Plot: per-molecule facets across AE (one PDF per molecule)
# Signature kept compatible: plot_dr_per_molecule_across_ae_facets(preds, models, outdir)
# =============================================================================
plot_dr_per_molecule_across_ae_facets <- function(preds, models = NULL, outdir,
                                                  max_ae_per_molecule = 16,
                                                  significant_only = FALSE,
                                                  p_threshold = 0.05) {
  .require_cols(preds, c("molecule", "ae_term", "dose_mg", "fit", "lwr", "upr"), where = "preds")
  .ensure_dir(outdir)
  
  # If asked, we can filter AE by significance *only if* models provides p-values.
  preds2 <- preds
  if (isTRUE(significant_only) && !is.null(models) && all(c("molecule","ae_term","pval") %in% names(models))) {
    keep <- models %>% filter(!is.na(pval), pval < p_threshold) %>% distinct(molecule, ae_term)
    preds2 <- preds2 %>% inner_join(keep, by = c("molecule","ae_term"))
  }
  
  mols <- sort(unique(preds2$molecule))
  if (!length(mols)) return(invisible(NULL))
  
  for (mol in mols) {
    df <- preds2 %>% filter(molecule == mol) %>% arrange(ae_term, dose_mg)
    if (!nrow(df)) next
    
    # limit number of AE facets (most frequent in preds)
    ae_rank <- df %>% count(ae_term, name = "npts") %>% arrange(desc(npts))
    ae_keep <- head(ae_rank$ae_term, max_ae_per_molecule)
    df <- df %>% filter(ae_term %in% ae_keep)
    
    p <- ggplot(df, aes(x = dose_mg, y = fit, group = ae_term)) +
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~ ae_term, scales = "free_y") +
      labs(
        title = paste0("Dose–response by AE (", mol, ")"),
        x = "Dose (mg)",
        y = "Effect size (log OR)"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        strip.text = element_text(size = 9),
        panel.spacing = unit(0.8, "lines")
      )
    
    ggsave(
      filename = file.path(outdir, paste0("dr_by_ae_facets_", .safe_filename(mol), ".pdf")),
      plot = p, width = 9.5, height = 7.2
    )
  }
  
  invisible(NULL)
}

# =============================================================================
# 3) Plot: per-AE curves on normalized dose axis (0–1 within molecule)
# Used in your pipeline: plot_dr_per_ae_normalized_dose(preds, outdir)
# =============================================================================
plot_dr_per_ae_normalized_dose <- function(preds, outdir,
                                           max_ae_total = 60,
                                           min_points_per_curve = 2) {
  .require_cols(preds, c("molecule", "ae_term", "dose_mg", "fit", "lwr", "upr"), where = "preds")
  .ensure_dir(outdir)
  
  # normalize dose within molecule
  df <- preds %>%
    group_by(molecule) %>%
    mutate(dose_norm = ifelse(max(dose_mg, na.rm = TRUE) > 0, dose_mg / max(dose_mg, na.rm = TRUE), NA_real_)) %>%
    ungroup() %>%
    filter(is.finite(dose_norm))
  
  if (!nrow(df)) return(invisible(NULL))
  
  # keep AE with enough points
  df <- df %>%
    group_by(molecule, ae_term) %>%
    mutate(npts = n()) %>%
    ungroup() %>%
    filter(npts >= min_points_per_curve)
  
  if (!nrow(df)) return(invisible(NULL))
  
  # choose top AE by number of points across all molecules
  top_ae <- df %>% count(ae_term, sort = TRUE) %>% slice_head(n = max_ae_total) %>% pull(ae_term)
  df <- df %>% filter(ae_term %in% top_ae)
  
  p <- ggplot(df, aes(x = dose_norm, y = fit, group = interaction(molecule, ae_term))) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.12) +
    geom_line(linewidth = 0.6) +
    facet_grid(molecule ~ ae_term, scales = "free_y") +
    labs(
      title = "Dose–response by AE on normalized dose axis (within molecule)",
      x = "Normalized dose (0–1 within molecule)",
      y = "Effect size (log OR)"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      strip.text.x = element_text(size = 7),
      strip.text.y = element_text(size = 8),
      panel.spacing = unit(0.4, "lines")
    )
  
  ggsave(
    filename = file.path(outdir, "dr_by_ae_normalized_dose.pdf"),
    plot = p, width = 16, height = 9
  )
  
  invisible(NULL)
}
