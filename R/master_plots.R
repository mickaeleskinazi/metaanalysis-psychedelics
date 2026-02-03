suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

.require_cols <- function(df, cols, where = "data") {
  miss <- setdiff(cols, names(df))
  if (length(miss)) stop("Missing columns in ", where, ": ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

.safe_filename <- function(x) {
  x %>%
    tolower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "_") %>%
    stringr::str_replace_all("^_+|_+$", "")
}

# =============================================================================
# Master DR by molecule (single figure)
# =============================================================================
plot_master_dr_by_molecule <- function(preds, outfile) {
  .require_cols(preds, c("molecule", "dose_mg", "fit", "lwr", "upr"), where = "preds")
  
  df <- preds %>% arrange(molecule, dose_mg)
  if (!nrow(df)) return(invisible(NULL))
  
  p <- ggplot(df, aes(x = dose_mg, y = fit, group = molecule)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.18) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ molecule, scales = "free_x") +
    labs(
      title = "Global dose–response by molecule",
      x = "Dose (mg)",
      y = "Effect size (log OR)"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.spacing = unit(1.0, "lines"))
  
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  ggsave(outfile, p, width = 10, height = 4.8)
  
  invisible(NULL)
}

# =============================================================================
# Master DR by AE (single figure)
# Your pipeline calls:
# plot_master_dr_by_ae(preds, models, outfile, max_ae_per_molecule, significant_only)
# We implement "significant_only" only if models has pval; otherwise fallback to top-k.
# =============================================================================
plot_master_dr_by_ae <- function(preds,
                                 models = NULL,
                                 outfile,
                                 max_ae_per_molecule = 20,
                                 significant_only = TRUE,
                                 p_threshold = 0.05) {
  
  .require_cols(preds, c("molecule", "ae_term", "dose_mg", "fit", "lwr", "upr"), where = "preds")
  
  df <- preds %>% arrange(molecule, ae_term, dose_mg)
  if (!nrow(df)) return(invisible(NULL))
  
  # Select AE terms per molecule
  df_sel <- df
  
  if (isTRUE(significant_only) &&
      !is.null(models) &&
      all(c("molecule","ae_term","pval") %in% names(models))) {
    
    keep <- models %>%
      filter(!is.na(pval), pval < p_threshold) %>%
      distinct(molecule, ae_term)
    
    df_sel <- df_sel %>% inner_join(keep, by = c("molecule","ae_term"))
  }
  
  # If still too many AEs, keep top max_ae_per_molecule by number of points per molecule
  df_sel <- df_sel %>%
    group_by(molecule, ae_term) %>% mutate(npts = n()) %>% ungroup()
  
  ae_keep <- df_sel %>%
    count(molecule, ae_term, wt = npts, name = "score") %>%
    group_by(molecule) %>%
    arrange(desc(score)) %>%
    slice_head(n = max_ae_per_molecule) %>%
    ungroup() %>%
    select(molecule, ae_term)
  
  df_sel <- df_sel %>% inner_join(ae_keep, by = c("molecule","ae_term"))
  
  if (!nrow(df_sel)) {
    # fallback: plot nothing but keep pipeline alive
    warning("plot_master_dr_by_ae(): no AE curves selected (check models/p-values or max_ae_per_molecule).")
    return(invisible(NULL))
  }
  
  p <- ggplot(df_sel, aes(x = dose_mg, y = fit, group = interaction(molecule, ae_term))) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.12) +
    geom_line(linewidth = 0.6) +
    facet_grid(molecule ~ ae_term, scales = "free_y") +
    labs(
      title = "Dose–response by AE (selected)",
      x = "Dose (mg)",
      y = "Effect size (log OR)"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      strip.text.x = element_text(size = 7),
      strip.text.y = element_text(size = 8),
      panel.spacing = unit(0.4, "lines")
    )
  
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  ggsave(outfile, p, width = 16, height = 9)
  
  invisible(NULL)
}
