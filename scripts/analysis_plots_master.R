# ---- SIGNIFICANT-ONLY DR FACETS (AEs) ---------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(stringr); library(purrr); library(tidyr) })

# Small helper for flexible column names in preds
.normalize_pred_cols <- function(df){
  nm <- names(df)
  pick <- function(cands) { cands[cands %in% nm][1] }
  .x  <- pick(c("dose_norm","xdose","dose_mg","dose"))
  .y  <- pick(c("fit","pred","estimate","mu","y"))
  .l  <- pick(c("lwr","ci_low","ci_lb","ylwr"))
  .u  <- pick(c("upr","ci_high","ci_ub","yupr"))
  stopifnot(all(c("molecule","ae_term") %in% nm))
  if (is.na(.x) || is.na(.y) || is.na(.l) || is.na(.u)) {
    stop("Predictions must have dose & (fit,lwr,upr)-like columns. Found: ",
         paste(nm, collapse=", "))
  }
  df %>% rename(xdose = all_of(.x), fit = all_of(.y), lwr = all_of(.l), upr = all_of(.u))
}

# Compute omnibus p for AE×molecule from the model table
.omnibus_p_ae_mol <- function(models_by_ae){
  if (!nrow(models_by_ae)) return(tibble(molecule=character(), ae_term=character(), p_overall=double()))
  models_by_ae %>%
    group_by(molecule, ae_term) %>%
    summarise(
      p_overall = {
        qmp <- suppressWarnings(min(QMp[is.finite(QMp)], na.rm = TRUE))
        if (is.finite(qmp)) qmp else {
          # fallback: min p among non-intercept terms
          suppressWarnings(min(pval[is.finite(pval) & !grepl("intrcpt", term, TRUE)], na.rm = TRUE))
        }
      },
      .groups = "drop"
    )
}

# Main plotting function
plot_master_dr_significant_only <- function(
    preds, models, outfile = "results/master/master_dr_significant_only.pdf",
    normalize_dose = TRUE, ncol = 3, dpi = 300, width = 13, height = 9
){
  if (!nrow(preds) || !nrow(models)) {
    warning("No preds/models to plot."); return(invisible(NULL))
  }
  
  # Which AE×molecule are significant?
  sig_tbl <- .omnibus_p_ae_mol(models) %>%
    mutate(sig = is.finite(p_overall) & p_overall < 0.05)
  
  # Keep only AEs with ≥1 significant molecule
  ae_keep <- sig_tbl %>%
    group_by(ae_term) %>%
    summarise(any_sig = any(sig, na.rm = TRUE), .groups = "drop") %>%
    filter(any_sig) %>% pull(ae_term)
  
  if (!length(ae_keep)) {
    warning("No significant AEs found (p<0.05). Nothing to plot.")
    return(invisible(NULL))
  }
  
  # Choose one model per AE×molecule (prefer spline_df3 if present)
  choose_model <- preds %>%
    group_by(molecule, ae_term) %>%
    summarise(use_spline = any(model == "spline_df3"), .groups = "drop")
  
  df_all <- preds %>%
    inner_join(choose_model, by = c("molecule","ae_term")) %>%
    filter((use_spline & model == "spline_df3") | (!use_spline & model == "linear")) %>%
    filter(ae_term %in% ae_keep)
  
  if (!nrow(df_all)) {
    warning("No predictions after model selection/filter."); return(invisible(NULL))
  }
  
  # Normalize columns
  df_all <- .normalize_pred_cols(df_all)
  
  # Optionally normalize dose to 0–1 per molecule (within each AE)
  if (normalize_dose) {
    df_all <- df_all %>%
      group_by(ae_term, molecule) %>%
      mutate(
        .min = suppressWarnings(min(xdose, na.rm = TRUE)),
        .max = suppressWarnings(max(xdose, na.rm = TRUE)),
        xdose = ifelse(is.finite(.min) & is.finite(.max) & (.max > .min),
                       (xdose - .min) / (.max - .min),
                       0)  # collapse if constant
      ) %>% ungroup() %>% select(-.min, -.max)
  }
  
  # Attach significance per AE×molecule
  dfp <- df_all %>%
    left_join(sig_tbl, by = c("molecule","ae_term")) %>%
    mutate(
      sig    = ifelse(is.na(sig), FALSE, sig),
      star   = case_when(!is.finite(p_overall) ~ "",
                         p_overall < 0.001 ~ "***",
                         p_overall < 0.01  ~ "**",
                         p_overall < 0.05  ~ "*",
                         TRUE ~ ""),
      lty    = ifelse(sig, "significant", "ns"),
      alphaL = ifelse(sig, 1.0, 0.55)
    )
  
  # build plot
  p <- ggplot(dfp, aes(x = xdose, y = fit, color = molecule, group = molecule)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = molecule), alpha = 0.12, color = NA, show.legend = FALSE) +
    geom_line(aes(linetype = lty, alpha = alphaL), linewidth = 0.9) +
    scale_linetype_manual(values = c(significant = "solid", ns = "dashed"), name = "Model") +
    scale_alpha_continuous(range = c(0.55, 1.0), guide = "none") +
    scale_color_discrete(name = "Molecule") + scale_fill_discrete(guide = "none") +
    facet_wrap(~ ae_term, scales = "free_y", ncol = ncol) +
    labs(
      title = "Dose–response (significant AEs only)",
      subtitle = if (normalize_dose) "Dose normalized per molecule (0–1); solid lines: AE×molecule p<0.05" else "Solid lines: AE×molecule p<0.05",
      x = if (normalize_dose) "Normalized dose (0–1)" else "Dose (mg)",
      y = "log(OR) vs reference"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "white")
    )
  
  # Save
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  ggsave(outfile, p, width = width, height = height, dpi = dpi)
  message("Saved: ", outfile)
  invisible(outfile)
}