# ---- Master plots for dose–response & pooled effects -------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(metafor)
})

# Small helper for flexible column names in preds (AE × molecule)
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
         paste(nm, collapse = ", "))
  }
  df %>% rename(xdose = all_of(.x), fit = all_of(.y), lwr = all_of(.l), upr = all_of(.u))
}

# Variant for molecule-only predictions
.normalize_pred_cols_molecule <- function(df){
  nm <- names(df)
  pick <- function(cands) { cands[cands %in% nm][1] }
  .x  <- pick(c("dose_mg","xdose","dose","dose_norm"))
  .y  <- pick(c("fit","pred","estimate","mu","y"))
  .l  <- pick(c("lwr","ci_low","ci_lb","ylwr"))
  .u  <- pick(c("upr","ci_high","ci_ub","yupr"))
  stopifnot("molecule" %in% nm)
  if (is.na(.x) || is.na(.y) || is.na(.l) || is.na(.u)) {
    stop("Predictions must have dose & (fit,lwr,upr)-like columns. Found: ",
         paste(nm, collapse = ", "))
  }
  df %>% rename(dose = all_of(.x), fit = all_of(.y), lwr = all_of(.l), upr = all_of(.u))
}

# Compute omnibus p for AE×molecule from the model table
.omnibus_p_ae_mol <- function(models_by_ae){
  if (!nrow(models_by_ae)) {
    return(tibble(molecule = character(), ae_term = character(), p_overall = double()))
  }
  models_by_ae %>%
    group_by(molecule, ae_term) %>%
    summarise(
      p_overall = {
        qmp <- suppressWarnings(min(QMp[is.finite(QMp)], na.rm = TRUE))
        if (is.finite(qmp)) qmp else {
          suppressWarnings(min(pval[is.finite(pval) & !grepl("intrcpt", term, TRUE)], na.rm = TRUE))
        }
      },
      .groups = "drop",
    )
}

.sig_stars_vec <- function(p){
  out <- rep("", length(p))
  out[!is.na(p) & p < 0.001] <- "***"
  out[!is.na(p) & p >= 0.001 & p < 0.01] <- "**"
  out[!is.na(p) & p >= 0.01  & p < 0.05] <- "*"
  out
}

.pooled_by_molecule_ae <- function(es){
  stopifnot(all(c("yi","vi","molecule","ae_term") %in% names(es)))
  es %>%
    filter(is.finite(yi), is.finite(vi)) %>%
    group_by(molecule, ae_term) %>%
    group_modify(function(.x, .key){
      k <- nrow(.x)
      if (k < 2) {
        tibble::tibble(or = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, pval = NA_real_, k = k)
      } else {
        fit <- tryCatch(metafor::rma(yi, vi, data = .x, method = "REML"), error = function(e) NULL)
        if (is.null(fit)) {
          tibble::tibble(or = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, pval = NA_real_, k = k)
        } else {
          tibble::tibble(
            or    = exp(as.numeric(fit$b[1])),
            ci_lo = exp(as.numeric(fit$ci.lb)),
            ci_hi = exp(as.numeric(fit$ci.ub)),
            pval  = suppressWarnings(as.numeric(fit$pval))[1],
            k     = fit$k
          )
        }
      }
    }) %>%
    ungroup() %>%
    mutate(
      stars = .sig_stars_vec(pval),
      sig   = !is.na(pval) & pval < 0.05
    )
}

.plot_master_dr_by_ae_impl <- function(
    preds,
    models,
    outfile,
    normalize_dose = TRUE,
    ncol = 3,
    dpi = 300,
    width = 13,
    height = 9,
    significant_only = TRUE,
    max_ae_per_molecule = Inf
){
  if (!nrow(preds)) {
    warning("No predictions to plot.")
    return(invisible(NULL))
  }

  combos <- preds %>% distinct(molecule, ae_term)
  if (!nrow(combos)) {
    warning("Predictions do not include molecule/AE columns.")
    return(invisible(NULL))
  }

  if (nrow(models)) {
    stats <- .omnibus_p_ae_mol(models)
    sig_tbl <- combos %>% left_join(stats, by = c("molecule","ae_term"))
  } else {
    if (significant_only) {
      warning("No models available to determine significance. Nothing to plot.")
      return(invisible(NULL))
    }
    sig_tbl <- combos %>% mutate(p_overall = NA_real_)
  }

  sig_tbl <- sig_tbl %>% mutate(sig = is.finite(p_overall) & p_overall < 0.05)

  if (significant_only) {
    ae_keep <- sig_tbl %>%
      group_by(ae_term) %>%
      summarise(any_sig = any(sig, na.rm = TRUE), .groups = "drop") %>%
      filter(any_sig) %>% pull(ae_term)
    if (!length(ae_keep)) {
      warning("No significant AEs found (p<0.05). Nothing to plot.")
      return(invisible(NULL))
    }
  } else {
    ae_keep <- unique(preds$ae_term)
  }

  choose_model <- preds %>%
    group_by(molecule, ae_term) %>%
    summarise(use_spline = any(model == "spline_df3"), .groups = "drop")

  df_all <- preds %>%
    inner_join(choose_model, by = c("molecule","ae_term")) %>%
    filter((use_spline & model == "spline_df3") | (!use_spline & model == "linear")) %>%
    filter(ae_term %in% ae_keep)

  if (!nrow(df_all)) {
    warning("No predictions after model selection/filter.")
    return(invisible(NULL))
  }

  if (is.finite(max_ae_per_molecule)) {
    keep_pairs <- sig_tbl %>%
      mutate(p_rank = ifelse(is.finite(p_overall), p_overall, Inf),
             sig_rank = ifelse(sig, 0, 1)) %>%
      arrange(molecule, sig_rank, p_rank, ae_term) %>%
      group_by(molecule) %>%
      mutate(idx = row_number()) %>%
      filter(idx <= max_ae_per_molecule) %>%
      ungroup() %>%
      select(molecule, ae_term) %>%
      distinct()

    if (!nrow(keep_pairs)) {
      keep_pairs <- df_all %>% distinct(molecule, ae_term)
    }

    df_all <- df_all %>% semi_join(keep_pairs, by = c("molecule","ae_term"))
    sig_tbl <- sig_tbl %>% semi_join(keep_pairs, by = c("molecule","ae_term"))
  }

  df_all <- .normalize_pred_cols(df_all)

  if (normalize_dose) {
    df_all <- df_all %>%
      group_by(ae_term, molecule) %>%
      mutate(
        .min = suppressWarnings(min(xdose, na.rm = TRUE)),
        .max = suppressWarnings(max(xdose, na.rm = TRUE)),
        xdose = ifelse(is.finite(.min) & is.finite(.max) & (.max > .min),
                       (xdose - .min) / (.max - .min),
                       0)
      ) %>% ungroup() %>% select(-.min, -.max)
  }

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

  title_txt <- if (significant_only) {
    "Dose–response (significant AEs only)"
  } else {
    "Dose–response by adverse event"
  }

  subtitle_txt <- if (normalize_dose) {
    if (significant_only) {
      "Dose normalized per molecule (0–1); solid lines: AE×molecule p<0.05"
    } else {
      "Dose normalized per molecule (0–1)"
    }
  } else if (significant_only) {
    "Solid lines: AE×molecule p<0.05"
  } else {
    NULL
  }

  p <- ggplot(dfp, aes(x = xdose, y = fit, color = molecule, group = molecule)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = molecule), alpha = 0.12, color = NA, show.legend = FALSE) +
    geom_line(aes(linetype = lty, alpha = alphaL), linewidth = 0.9) +
    scale_linetype_manual(values = c(significant = "solid", ns = "dashed"), name = "Model") +
    scale_alpha_continuous(range = c(0.55, 1.0), guide = "none") +
    scale_color_discrete(name = "Molecule") + scale_fill_discrete(guide = "none") +
    facet_wrap(~ ae_term, scales = "free_y", ncol = ncol) +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = if (normalize_dose) "Normalized dose (0–1)" else "Dose (mg)",
      y = "log(OR) vs reference"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "white")
    )

  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  ggsave(outfile, p, width = width, height = height, dpi = dpi)
  message("Saved: ", outfile)
  invisible(outfile)
}

plot_master_dr_significant_only <- function(
    preds,
    models,
    outfile = "results/master/master_dr_significant_only.pdf",
    normalize_dose = TRUE,
    ncol = 3,
    dpi = 300,
    width = 13,
    height = 9,
    max_ae_per_molecule = Inf
){
  .plot_master_dr_by_ae_impl(
    preds = preds,
    models = models,
    outfile = outfile,
    normalize_dose = normalize_dose,
    ncol = ncol,
    dpi = dpi,
    width = width,
    height = height,
    significant_only = TRUE,
    max_ae_per_molecule = max_ae_per_molecule
  )
}

plot_master_dr_by_ae <- function(
    preds,
    models,
    outfile = "results/master/master_dr_by_ae.pdf",
    max_ae_per_molecule = Inf,
    significant_only = TRUE,
    normalize_dose = TRUE,
    ncol = 3,
    dpi = 300,
    width = 13,
    height = 9
){
  .plot_master_dr_by_ae_impl(
    preds = preds,
    models = models,
    outfile = outfile,
    normalize_dose = normalize_dose,
    ncol = ncol,
    dpi = dpi,
    width = width,
    height = height,
    significant_only = significant_only,
    max_ae_per_molecule = max_ae_per_molecule
  )
}

plot_master_dr_by_molecule <- function(
    preds_mol,
    outfile = "results/master/master_dr_by_molecule.pdf",
    width = 14,
    height = 4.5,
    dpi = 300
){
  if (is.null(preds_mol) || !nrow(preds_mol)) {
    warning("No predictions supplied for master dose–response by molecule plot.")
    return(invisible(NULL))
  }

  preds_use <- preds_mol
  if ("model" %in% names(preds_use)) {
    choose <- preds_use %>%
      group_by(molecule) %>%
      summarise(use_spline = any(model == "spline_df3"), .groups = "drop")
    preds_use <- preds_use %>%
      inner_join(choose, by = "molecule") %>%
      filter((use_spline & model == "spline_df3") | (!use_spline & model == "linear"))
  }

  preds_use <- .normalize_pred_cols_molecule(preds_use)

  p <- ggplot(preds_use, aes(x = dose, y = fit)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.18, fill = "grey80", color = NA) +
    geom_line(linewidth = 0.9, color = "black") +
    facet_wrap(~ molecule, scales = "free", nrow = 1) +
    labs(
      x = "Dose (mg)",
      y = "log(OR) vs reference",
      title = "Dose–response by molecule (all adverse events)"
    ) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"))

  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

  pooled <- .pooled_by_molecule_ae(es) %>% filter(k >= min_k)
  if (!nrow(pooled)) {
    warning("No adverse events with k >= ", min_k, ".")
    return(invisible(NULL))
  }

  pooled <- pooled %>%
    group_by(molecule) %>%
    arrange(desc(sig), pval, .by_group = TRUE) %>%
    mutate(ae_term = factor(ae_term, levels = unique(ae_term))) %>%
    ungroup()

  p <- ggplot(pooled, aes(x = or, y = ae_term)) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi, color = sig), height = 0, na.rm = TRUE) +
    geom_point(aes(color = sig), size = 2, na.rm = TRUE) +
    geom_text(aes(label = stars, color = sig), nudge_x = 0.02, hjust = 0, size = 4, na.rm = TRUE) +
    scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "grey30"), guide = "none") +
    scale_x_log10() +
    facet_wrap(~ molecule, scales = "free_y", nrow = 1) +
    labs(x = "Odds ratio", y = NULL, title = "Forest plot by molecule (pooled across AEs)") +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"))

  ggsave(outfile, p, width = width, height = height, dpi = dpi)
  message("Saved: ", outfile)
  invisible(outfile)
}

plot_master_forest_by_molecule <- function(
    es,
    outfile = "results/master/master_forest_by_molecule.pdf",
    min_k = 2,
    width = 14,
    height = 5.5,
    dpi = 300
){
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

  pooled <- .pooled_by_molecule_ae(es) %>% filter(k >= min_k)
  if (!nrow(pooled)) {
    warning("No adverse events with k >= ", min_k, ".")
    return(invisible(NULL))
  }

  pooled <- pooled %>%
    group_by(molecule) %>%
    arrange(desc(sig), pval, .by_group = TRUE) %>%
    mutate(ae_term = factor(ae_term, levels = unique(ae_term))) %>%
    ungroup()

  p <- ggplot(pooled, aes(x = or, y = ae_term)) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi, color = sig), height = 0, na.rm = TRUE) +
    geom_point(aes(color = sig), size = 2, na.rm = TRUE) +
    geom_text(aes(label = stars, color = sig), nudge_x = 0.02, hjust = 0, size = 4, na.rm = TRUE) +
    scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "grey30"), guide = "none") +
    scale_x_log10() +
    facet_wrap(~ molecule, scales = "free_y", nrow = 1) +
    labs(x = "Odds ratio", y = NULL, title = "Forest plot by molecule (pooled across AEs)") +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"))

  ggsave(outfile, p, width = width, height = height, dpi = dpi)
  message("Saved: ", outfile)
  invisible(outfile)
}

