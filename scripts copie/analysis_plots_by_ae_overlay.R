# scripts/analysis_plots_by_ae_overlay.R
# ------------------------------------------------------------
# Overlay plots "by AE":
#  - Dose–response curves per AE, molecules overlaid
#  - Forest plots per AE, one line per molecule
# Significance:
#  - From model terms: any non-intercept "dose" term with p < 0.05
#  - DR overlay: solid line if significant; stars at end of curve
#  - Forest overlay: points/CI in red if significant
# Dependencies: dplyr, ggplot2, purrr, stringr, metafor
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(stringr)
  library(metafor)
})

# ------------------------------------------------------------
# Helper: significance stars
# ------------------------------------------------------------
sig_stars_vec <- function(p){
  out <- rep("", length(p))
  out[!is.na(p) & p < 0.001] <- "***"
  out[!is.na(p) & p >= 0.001 & p < 0.01] <- "**"
  out[!is.na(p) & p >= 0.01  & p < 0.05] <- "*"
  out
}

# ------------------------------------------------------------
# DR overlay per AE:
#   preds: dr_ae$preds with cols {molecule, ae_term, dose_mg, fit, lwr, upr}
#   models: dr_ae$models with cols {molecule, ae_term, term, pval, ...}
# Behavior:
#   - keeps only AE present in >= min_molecules molecules
#   - optional dose normalization to [0,1] within AE×molecule for comparability
#   - per AE -> one PDF (molecules colored, solid if significant, stars at end)
# ------------------------------------------------------------
plot_dr_overlay_per_ae <- function(
    preds,
    models,
    outdir = "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/results/master/dr_overlay_by_ae",
    min_molecules = 2,
    normalize_dose = TRUE,
    width = 8, height = 5, dpi = 300
){
  stopifnot(all(c("molecule","ae_term","dose_mg","fit","lwr","upr") %in% names(preds)))
  stopifnot(all(c("molecule","ae_term","term","pval") %in% names(models)))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # 1) Significance per AE×molecule from model terms (exclude intercept-like terms)
  sig_tbl <- models %>%
    mutate(is_dose = !grepl("intrcpt", term, ignore.case = TRUE)) %>%
    group_by(ae_term, molecule) %>%
    summarise(
      min_p = suppressWarnings(min(pval[is_dose], na.rm = TRUE)),
      sig   = is.finite(min_p) & (min_p < 0.05),
      stars = sig_stars_vec(min_p),
      .groups = "drop"
    )
  
  # 2) Keep only AE that occur in >= min_molecules molecules
  ae_common <- preds %>%
    distinct(ae_term, molecule) %>%
    count(ae_term, name = "n_mol") %>%
    filter(n_mol >= min_molecules)
  
  # 3) Dose normalization to [0,1] within AE×molecule (optional)
  if (normalize_dose) {
    dose_ref <- preds %>%
      inner_join(ae_common, by = "ae_term") %>%
      group_by(ae_term, molecule) %>%
      summarise(
        dose_max = max(dose_mg, na.rm = TRUE),
        dose_min = min(dose_mg, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(span = ifelse((dose_max - dose_min) > 0, dose_max - dose_min, NA_real_))
    
    preds_use <- preds %>%
      inner_join(ae_common, by = "ae_term") %>%
      left_join(dose_ref, by = c("ae_term","molecule")) %>%
      mutate(
        xdose = ifelse(!is.na(span), (dose_mg - dose_min) / span, dose_mg),
        xlab  = "Dose normalised (0-1, per molecule)"
      ) %>%
      select(-dose_max, -dose_min, -span)
  } else {
    preds_use <- preds %>%
      inner_join(ae_common, by = "ae_term") %>%
      mutate(xdose = dose_mg, xlab = "Dose (mg)")
  }
  
  # 4) Join stars/significance
  preds_use <- preds_use %>% left_join(sig_tbl, by = c("ae_term","molecule"))
  
  # 5) One figure per AE
  ae_list <- sort(unique(preds_use$ae_term))
  walk(ae_list, function(ae){
    dat <- preds_use %>% filter(ae_term == ae) %>% arrange(molecule, xdose)
    if (n_distinct(dat$molecule) < min_molecules) return(NULL)
    
    # order molecules: significant first, then by p, then name
    mol_order <- dat %>%
      distinct(molecule, sig, min_p) %>%
      arrange(desc(sig), min_p, molecule) %>%
      pull(molecule)
    dat$molecule <- factor(dat$molecule, levels = mol_order)
    
    # star positions = last x for each molecule (end of curve), nudged right
    xr <- range(dat$xdose, na.rm = TRUE)
    dx <- if (is.finite(diff(xr)) && diff(xr) > 0) 0.02 * diff(xr) else 0.02
    
    stars_df <- dat %>%
      group_by(molecule) %>%
      slice_max(order_by = xdose, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      filter(sig, !is.na(fit), !is.na(xdose)) %>%
      mutate(
        label  = stars,
        x_text = pmin(xdose + dx, xr[2])  # keep inside panel
      )
    
    p <- ggplot(dat, aes(x = xdose, y = fit, color = molecule, linetype = sig, group = molecule)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_ribbon(aes(ymin = lwr, ymax = upr, fill = molecule),
                  alpha = 0.12, color = NA, show.legend = FALSE) +
      geom_line(size = 1.05) +
      # stars on curve ends
      geom_text(
        data = stars_df,
        aes(x = x_text, y = fit, label = label),
        inherit.aes = FALSE,
        size = 5,
        show.legend = FALSE
      ) +
      scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dashed"), guide = "none") +
      scale_color_discrete() + scale_fill_discrete(guide = "none") +
      labs(
        title    = paste0("Dose-response for AE: ", ae),
        x        = unique(dat$xlab),
        y        = "log(OR) vs reference",
        subtitle = "Stars mark AE×molecule with p<0.05 for >=1 non-intercept dose term (solid if significant)"
      ) +
      theme_bw() +
      theme(
        legend.position  = "bottom",
        strip.background = element_rect(fill = "white")
      )
    
    ggsave(
      file.path(outdir, paste0("DR_overlay_", stringr::str_replace_all(ae, "[^A-Za-z0-9]+","_"), ".pdf")),
      p, width = width, height = height, dpi = dpi
    )
  })
  
  invisible(TRUE)
}

# ------------------------------------------------------------
# Forest overlay per AE:
#   es: escalc data with cols {yi, vi, molecule, ae_term}
#   Pooled random-effects per AE×molecule; red if p<0.05
#   One PDF per AE
# ------------------------------------------------------------
plot_forest_overlay_per_ae <- function(
    es,
    outdir = "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/results/master/forest_overlay_by_ae",
    min_k = 2,
    min_molecules = 2,
    width = 7, height = 4, dpi = 300
){
  stopifnot(all(c("yi","vi","molecule","ae_term") %in% names(es)))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  pooled <- es %>%
    filter(is.finite(yi), is.finite(vi)) %>%
    group_by(ae_term, molecule) %>%
    group_modify(function(.x, .key){
      if (nrow(.x) < min_k) {
        return(tibble(or = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, k = nrow(.x)))
      }
      fit <- tryCatch(rma(yi, vi, data = .x, method = "REML"), error = function(e) NULL)
      if (is.null(fit)) {
        return(tibble(or = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, k = nrow(.x)))
      }
      tibble(
        or = exp(as.numeric(fit$b[1])),
        lo = exp(as.numeric(fit$ci.lb)),
        hi = exp(as.numeric(fit$ci.ub)),
        p  = suppressWarnings(as.numeric(fit$pval))[1],
        k  = fit$k
      )
    }) %>%
    ungroup() %>%
    mutate(sig = !is.na(p) & p < 0.05)
  
  # keep only AE present in >= min_molecules molecules
  ae_keep <- pooled %>%
    filter(!is.na(or)) %>%
    distinct(ae_term, molecule) %>%
    count(ae_term, name = "n_mol") %>%
    filter(n_mol >= min_molecules)
  
  pooled <- pooled %>% semi_join(ae_keep, by = "ae_term")
  
  ae_list <- pooled %>% pull(ae_term) %>% unique() %>% sort()
  walk(ae_list, function(ae){
    dat <- pooled %>% filter(ae_term == ae, !is.na(or))
    if (!nrow(dat)) return(NULL)
    
    # order: significant first, then by p, then name
    dat <- dat %>% arrange(desc(sig), p, molecule)
    dat$molecule <- factor(dat$molecule, levels = unique(dat$molecule))
    
    p <- ggplot(dat, aes(x = or, y = molecule, color = sig)) +
      geom_vline(xintercept = 1, linetype = "dashed") +
      geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, size = 0.9) +
      geom_point(size = 2) +
      scale_x_log10() +
      scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "grey30"), guide = "none") +
      labs(
        title    = paste0("Forest per AE: ", ae),
        x        = "OR (log-scale)",
        y        = NULL,
        subtitle = "Red = p<0.05 (random-effects pooled per molecule)"
      ) +
      theme_bw() +
      theme(strip.background = element_rect(fill = "white"))
    
    ggsave(
      file.path(outdir, paste0("FOREST_overlay_", stringr::str_replace_all(ae, "[^A-Za-z0-9]+","_"), ".pdf")),
      p, width = width, height = height, dpi = dpi
    )
  })
  
  invisible(TRUE)
}