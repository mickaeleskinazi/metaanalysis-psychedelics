# ============================================================
# analysis_plots_master.R
# Master plots: dose–response global, forest globaux, DR par AE
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(metafor)
})



# ---------- A) MASTER DR: global dose–response by molecule ----------
# preds_mol : dr_mol$preds
# colonnes attendues: molecule, dose_mg, fit, lwr, upr
plot_master_dr_by_molecule <- function(preds_mol, outfile = "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/results/master/master_dr_by_molecule.pdf"){
  stopifnot(all(c("molecule","dose_mg","fit","lwr","upr") %in% names(preds_mol)))
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  
  p <- ggplot2::ggplot(preds_mol, ggplot2::aes(x = dose_mg, y = fit)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr),
                         alpha = 0.18, fill = "grey80", color = NA) +
    ggplot2::geom_line(size = 1, color = "black") +
    ggplot2::facet_wrap(~ molecule, scales = "free", nrow = 1) +   # côte à côte
    ggplot2::labs(x = "Dose (mg)", y = "log(OR) vs ref",
                  title = "Dose-response globale par molécule (tous AE confondus)") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"))
  
  ggplot2::ggsave(outfile, p, width = 14, height = 4.5, dpi = 300)
  invisible(p)
}



# ---------- B) FOREST GLOBAUX par molécule ----------
# Helper: pooled OR per (molecule × AE)
sig_stars_vec <- function(p){
  out <- rep("", length(p))
  out[!is.na(p) & p < 0.001] <- "***"
  out[!is.na(p) & p >= 0.001 & p < 0.01] <- "**"
  out[!is.na(p) & p >= 0.01  & p < 0.05] <- "*"
  out
}

pooled_by_molecule_ae <- function(es){
  stopifnot(all(c("yi","vi","molecule","ae_term") %in% names(es)))
  es %>%
    dplyr::filter(is.finite(yi), is.finite(vi)) %>%
    dplyr::group_by(molecule, ae_term) %>%
    dplyr::group_modify(function(.x, .key){
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
    dplyr::ungroup() %>%
    dplyr::mutate(
      stars = sig_stars_vec(pval),
      sig   = !is.na(pval) & pval < 0.05
    )
}

plot_master_forest_by_molecule <- function(
    es,
    outfile = "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/results/master/master_forest_by_molecule.pdf",
    min_k = 2
){
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  
  pooled <- pooled_by_molecule_ae(es) %>% dplyr::filter(k >= min_k)
  if (!nrow(pooled)) { warning("Aucun AE avec k >= ", min_k); return(invisible(NULL)) }
  
  # Ordre : significatifs d'abord, puis p croissante
  pooled <- pooled %>%
    dplyr::group_by(molecule) %>%
    dplyr::arrange(dplyr::desc(sig), pval, .by_group = TRUE) %>%
    dplyr::mutate(ae_term = factor(ae_term, levels = unique(ae_term))) %>%
    dplyr::ungroup()
  
  p <- ggplot2::ggplot(pooled, ggplot2::aes(x = or, y = ae_term)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lo, xmax = ci_hi, color = sig),
                            height = 0, na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(color = sig), size = 2, na.rm = TRUE) +
    ggplot2::geom_text(ggplot2::aes(label = stars, color = sig),
                       nudge_x = 0.02, hjust = 0, size = 4, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "grey30"), guide = "none") +
    ggplot2::scale_x_log10() +
    ggplot2::facet_wrap(~ molecule, scales = "free_y", nrow = 1) +   # côte à côte
    ggplot2::labs(x = "OR (log)", y = NULL, title = "Forest plot global (pooled par AE)") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"))
  
  ggplot2::ggsave(outfile, p, width = 14, height = 5.5, dpi = 300)
  invisible(p)
}

