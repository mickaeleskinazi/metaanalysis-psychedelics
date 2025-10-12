# scripts/forest_compare_session_followup.R
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(metafor)
  library(stringr)
  library(patchwork)  # pour combiner les figures et partager la légende
})

# ---- Utilitaires internes ----

# Résume (pooled) par AE x time_window pour un sous-ensemble (une molécule)
.summarise_ae_by_window <- function(dmol, min_k_per_window = 2){
  dmol %>%
    filter(is.finite(yi), is.finite(vi)) %>%
    group_by(ae_term, time_window) %>%
    group_map(~{
      g <- .x
      if (nrow(g) < min_k_per_window) return(NULL)
      fit <- try(rma(yi, vi, data = g, method = "REML"), silent = TRUE)
      if (inherits(fit, "try-error")) return(NULL)
      tibble(
        ae_term    = g$ae_term[[1]],
        time_window= g$time_window[[1]],
        k          = nrow(g),
        est        = as.numeric(fit$b[1,1]),
        ci_lb      = as.numeric(fit$ci.lb),
        ci_ub      = as.numeric(fit$ci.ub),
        pval       = as.numeric(fit$pval)
      )
    }, .keep = TRUE) %>% list_rbind()
}

# Ordonne les AE pour lisibilité
.order_ae <- function(sumdf, order_by = c("pval","|est|","name")){
  order_by <- match.arg(order_by)
  ord <- sumdf %>%
    group_by(ae_term) %>%
    summarise(
      p_best = suppressWarnings(min(pval, na.rm = TRUE)),
      est_max = max(abs(est), na.rm = TRUE),
      .groups = "drop"
    )
  if (order_by == "pval") {
    ord <- ord %>% arrange(is.na(p_best), p_best, desc(est_max))
  } else if (order_by == "|est|") {
    ord <- ord %>% arrange(desc(est_max), p_best)
  } else {
    ord <- ord %>% arrange(ae_term)
  }
  ord$ae_term
}

# ---- 1A. Forest comparatif (un PDF par molécule) ----
forest_compare_all_molecules <- function(es_all,
                                         outdir = "results_compare/forest_by_ae_side_by_side",
                                         min_k_per_window = 2,
                                         order_by = c("pval","|est|","name"),
                                         star_offset = 0.06) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  order_by <- match.arg(order_by)
  
  stopifnot(all(c("yi","vi","study_id","molecule","ae_term","time_window") %in% names(es_all)))
  
  es_all <- es_all %>%
    filter(is.finite(yi), is.finite(vi)) %>%
    mutate(time_window = factor(time_window, levels = c("session","follow_up")))
  
  mols <- sort(unique(es_all$molecule))
  
  for (mol in mols) {
    dmol <- es_all %>% filter(molecule == !!mol)
    if (!nrow(dmol)) next
    
    sumdf <- .summarise_ae_by_window(dmol, min_k_per_window = min_k_per_window)
    if (is.null(sumdf) || !nrow(sumdf)) next
    
    sumdf <- sumdf %>%
      mutate(
        stars = case_when(
          is.na(pval) ~ "",
          pval < 0.001 ~ "***",
          pval < 0.01  ~ "**",
          pval < 0.05  ~ "*",
          TRUE ~ ""
        ),
        # décalage horizontal dynamique des étoiles (différent par fenêtre)
        x_star = ci_ub + ifelse(time_window == "session",  star_offset, -star_offset)
      )
    
    # Ordre des AE
    levels_ae <- .order_ae(sumdf, order_by = order_by)
    sumdf <- sumdf %>% mutate(ae_term = factor(ae_term, levels = rev(levels_ae)))
    
    p <- ggplot(sumdf, aes(y = ae_term, x = est, xmin = ci_lb, xmax = ci_ub, color = time_window)) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.45) +
      geom_errorbarh(height = 0.22, position = position_dodge(width = 0.6)) +
      geom_point(size = 2.4, position = position_dodge(width = 0.6)) +
      # étoiles décalées (pas de position + nudge ensemble)
      geom_text(aes(x = x_star, label = stars), color = "red3", size = 4, show.legend = FALSE) +
      scale_color_manual(values = c(session = "#1f77b4", follow_up = "#ff7f0e"),
                         name = NULL, labels = c("Session","Follow-up")) +
      labs(title = paste0("Forest (comparatif) — ", mol),
           x = "log(OR) pooled", y = "Adverse event (AE)") +
      theme_bw() +
      theme(legend.position = "bottom", panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5))
    
    outfile <- file.path(outdir, paste0("forest_compare_", gsub("[^A-Za-z0-9]+","_", mol), ".pdf"))
    ggsave(outfile, p, width = 8.5, height = max(4.5, 0.25*length(levels(sumdf$ae_term)) + 2), dpi = 300)
    message("Saved: ", outfile, "  (", length(levels(sumdf$ae_term)), " AE)")
  }
  
  invisible(TRUE)
}

# ---- 1B. Figure composite avec légende commune (toutes molécules) ----
forest_compare_all_molecules_combined <- function(es_all,
                                                  outfile = "results_compare/forest_combined_all_molecules.pdf",
                                                  min_k_per_window = 2,
                                                  order_by = c("pval","|est|","name"),
                                                  star_offset = 0.06,
                                                  ncol = 2) {
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  order_by <- match.arg(order_by)
  
  es_all <- es_all %>%
    filter(is.finite(yi), is.finite(vi)) %>%
    mutate(time_window = factor(time_window, levels = c("session","follow_up")))
  
  mols  <- sort(unique(es_all$molecule))
  plist <- list()
  
  for (mol in mols) {
    dmol <- es_all %>% filter(molecule == !!mol)
    if (!nrow(dmol)) next
    
    sumdf <- .summarise_ae_by_window(dmol, min_k_per_window = min_k_per_window)
    if (is.null(sumdf) || !nrow(sumdf)) next
    
    sumdf <- sumdf %>%
      mutate(
        stars = case_when(
          is.na(pval) ~ "",
          pval < 0.001 ~ "***",
          pval < 0.01  ~ "**",
          pval < 0.05  ~ "*",
          TRUE ~ ""
        ),
        x_star = ci_ub + ifelse(time_window == "session",  star_offset, -star_offset)
      )
    
    levels_ae <- .order_ae(sumdf, order_by = order_by)
    sumdf <- sumdf %>% mutate(ae_term = factor(ae_term, levels = rev(levels_ae)))
    
    p <- ggplot(sumdf, aes(y = ae_term, x = est, xmin = ci_lb, xmax = ci_ub, color = time_window)) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.45) +
      geom_errorbarh(height = 0.22, position = position_dodge(width = 0.6)) +
      geom_point(size = 2.3, position = position_dodge(width = 0.6)) +
      geom_text(aes(x = x_star, label = stars), color = "red3", size = 4, show.legend = FALSE) +
      scale_color_manual(values = c(session = "#1f77b4", follow_up = "#ff7f0e"),
                         name = "Window", labels = c("Session","Follow-up")) +
      labs(title = mol, x = "log(OR)", y = "AE") +
      theme_bw(base_size = 10) +
      theme(legend.position = "bottom",
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5))
    plist[[mol]] <- p
  }
  
  if (!length(plist)) { warning("Aucune figure à combiner."); return(invisible(NULL)) }
  
  combined <- wrap_plots(plist, ncol = ncol, guides = "collect") & theme(legend.position = "bottom")
  ggsave(outfile, combined, width = 14, height = 10, dpi = 300)
  message("✅ Combined forest plots saved to: ", outfile)
  invisible(combined)
}