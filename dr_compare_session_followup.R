# scripts/dr_compare_session_followup.R
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

# Détecte la colonne dose à utiliser
.detect_dose_col <- function(df){
  if ("dose_eff" %in% names(df)) return("dose_eff")
  if ("dose_mg"  %in% names(df)) return("dose_mg")
  stop("No dose column found (expected 'dose_eff' or 'dose_mg').")
}

# Compare session vs follow_up en UNE figure (facets par molécule), avec IC
dr_compare_all_molecules <- function(preds_session, preds_followup,
                                     outfile = "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/results_compare/dr_by_molecule_session_vs_followup.pdf",
                                     scales = "free_x") {
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  
  ps <- preds_session  %>% mutate(time_window = "session")
  pf <- preds_followup %>% mutate(time_window = "follow_up")
  pdat <- bind_rows(ps, pf)
  
  stopifnot(all(c("molecule","pred","lwr","upr") %in% names(pdat)))
  dose_col <- .detect_dose_col(pdat)
  
  p <- ggplot(pdat, aes(x = .data[[dose_col]], y = pred, color = time_window, fill = time_window)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.45) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.12, color = NA, show.legend = FALSE) +
    geom_line(size = 1) +
    facet_wrap(~ molecule, scales = scales) +
    scale_color_manual(values = c(session = "#1f77b4", follow_up = "#ff7f0e"),
                       name = "Window", labels = c("Session","Follow-up")) +
    scale_fill_manual(values  = c(session = "#1f77b4", follow_up = "#ff7f0e"),
                      guide = "none") +
    labs(
      title = "Dose–response: session vs follow-up (par molécule)",
      x = if (dose_col=="dose_eff") "Dose difference (mg)" else "Dose (mg)",
      y = "log(OR) vs référence"
    ) +
    theme_bw() +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white"),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  ggsave(outfile, p, width = 12, height = 8, dpi = 300)
  message("✅ Saved: ", outfile)
  invisible(p)
}