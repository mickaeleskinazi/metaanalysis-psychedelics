# =========================
# FILE: R/bubble_weight_plots.R
# =========================
suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(scales)
})

paper_theme <- function(base_size = 12){
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey30"),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", color = NA)
    )
}

mol_colors <- function(){
  c(
    "MDMA"       = "#C0392B",
    "LSD"        = "#2E86C1",
    "PSILOCYBIN" = "#239B56",
    "AYAHUASCA"  = "#8E44AD"
  )
}

# Bubble plot: study-level effect sizes with bubble size ~ inverse-variance weight
# Shows "which studies drive the signal" per molecule for the densest AEs
plot_bubble_weights_by_molecule <- function(
    es,
    outfile,
    molecules = c("MDMA","LSD","PSILOCYBIN","AYAHUASCA"),
    top_n_ae = 16,
    min_k_ae = 3,
    title = "Study weights by molecule (bubble plot)",
    subtitle = "Each point = study contrast; bubble size ∝ inverse-variance weight (1/vi). X axis = odds ratio (log scale).",
    xlab = "Odds ratio (log scale)",
    ylab = NULL
){
  if (is.null(es) || !nrow(es)) return(invisible(NULL))
  
  df <- es %>%
    filter(
      molecule %in% molecules,
      is.finite(yi), is.finite(vi), vi > 0,
      !is.na(ae_term), !is.na(study_id)
    ) %>%
    mutate(
      OR = exp(yi),
      w_raw = 1 / vi
    )
  
  if (!nrow(df)) return(invisible(NULL))
  
  # keep densest AEs overall (or per molecule if you prefer — here overall)
  ae_keep <- df %>%
    group_by(ae_term) %>%
    summarise(k_total = n(), .groups = "drop") %>%
    filter(k_total >= min_k_ae) %>%
    arrange(desc(k_total)) %>%
    slice_head(n = top_n_ae) %>%
    pull(ae_term)
  
  df <- df %>%
    filter(ae_term %in% ae_keep) %>%
    mutate(
      molecule = factor(molecule, levels = molecules),
      ae_term  = fct_reorder(ae_term, OR, .fun = median, na.rm = TRUE)
    )
  
  # normalize weights within each molecule (so sizes are comparable inside facet)
  df <- df %>%
    group_by(molecule) %>%
    mutate(
      w = w_raw / max(w_raw, na.rm = TRUE)  # 0–1
    ) %>%
    ungroup()
  
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  
  p <- ggplot(df, aes(x = OR, y = ae_term)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
    geom_point(aes(size = w, color = molecule), alpha = 0.70) +
    facet_wrap(~ molecule, nrow = 1, scales = "free_y") +
    scale_x_log10(labels = label_number(accuracy = 0.01)) +
    scale_color_manual(values = mol_colors()) +
    scale_size(range = c(1.2, 7.5), breaks = c(0.25, 0.5, 0.75, 1.0),
               labels = c("25%", "50%", "75%", "100%")) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab, size = "Relative weight") +
    paper_theme() +
    theme(legend.box = "vertical")
  
  ggsave(outfile, p, width = 12.8, height = 6.6, units = "in", dpi = 300)
  invisible(p)
}