suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(ggrepel) })

# 1) Par molécule – tous AE confondus
plot_dr_by_molecule_split <- function(preds_by_molecule, models_by_molecule,
                                      outdir = "results/dose_response/by_molecule_split"){
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(preds_by_molecule) || !nrow(preds_by_molecule)) return(invisible(NULL))
  mols <- sort(unique(preds_by_molecule$molecule))
  for (mol in mols){
    pdat <- preds_by_molecule %>% filter(molecule == mol)
    if (!nrow(pdat)) next
    pref_model <- if (any(pdat$model == "spline_df3")) "spline_df3" else "linear"
    pdat <- pdat %>% filter(model == pref_model)
    
    mm <- models_by_molecule %>% filter(molecule == mol, model == pref_model)
    p_overall <- suppressWarnings(if (any(is.finite(mm$QMp))) min(mm$QMp, na.rm=TRUE) else
      min(mm$pval[!grepl("intrcpt", mm$term, TRUE)], na.rm=TRUE))
    stars <- if (is.finite(p_overall)) if (p_overall<0.001) "***" else if (p_overall<0.01) "**" else if (p_overall<0.05) "*" else "" else ""
    
    p <- ggplot(pdat, aes(dose_mg, fit)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "grey85") +
      geom_line(linewidth = 1, color = "#1f77b4") +
      labs(x = "Dose (mg)", y = "Log(OR) vs ref", title = paste0("Dose–response — ", mol, " ", stars),
           subtitle = "All AEs, all time windows") +
      theme_bw()
    ggsave(file.path(outdir, paste0("dr_", gsub("[^A-Za-z0-9._-]+","_", mol), ".png")),
           p, width = 8, height = 5.5, dpi = 150)
  }
  invisible(NULL)
}

# 2) Par molécule – facettes par AE
plot_dr_per_molecule_across_ae_facets <- function(preds_by_ae, models_by_ae,
                                                  outdir = "results/dose_response/by_molecule_facets"){
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (!nrow(preds_by_ae)) return(invisible(NULL))
  
  choose_model <- preds_by_ae %>%
    group_by(molecule, ae_term) %>% summarise(use_spline = any(model == "spline_df3"), .groups = "drop")
  df_all <- preds_by_ae %>%
    left_join(choose_model, by = c("molecule","ae_term")) %>%
    filter((use_spline & model == "spline_df3") | (!use_spline & model == "linear"))
  
  mols <- sort(unique(df_all$molecule))
  for (mol in mols){
    df <- df_all %>% filter(molecule == mol)
    if (!nrow(df)) next
    star_per_ae <- models_by_ae %>%
      filter(molecule == mol) %>%
      group_by(ae_term) %>%
      summarise(p = {
        v <- unique(QMp[is.finite(QMp)])
        if (length(v)) min(v) else suppressWarnings(min(pval[!grepl("intrcpt", term, TRUE)], na.rm = TRUE))
      }, .groups = "drop") %>%
      mutate(stars = case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ ""))
    
    df <- df %>%
      left_join(star_per_ae, by = "ae_term") %>%
      mutate(ae_lab = ifelse(is.na(stars) | stars=="", ae_term, paste0(ae_term, " ", stars)))
    
    p <- ggplot(df, aes(dose_mg, fit)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "grey90") +
      geom_line(linewidth = 0.9, color = "#1f77b4") +
      facet_wrap(~ ae_lab, scales = "free_y") +
      labs(x = "Dose (mg)", y = "Log(OR) vs ref",
           title = paste0("Dose–response — ", mol, " (facets by AE)")) +
      theme_bw()
    ggsave(file.path(outdir, paste0("dr_", gsub("[^A-Za-z0-9._-]+","_", mol), "_byAE_facets.png")),
           p, width = 12, height = 8, dpi = 150)
  }
  invisible(NULL)
}

# 3) Par AE – superposition molécules avec dose normalisée 0–100%
plot_dr_per_ae_normalized_dose <- function(preds_by_ae, outdir = "results/dose_response/by_ae_normalized"){
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (!nrow(preds_by_ae)) return(invisible(NULL))
  
  choose_model <- preds_by_ae %>%
    group_by(ae_term, molecule) %>% summarise(use_spline = any(model == "spline_df3"), .groups = "drop")
  df <- preds_by_ae %>%
    left_join(choose_model, by = c("ae_term","molecule")) %>%
    filter((use_spline & model == "spline_df3") | (!use_spline & model == "linear"))
  
  maxdose <- df %>% group_by(molecule) %>% summarise(max_dose = max(dose_mg, na.rm = TRUE), .groups = "drop")
  df <- df %>%
    left_join(maxdose, by = "molecule") %>%
    mutate(dose_pct = ifelse(is.finite(max_dose) & max_dose > 0, 100*dose_mg/max_dose, NA_real_)) %>%
    filter(is.finite(dose_pct))
  
  for (ae in sort(unique(df$ae_term))){
    d <- df %>% filter(ae_term == ae)
    if (!nrow(d)) next
    p <- ggplot(d, aes(dose_pct, fit, color = molecule, fill = molecule)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.12, color = NA) +
      geom_line(linewidth = 1) +
      scale_x_continuous(labels = function(x) paste0(x, "%")) +
      labs(x = "Normalized dose (% of molecule max)", y = "Log(OR) vs ref",
           title = paste0("AE: ", ae, " — molecules overlay")) +
      theme_bw()
    ggsave(file.path(outdir, paste0("dr_ae_", gsub("[^A-Za-z0-9._-]+","_", ae), "_doseNorm.png")),
           p, width = 10, height = 6, dpi = 150)
  }
  invisible(NULL)
}

# --- Diagnostics LSD (optionnels) ---
plot_lsd_ref_to_active <- function(contr, outdir = "results/diagnostics/lsd"){
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  df <- contr %>%
    filter(molecule == "LSD") %>%
    distinct(study_id, ref_arm_type, ref_dose_mg, dose_mg) %>%
    arrange(study_id, ref_dose_mg, dose_mg) %>%
    mutate(y = as.integer(factor(paste(study_id, round(dose_mg,3)))))
  if (!nrow(df)) return(invisible(NULL))
  p <- ggplot(df) +
    geom_segment(aes(x = ref_dose_mg, xend = dose_mg, y = y, yend = y),
                 arrow = arrow(length = unit(0.15,"cm")), linewidth = 0.6, alpha = 0.8) +
    geom_point(aes(x = ref_dose_mg, y = y, color = ref_arm_type), size = 2) +
    geom_point(aes(x = dose_mg,     y = y), size = 2, shape = 21, fill = "white") +
    scale_y_continuous(breaks = unique(df$y),
                       labels = df %>% group_by(y) %>% summarise(lab = first(paste0(study_id," | dose=",dose_mg)), .groups="drop") %>% pull(lab)) +
    labs(x = "Dose LSD (mg) — absolute", y = "Study | active dose",
         title = "LSD — Reference → Active dose per study",
         subtitle = "Left point = reference (inactive/active placebo or lowest dose)") +
    theme_bw() + theme(axis.text.y = element_text(size = 8))
  ggsave(file.path(outdir, "lsd_ref_to_active_segments.png"), p, width = 11, height = max(5, nrow(df)*0.18), dpi = 150)
  invisible(p)
}

plot_lsd_abs_vs_diff <- function(es, outdir = "results/diagnostics/lsd"){
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  df <- es %>%
    filter(molecule == "LSD") %>%
    distinct(study_id, dose_mg, ref_dose_mg, dose_diff)
  if (!nrow(df)) return(invisible(NULL))
  p <- ggplot(df, aes(dose_mg, dose_diff, label = study_id)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(size = 2) +
    ggrepel::geom_text_repel(min.segment.length = 0, size = 3, seed = 1) +
    labs(x = "Dose LSD (mg) — absolute", y = "Dose difference (dose - ref)",
         title = "LSD — absolute dose vs dose difference",
         subtitle = "dose_diff aligns studies with active placebo or no true placebo") +
    theme_bw()
  ggsave(file.path(outdir, "lsd_abs_vs_diff_scatter.png"), p, width = 8, height = 5.5, dpi = 150)
  invisible(p)
}
