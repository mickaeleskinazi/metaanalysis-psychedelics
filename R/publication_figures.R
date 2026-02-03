# =========================
# FILE: R/publication_figures.R
# =========================
# Drop-in replacement (publication-grade, consistent across Fig1/Fig2/Fig3)
# - White background everywhere (theme_classic)
# - No titles/subtitles (captions handled in manuscript)
# - High-quality export (PDF vector, dpi=600, bg=white)
# - time_window normalized to session / follow_up everywhere
# - Fig2: include ALL significant AEs per molecule (p<0.05), with optional per-molecule fallback fill
# - Fig2: WINDOW COLORS IN COLOR (not black/grey), legend simplified
#
# NOTE: This file is self-contained. Replace your current R/publication_figures.R with this.

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(scales)
  library(metafor)
  library(purrr)
  library(tibble)
  library(tidyr)
  library(readr)
})

# -------------------------
# HELPERS
# -------------------------
safe_dir <- function(path){
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

# Publication theme: white background, strong axes, no grids
paper_theme_publication <- function(base_size = 12){
  theme_classic(base_size = base_size) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      
      axis.line  = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      
      legend.position = "bottom",
      legend.title    = element_text(face = "bold"),
      legend.key      = element_rect(fill = "white", color = NA),
      
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold"),
      
      plot.title    = element_blank(),
      plot.subtitle = element_blank(),
      
      plot.margin = margin(12, 40, 12, 12)
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

# Session vs follow-up — COLORED (publication-friendly)
# (Okabe-Ito-ish: good contrast, print-safe enough)
window_cols <- function(){
  c(
    "session"   = "#0072B2",
    "follow_up" = "#D55E00"
  )
}

p_to_stars <- function(p){
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

# Normalize any time_window to "session" / "follow_up"
.normalize_window <- function(x){
  x <- as.character(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- sub("^_", "", x)
  x <- sub("_$", "", x)
  x <- gsub("^followup$", "follow_up", x)
  x <- gsub("^follow_up$", "follow_up", x)
  x
}

make_plot_dose <- function(df, dose_col = "dose_mg"){
  df %>%
    mutate(
      molecule  = toupper(as.character(molecule)),
      dose_plot = .data[[dose_col]],
      dose_unit = ifelse(molecule == "LSD", "µg", "mg")
    )
}

# High-quality save: PDF vector + 600 dpi + white bg
save_pub <- function(outfile, plot, width, height){
  safe_dir(dirname(outfile))
  ggsave(
    filename = outfile,
    plot     = plot,
    width    = width,
    height   = height,
    units    = "in",
    dpi      = 600,
    bg       = "white"
  )
  invisible(TRUE)
}

.fmt_p <- function(p){
  ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))
}

.sig_stars_chr <- function(p){
  out <- rep("", length(p))
  out[!is.na(p) & p < 0.001] <- "***"
  out[!is.na(p) & p >= 0.001 & p < 0.01] <- "**"
  out[!is.na(p) & p >= 0.01  & p < 0.05] <- "*"
  out
}

# -------------------------
# FIG 1 — Global DR session+follow-up (per model)
# -------------------------
plot_fig1_global_dr_session_followup <- function(
    preds_session,
    preds_followup,
    outfile,
    molecules = c("PSILOCYBIN","MDMA","LSD"),
    xlab = "Dose",
    ylab = "Log odds ratio (adverse events)"
){
  if (is.null(preds_session) || !nrow(preds_session)) return(invisible(NULL))
  
  df_s <- preds_session %>%
    mutate(molecule = toupper(as.character(molecule))) %>%
    filter(molecule %in% molecules) %>%
    mutate(window = "session")
  
  df_f <- tibble()
  if (!is.null(preds_followup) && nrow(preds_followup)) {
    df_f <- preds_followup %>%
      mutate(molecule = toupper(as.character(molecule))) %>%
      filter(molecule %in% c("LSD","MDMA")) %>%
      mutate(window = "follow_up")
  }
  
  df <- bind_rows(df_s, df_f) %>%
    mutate(
      molecule = factor(molecule, levels = molecules),
      window   = factor(window, levels = c("session","follow_up"))
    ) %>%
    make_plot_dose("dose_mg")
  
  if (!nrow(df)) return(invisible(NULL))
  
  p <- ggplot(df, aes(x = dose_plot, y = fit, color = molecule, fill = molecule, linetype = window)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.12, linewidth = 0) +
    geom_line(linewidth = 1.15) +
    facet_wrap(~ molecule, nrow = 1, scales = "free_x") +
    scale_color_manual(values = mol_colors()) +
    scale_fill_manual(values = mol_colors()) +
    scale_linetype_manual(values = c("session" = "solid", "follow_up" = "22")) +
    labs(x = xlab, y = ylab) +
    paper_theme_publication(base_size = 12) +
    guides(
      fill = "none",
      linetype = guide_legend(title = "Window", nrow = 1),
      color = guide_legend(title = "Molecule", nrow = 1)
    )
  
  save_pub(outfile, p, width = 12.0, height = 4.6)
  invisible(p)
}

# -------------------------
# FIG 2 — Forest plot overlay (session vs follow-up)
# -------------------------
summarise_or_by_ae_molecule_window <- function(es, min_k = 2){
  es2 <- es %>%
    mutate(
      molecule    = toupper(as.character(molecule)),
      time_window = .normalize_window(time_window)
    ) %>%
    filter(is.finite(yi), is.finite(vi)) %>%
    group_by(time_window, molecule, ae_term) %>%
    filter(n() >= min_k) %>%
    ungroup()
  
  if (!nrow(es2)) return(tibble())
  
  es2 %>%
    group_by(time_window, molecule, ae_term) %>%
    group_modify(~{
      m <- tryCatch(metafor::rma(yi, vi, data = .x, method = "REML"),
                    error = function(e) NULL)
      if (is.null(m)) return(tibble())
      tibble(
        k       = m$k,
        yi_hat  = as.numeric(m$b),
        ci_low  = as.numeric(m$ci.lb),
        ci_high = as.numeric(m$ci.ub),
        pval    = as.numeric(m$pval),
        OR      = exp(yi_hat),
        OR_low  = exp(ci_low),
        OR_high = exp(ci_high),
        stars   = p_to_stars(pval)
      )
    }) %>%
    ungroup()
}

plot_fig2_forest_main_aes_by_molecule_both_windows <- function(
    es_both,
    outfile,
    min_k = 2,
    molecules = c("MDMA","LSD","PSILOCYBIN","AYAHUASCA"),
    windows   = c("session","follow_up"),
    include_all_significant = TRUE,  # include ALL p<0.05 AEs per molecule
    fill_if_none = TRUE,             # if a molecule has 0 significant, show top-k by k_total
    fill_n = 8,                      # top-k count used if none significant for that molecule
    xlab = "Odds ratio (log scale)"
){
  sum_df <- summarise_or_by_ae_molecule_window(es_both, min_k = min_k) %>%
    filter(molecule %in% molecules, time_window %in% windows)
  
  if (!nrow(sum_df)) return(invisible(NULL))
  
  # -------- AE selection --------
  if (isTRUE(include_all_significant)) {
    
    keep_df <- sum_df %>%
      group_by(molecule) %>%
      summarise(
        ae_sig = list(unique(ae_term[!is.na(pval) & pval < 0.05])),
        .groups = "drop"
      )
    
    if (isTRUE(fill_if_none)) {
      
      topk_df <- sum_df %>%
        group_by(molecule, ae_term) %>%
        summarise(k_total = sum(k, na.rm = TRUE), .groups = "drop") %>%
        arrange(molecule, desc(k_total), ae_term) %>%
        group_by(molecule) %>%
        summarise(ae_topk = list(head(ae_term, fill_n)), .groups = "drop")
      
      keep_df <- keep_df %>%
        left_join(topk_df, by = "molecule") %>%
        mutate(
          ae_keep = purrr::map2(ae_sig, ae_topk, ~{
            sig <- .x
            top <- .y
            if (is.null(sig) || length(sig) == 0) top else sig
          })
        )
      
    } else {
      keep_df <- keep_df %>% mutate(ae_keep = ae_sig)
    }
    
    ae_keep <- keep_df %>% pull(ae_keep) %>% unlist() %>% unique()
    
  } else {
    top_n_ae <- 12
    ae_sig <- sum_df %>%
      filter(!is.na(pval), pval < 0.05) %>%
      distinct(ae_term) %>%
      pull(ae_term)
    
    ae_topk <- sum_df %>%
      group_by(ae_term) %>%
      summarise(k_total = sum(k, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(k_total)) %>%
      pull(ae_term)
    
    ae_keep <- ae_sig
    if (length(ae_keep) < top_n_ae) {
      ae_fill <- setdiff(ae_topk, ae_keep)
      ae_keep <- c(ae_keep, head(ae_fill, top_n_ae - length(ae_keep)))
    } else {
      ae_keep <- ae_keep[seq_len(top_n_ae)]
    }
  }
  
  df <- sum_df %>%
    filter(ae_term %in% ae_keep) %>%
    mutate(
      molecule    = factor(molecule, levels = molecules),
      time_window = factor(time_window, levels = windows)
    )
  
  # Order AEs within each facet: significant first, then by k_total
  df <- df %>%
    group_by(molecule, ae_term) %>%
    mutate(
      k_total    = sum(k, na.rm = TRUE),
      is_sig_any = any(!is.na(pval) & pval < 0.05)
    ) %>%
    ungroup() %>%
    arrange(molecule, desc(is_sig_any), desc(k_total), ae_term) %>%
    group_by(molecule) %>%
    mutate(ae_term = factor(ae_term, levels = unique(ae_term))) %>%
    ungroup()
  
  pd <- position_dodge(width = 0.60)
  win_cols <- window_cols()
  
  p <- ggplot(df, aes(y = ae_term, x = OR, xmin = OR_low, xmax = OR_high)) +
    geom_vline(xintercept = 1, linewidth = 0.6, color = "black") +
    geom_errorbarh(aes(color = time_window), height = 0.22, linewidth = 1.0, position = pd) +
    geom_point(aes(shape = time_window, color = time_window), size = 2.8, stroke = 0.4, position = pd) +
    geom_text(aes(label = stars, color = time_window),
              size = 3.8, position = pd, hjust = -0.35, show.legend = FALSE) +
    facet_wrap(~ molecule, nrow = 1, scales = "free_y") +
    scale_color_manual(values = win_cols) +
    scale_shape_manual(values = c("session" = 16, "follow_up" = 17)) +
    scale_x_log10(labels = scales::label_number(accuracy = 0.01)) +
    coord_cartesian(clip = "off") +
    labs(x = xlab, y = NULL) +
    paper_theme_publication(base_size = 12) +
    guides(
      color = guide_legend(title = "Window", nrow = 1),
      shape = "none"
    )
  
  n_ae <- length(unique(df$ae_term))
  height <- max(6.5, 0.22 * n_ae + 2.6)
  
  save_pub(outfile, p, width = 13.0, height = height)
  invisible(p)
}

# -------------------------
# FIG 3 — Shared AEs DR (normalized dose)
# -------------------------
summarise_dr_slope_by_ae_molecule_linear <- function(es, min_k = 4){
  if (is.null(es) || !nrow(es)) return(tibble())
  
  es2 <- es %>%
    filter(is.finite(yi), is.finite(vi), is.finite(dose_mg)) %>%
    mutate(molecule = toupper(as.character(molecule))) %>%
    group_by(molecule, ae_term) %>%
    filter(n() >= min_k, n_distinct(dose_mg) >= 2) %>%
    ungroup()
  
  if (!nrow(es2)) return(tibble())
  
  groups <- es2 %>% group_by(molecule, ae_term) %>% group_split()
  
  purrr::map_dfr(groups, function(g){
    m <- tryCatch(metafor::rma(yi, vi, mods = ~ dose_mg, data = g, method = "REML"),
                  error = function(e) NULL)
    if (is.null(m)) return(tibble())
    p <- tryCatch(as.numeric(m$QMp), error = function(e) NA_real_)
    tibble(
      molecule = g$molecule[[1]],
      ae_term  = g$ae_term[[1]],
      p_slope  = p
    )
  })
}

plot_fig3_dr_shared_aes_normalized <- function(
    preds_ae,
    es,
    outfile,
    min_shared_molecules = 2,
    top_n_ae = 9,
    min_k_slope = 4,
    molecules = c("MDMA","LSD","PSILOCYBIN","AYAHUASCA"),
    xlab = "Normalized dose within molecule (0–1)",
    ylab = "Log odds ratio (adverse event)"
){
  if (is.null(preds_ae) || !nrow(preds_ae)) return(invisible(NULL))
  if (is.null(es) || !nrow(es)) return(invisible(NULL))
  
  preds_ae <- preds_ae %>%
    mutate(molecule = toupper(as.character(molecule))) %>%
    filter(molecule %in% molecules)
  
  if (!nrow(preds_ae)) return(invisible(NULL))
  
  shared_ae <- preds_ae %>%
    distinct(ae_term, molecule) %>%
    count(ae_term, name = "n_mol") %>%
    filter(n_mol >= min_shared_molecules) %>%
    pull(ae_term)
  
  if (!length(shared_ae)) return(invisible(NULL))
  
  ae_keep <- preds_ae %>%
    filter(ae_term %in% shared_ae) %>%
    count(ae_term, name = "n_points") %>%
    arrange(desc(n_points)) %>%
    slice_head(n = top_n_ae) %>%
    pull(ae_term)
  
  df <- preds_ae %>%
    filter(ae_term %in% ae_keep) %>%
    make_plot_dose("dose_mg") %>%
    group_by(molecule) %>%
    mutate(
      mx = max(dose_plot, na.rm = TRUE),
      dose_norm = if_else(is.finite(mx) & mx > 0, dose_plot / mx, NA_real_)
    ) %>%
    ungroup() %>%
    filter(is.finite(dose_norm)) %>%
    mutate(
      molecule = factor(molecule, levels = molecules),
      ae_term  = factor(ae_term, levels = ae_keep)
    )
  
  if (!nrow(df)) return(invisible(NULL))
  
  sig_df <- summarise_dr_slope_by_ae_molecule_linear(es, min_k = min_k_slope) %>%
    mutate(
      molecule  = toupper(as.character(molecule)),
      sig_label = if_else(!is.na(p_slope) & p_slope < 0.05, "Significant", "Not significant")
    ) %>%
    filter(molecule %in% molecules, ae_term %in% ae_keep) %>%
    select(molecule, ae_term, sig_label)
  
  df <- df %>%
    left_join(sig_df, by = c("molecule","ae_term")) %>%
    mutate(
      sig_label = if_else(is.na(sig_label), "Not significant", sig_label),
      sig_label = factor(sig_label, levels = c("Significant", "Not significant"))
    )
  
  p <- ggplot(
    df,
    aes(
      x = dose_norm, y = fit,
      color = molecule,
      linetype = sig_label,
      group = interaction(molecule, ae_term)
    )
  ) +
    geom_line(linewidth = 1.15) +
    facet_wrap(~ ae_term, ncol = 3) +
    scale_color_manual(values = mol_colors()) +
    scale_linetype_manual(values = c("Significant" = "solid", "Not significant" = "22")) +
    scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
    labs(x = xlab, y = ylab) +
    paper_theme_publication(base_size = 11) +
    guides(
      color    = guide_legend(title = "Molecule", nrow = 1, override.aes = list(linetype = "solid", linewidth = 2)),
      linetype = guide_legend(title = "Dose–slope", nrow = 1, override.aes = list(color = "black", linewidth = 2))
    )
  
  save_pub(outfile, p, width = 12.2, height = 9.2)
  invisible(p)
}

# -------------------------
# WRAPPERS CALLED BY run_main_analysis.R
# -------------------------
make_paper_figures_for_window <- function(
    window_value,
    es,
    dr_mol_plot,
    dr_ae_plot,
    out_dir,
    model,
    df_spline = 3,
    only_session = TRUE,
    ...
){
  if (isTRUE(only_session) && window_value != "session") return(invisible(TRUE))
  
  fig_dir <- file.path(out_dir, "paper_figures")
  safe_dir(fig_dir)
  
  tag <- ifelse(is.null(model), "model", model)
  
  # Fig3 (session only)
  if (!is.null(dr_ae_plot) && !is.null(dr_ae_plot$preds) && nrow(dr_ae_plot$preds)) {
    plot_fig3_dr_shared_aes_normalized(
      preds_ae = dr_ae_plot$preds,
      es       = es,
      outfile  = file.path(fig_dir, paste0("Fig3_DR_sharedAEs_", tag, ".pdf")),
      top_n_ae = 9
    )
  }
  
  invisible(TRUE)
}

make_paper_figures_comparison <- function(results, out_dir, df_spline = 3){
  if (is.null(results$session) || is.null(results$follow_up)) return(invisible(NULL))
  
  model_tags <- intersect(names(results$session$dr), names(results$follow_up$dr))
  if (!length(model_tags)) return(invisible(NULL))
  
  es_both <- bind_rows(results$session$es, results$follow_up$es) %>%
    mutate(time_window = .normalize_window(time_window))
  
  purrr::walk(model_tags, function(tag){
    
    fig_dir2 <- file.path(out_dir, "compare", "paper_figures_mixed")
    safe_dir(fig_dir2)
    
    # Fig1 mixed: session spline + follow-up (LSD spline, MDMA linear)
    preds_session_mix <- results$session$dr$spline$mol_plot$preds
    preds_followup_mix <- dplyr::bind_rows(
      results$follow_up$dr$spline$mol_plot$preds %>% dplyr::filter(toupper(molecule) == "LSD"),
      results$follow_up$dr$linear$mol_plot$preds %>% dplyr::filter(toupper(molecule) == "MDMA")
    )
    
    plot_fig1_global_dr_session_followup(
      preds_session  = preds_session_mix,
      preds_followup = preds_followup_mix,
      outfile        = file.path(fig_dir2, "Fig1_Global_DR_session_spline_followup_mixed.pdf")
    )
    
    # Fig2 forest: include ALL significant AEs per molecule
    plot_fig2_forest_main_aes_by_molecule_both_windows(
      es_both  = es_both,
      outfile  = file.path(fig_dir2, paste0("Fig2_Forest_mainAEs_session_followup_", tag, ".pdf")),
      include_all_significant = TRUE,
      fill_if_none = TRUE,
      fill_n = 8
    )
  })
  
  invisible(TRUE)
}

# -------------------------
# FOREST TABLE EXPORTS (publication tables)
# -------------------------
make_forest_tables <- function(es,
                               out_dir = file.path("results", "forest_tables"),
                               min_k = 2){
  
  required_cols <- c("yi", "vi", "molecule", "ae_term", "time_window")
  if (!all(required_cols %in% names(es))) {
    stop("Effect size table is missing required columns: ",
         paste(setdiff(required_cols, names(es)), collapse = ", "))
  }
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (!nrow(es)) {
    empty_long <- tibble::tibble(
      molecule = character(),
      time_window = character(),
      ae_term = character(),
      k = integer(),
      log_or = double(),
      se = double(),
      z = double(),
      pval = double(),
      ci_lb = double(),
      ci_ub = double(),
      or = double(),
      or_ci_lb = double(),
      or_ci_ub = double(),
      tau2 = double(),
      I2 = double(),
      QE = double(),
      QEp = double(),
      stars = character(),
      log_or_ci = character(),
      or_ci = character(),
      pval_fmt = character(),
      QEp_fmt = character(),
      status = character()
    )
    readr::write_csv(empty_long, file.path(out_dir, "forest_by_molecule_ae_window.csv"))
    readr::write_csv(empty_long, file.path(out_dir, "forest_by_molecule_ae_window_publication.csv"))
    readr::write_csv(tibble::tibble(), file.path(out_dir, "forest_by_molecule_ae_window_wide.csv"))
    return(invisible(list(long = empty_long, wide = tibble::tibble(), publication = empty_long)))
  }
  
  es <- es %>%
    mutate(
      molecule    = toupper(as.character(molecule)),
      time_window = .normalize_window(time_window)
    )
  
  grouped <- es %>%
    filter(is.finite(yi), is.finite(vi)) %>%
    group_by(molecule, time_window, ae_term)
  
  pooled <- grouped %>%
    group_modify(~{
      dat <- .x
      k_obs <- nrow(dat)
      base <- tibble::tibble(
        k = as.integer(k_obs),
        log_or = NA_real_,
        se = NA_real_,
        z = NA_real_,
        pval = NA_real_,
        ci_lb = NA_real_,
        ci_ub = NA_real_,
        tau2 = NA_real_,
        I2 = NA_real_,
        QE = NA_real_,
        QEp = NA_real_
      )
      if (k_obs < min_k) return(base)
      
      fit <- tryCatch(
        metafor::rma(yi, vi, data = dat, method = "REML"),
        error = function(e) NULL
      )
      if (is.null(fit)) return(base)
      
      tibble::tibble(
        k = as.integer(fit$k),
        log_or = as.numeric(fit$b[1]),
        se = suppressWarnings(as.numeric(fit$se[1])),
        z = suppressWarnings(as.numeric(fit$zval[1])),
        pval = suppressWarnings(as.numeric(fit$pval[1])),
        ci_lb = as.numeric(fit$ci.lb),
        ci_ub = as.numeric(fit$ci.ub),
        tau2 = suppressWarnings(tryCatch(as.numeric(fit$tau2), error = function(e) NA_real_)),
        I2 = suppressWarnings(tryCatch(as.numeric(fit$I2), error = function(e) NA_real_)),
        QE = suppressWarnings(tryCatch(as.numeric(fit$QE), error = function(e) NA_real_)),
        QEp = suppressWarnings(tryCatch(as.numeric(fit$QEp), error = function(e) NA_real_))
      )
    }) %>%
    ungroup()
  
  pooled <- pooled %>%
    mutate(
      or = exp(log_or),
      or_ci_lb = exp(ci_lb),
      or_ci_ub = exp(ci_ub),
      stars = .sig_stars_chr(pval),
      pval_fmt = .fmt_p(pval),
      QEp_fmt = .fmt_p(QEp),
      status = dplyr::case_when(
        k < min_k ~ "insufficient_k",
        !is.na(log_or) ~ "ok",
        TRUE ~ "model_error"
      ),
      log_or_ci = dplyr::if_else(
        !is.na(log_or) & !is.na(ci_lb) & !is.na(ci_ub),
        sprintf("%.3f [%.3f, %.3f]%s", log_or, ci_lb, ci_ub,
                dplyr::if_else(stars == "", "", paste0(" ", stars))),
        ""
      ),
      or_ci = dplyr::if_else(
        !is.na(or) & !is.na(or_ci_lb) & !is.na(or_ci_ub),
        sprintf("%.2f [%.2f, %.2f]%s", or, or_ci_lb, or_ci_ub,
                dplyr::if_else(stars == "", "", paste0(" ", stars))),
        ""
      )
    ) %>%
    mutate(
      I2 = suppressWarnings(as.numeric(I2)),
      tau2 = suppressWarnings(as.numeric(tau2))
    ) %>%
    arrange(molecule, time_window, ae_term)
  
  readr::write_csv(pooled, file.path(out_dir, "forest_by_molecule_ae_window.csv"))
  
  wide <- pooled %>%
    tidyr::pivot_wider(
      id_cols = c(molecule, ae_term),
      names_from = time_window,
      values_from = c(log_or, ci_lb, ci_ub, or, or_ci_lb, or_ci_ub,
                      pval, pval_fmt, stars, k, I2, tau2, QE, QEp,
                      log_or_ci, or_ci, status),
      names_sep = "."
    ) %>%
    arrange(molecule, ae_term)
  
  readr::write_csv(wide, file.path(out_dir, "forest_by_molecule_ae_window_wide.csv"))
  
  publication <- pooled %>%
    transmute(
      Molecule = molecule,
      `Adverse event` = ae_term,
      Window = time_window,
      k = k,
      `OR [95% CI]` = or_ci,
      `log(OR) [95% CI]` = log_or_ci,
      `p` = pval_fmt,
      Stars = stars,
      `τ²` = tau2,
      `I²` = I2,
      `Q` = QE,
      `p(Q)` = QEp_fmt,
      Status = status
    )
  
  readr::write_csv(publication, file.path(out_dir, "forest_by_molecule_ae_window_publication.csv"))
  
  invisible(list(long = pooled, wide = wide, publication = publication))
}