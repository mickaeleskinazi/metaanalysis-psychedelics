suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(metafor); library(purrr) })

# PNGs unitaires (optionnels, si tu as déjà une fonction existante)
make_forest_plots <- function(es, outdir = "results/forest_plots"){
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (!nrow(es)) return(invisible(NULL))
  # Exemple minimal : par (molecule, time_window, ae_term)
  groups <- es %>% group_by(molecule, time_window, ae_term) %>% group_split()
  for (g in groups){
    if (nrow(g) < 2) next
    m <- try(rma(yi, vi, data = g, method = "REML"), silent = TRUE)
    if (inherits(m,"try-error")) next
    fname <- file.path(outdir, paste0("forest_",
                                      g$molecule[[1]], "_", gsub("[^A-Za-z0-9]+","_", g$time_window[[1]]), "_",
                                      gsub("[^A-Za-z0-9]+","_", g$ae_term[[1]]), ".png"))
    png(fname, width = 960, height = 720)
    metafor::forest(m, xlab = "Log(OR)", slab = g$study_id)
    title(main = paste(g$molecule[[1]], "—", g$time_window[[1]], "—", g$ae_term[[1]]))
    dev.off()
  }
  invisible(NULL)
}

# Un PDF par molécule (une page par AE) — robuste
make_forest_plots_per_molecule_pdf <- function(es, outdir = "results/forest_plots_by_molecule"){
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (!nrow(es)) return(invisible(NULL))
  mols <- sort(unique(es$molecule))
  for (mol in mols){
    dmol <- es %>% filter(molecule == mol)
    aes_list <- unique(dmol$ae_term)
    f <- file.path(outdir, paste0("forest_", gsub("[^A-Za-z0-9._-]+","_", mol), ".pdf"))
    opened <- FALSE; pages <- 0L
    for (ae in aes_list){
      g <- dmol %>% filter(ae_term == ae)
      if (nrow(g) < 2) next
      m <- try(metafor::rma(yi, vi, data = g, method = "REML"), silent = TRUE)
      if (inherits(m, "try-error")) next
      if (!opened) { grDevices::pdf(f, width = 10, height = 7); opened <- TRUE }
      metafor::forest(m, slab = g$study_id, xlab = "Log(OR) for AE"); title(main = paste(mol, "–", ae)); pages <- pages + 1L
    }
    if (opened && pages > 0L) grDevices::dev.off()
    else if (opened && pages == 0L) { grDevices::dev.off(); if (file.exists(f) && file.info(f)$size == 0) unlink(f) }
  }
  invisible(NULL)
}

# Résumé par molécule : 1 ligne = 1 AE (toutes études/bras confondus)
make_forest_summary_per_molecule <- function(es, outdir = file.path("results", "forest_plots_summary")){
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (!nrow(es)) return(invisible(NULL))
  
  mols <- sort(unique(es$molecule))
  for (mol in mols){
    dmol <- es %>% filter(molecule == mol)
    if (!nrow(dmol)) next
    ae_summ <- dmol %>%
      group_by(ae_term) %>%
      group_map(~{
        g <- .x; if (nrow(g) < 2) return(NULL)
        m <- try(metafor::rma(yi, vi, data = g, method = "REML"), silent = TRUE)
        if (inherits(m,"try-error")) return(NULL)
        tibble(
          ae_term = g$ae_term[[1]],
          est     = as.numeric(m$b[1,1]),
          ci_lb   = as.numeric(m$ci.lb),
          ci_ub   = as.numeric(m$ci.ub),
          pval    = as.numeric(m$pval)
        )
      }, .keep = TRUE) %>% purrr::list_rbind()
    
    if (is.null(ae_summ) || !nrow(ae_summ)) next
    
    df <- ae_summ %>%
      mutate(stars = case_when(is.na(pval) ~ "", pval < 0.001 ~ "***", pval < 0.01 ~ "**", pval < 0.05 ~ "*", TRUE ~ "")) %>%
      arrange(est) %>% mutate(ae_term = factor(ae_term, levels = ae_term))
    
    p <- ggplot(df, aes(x = est, y = ae_term)) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_errorbarh(aes(xmin = ci_lb, xmax = ci_ub), height = 0.22) +
      geom_point(aes(color = stars != ""), size = 2) +
      scale_color_manual(values = c("FALSE"="black","TRUE"="#d62728"), guide = "none") +
      geom_text(aes(label = stars), nudge_x = 0.06, size = 4) +
      labs(x = "Log(OR) (pooled, all studies)", y = "AE", title = paste0("Forest (summary) — ", mol)) +
      theme_bw()
    ggsave(file.path(outdir, paste0("forest_summary_", gsub("[^A-Za-z0-9._-]+","_", mol), ".png")),
           p, width = 8, height = 10, dpi = 150)
  }
  invisible(NULL)
}

.normalize_window <- function(x){
  x <- as.character(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- sub("^_", "", x)
  sub("_$", "", x)
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

#' Create tabular exports of pooled forest statistics per molecule × AE × window
#'
#' @param es Effect-size data frame containing yi/vi columns and grouping variables.
#' @param out_dir Output directory where CSV files will be written.
#' @param min_k Minimum number of contrasts required to fit a pooled model.
#'
#' @return Invisibly returns a list with the long, wide, and publication-style tables.
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

  es <- es %>% mutate(time_window = .normalize_window(time_window))

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
      if (k_obs < min_k) {
        return(base)
      }
      fit <- tryCatch(
        metafor::rma(yi, vi, data = dat, method = "REML"),
        error = function(e) NULL
      )
      if (is.null(fit)) {
        return(base)
      }
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
    mutate(I2 = suppressWarnings(as.numeric(I2)),
           tau2 = suppressWarnings(as.numeric(tau2))) %>%
    arrange(molecule, time_window, ae_term)

  readr::write_csv(pooled, file.path(out_dir, "forest_by_molecule_ae_window.csv"))

  if (nrow(pooled)) {
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
  } else {
    wide <- tibble::tibble()
  }

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
