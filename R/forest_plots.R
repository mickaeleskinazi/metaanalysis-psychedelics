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
