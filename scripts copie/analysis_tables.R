<<<<<<< Updated upstream
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr) })

save_significance_table <- function(results, outfile){
  if (is.null(results) || !nrow(results)) { warning("No results to save for: ", outfile); return(invisible(NULL)) }
  tab <- results %>%
    mutate(significance = dplyr::case_when(
      is.na(pval)        ~ "",
      pval < 0.001       ~ "***",
      pval < 0.01        ~ "**",
      pval < 0.05        ~ "*",
      TRUE               ~ ""
    )) %>%
    select(
      molecule,
      dplyr::any_of("ae_term"),
      dplyr::any_of("time_window"),
      model, term, estimate, ci_low, ci_high, z, pval, I2, tau2, k, QE,
      dplyr::any_of(c("QM","QMp")),
      significance
    )
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(tab, outfile)
  invisible(tab)
}

save_agg_significance_by_molecule <- function(models_by_molecule, outfile){
  if (is.null(models_by_molecule) || !nrow(models_by_molecule)) return(invisible(NULL))
  tab <- models_by_molecule %>%
    filter(model %in% c("spline_df3","linear")) %>%
    group_by(molecule) %>%
    summarise(
      k_total   = suppressWarnings(max(k, na.rm = TRUE)),
      QM        = suppressWarnings(min(QM,  na.rm = TRUE)),
      QMp       = suppressWarnings(min(QMp, na.rm = TRUE)),
      p_overall = ifelse(is.finite(QMp), QMp,
                         suppressWarnings(min(pval[!grepl("intrcpt", term, TRUE)], na.rm = TRUE))),
      stars     = dplyr::case_when(
        is.na(p_overall)   ~ "",
        p_overall < 0.001  ~ "***",
        p_overall < 0.01   ~ "**",
        p_overall < 0.05   ~ "*",
        TRUE               ~ ""
      ),
      .groups   = "drop"
    )
  readr::write_csv(tab, outfile); invisible(tab)
}

save_agg_significance_by_ae_molecule <- function(models_by_ae, outfile){
  if (is.null(models_by_ae) || !nrow(models_by_ae)) return(invisible(NULL))
  tab <- models_by_ae %>%
    filter(model %in% c("spline_df3","linear")) %>%
    group_by(ae_term, molecule) %>%
    summarise(
      k_total   = suppressWarnings(max(k, na.rm = TRUE)),
      QM        = suppressWarnings(min(QM,  na.rm = TRUE)),
      QMp       = suppressWarnings(min(QMp, na.rm = TRUE)),
      p_overall = ifelse(is.finite(QMp), QMp,
                         suppressWarnings(min(pval[!grepl("intrcpt", term, TRUE)], na.rm = TRUE))),
      stars     = dplyr::case_when(
        is.na(p_overall)   ~ "",
        p_overall < 0.001  ~ "***",
        p_overall < 0.01   ~ "**",
        p_overall < 0.05   ~ "*",
        TRUE               ~ ""
      ),
      .groups   = "drop"
    )
  readr::write_csv(tab, outfile); invisible(tab)
}
=======
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

add_sig_stars <- function(p){
  dplyr::case_when(
    is.na(p)           ~ "",
    p < 0.001          ~ "***",
    p < 0.01           ~ "**",
    p < 0.05           ~ "*",
    TRUE               ~ ""
  )
}

save_significance_tables <- function(results, outfile = "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/results/tables/significance_results.csv"){
  if (!nrow(results)) {
    warning("No results to save.")
    return(invisible(NULL))
  }
  tab <- results |>
    mutate(significance = add_sig_stars(pval)) |>
    select(molecule, time_window, ae_term, model, term,
           estimate, ci_low, ci_high, z, pval, I2, tau2, k, QE, significance)
  
  readr::write_csv(tab, outfile)
  invisible(tab)
}
>>>>>>> Stashed changes
