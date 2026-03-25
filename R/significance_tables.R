suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr) })

normalize_model_results <- function(df){
  if (is.null(df) || !nrow(df)) return(df)

  needed <- c("df_spline","term","estimate","beta","ci_low","ci_high","z","tau2","QE","I2","QM","QMp")
  for (nm in needed) {
    if (!nm %in% names(df)) df[[nm]] <- NA_real_
  }
  if (!"model" %in% names(df)) df$model <- NA_character_

  df %>%
    mutate(
      model = dplyr::case_when(
        !is.na(df_spline) & grepl("^spline$", model, ignore.case = TRUE) ~ paste0("spline_df", as.integer(df_spline)),
        TRUE ~ as.character(model)
      ),
      term = dplyr::coalesce(as.character(term), "dose_diff"),
      estimate = dplyr::coalesce(as.numeric(estimate), as.numeric(beta))
    )
}

save_significance_table <- function(results, outfile){
  if (is.null(results) || !nrow(results)) { warning("No results to save for: ", outfile); return(invisible(NULL)) }
  tab <- normalize_model_results(results) %>%
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
  tab <- normalize_model_results(models_by_molecule) %>%
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
  tab <- normalize_model_results(models_by_ae) %>%
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
