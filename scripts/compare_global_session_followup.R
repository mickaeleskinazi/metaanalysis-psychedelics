#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(metafor)
  library(ggplot2)
})

.map_groups_dfr <- function(grouped_data, fn) {
  keys   <- dplyr::group_keys(grouped_data)
  splits <- dplyr::group_split(grouped_data, .keep = FALSE)
  purrr::imap_dfr(splits, function(dat, idx) {
    key <- keys[idx, , drop = FALSE]
    res <- fn(dat, key)
    if (!inherits(res, "data.frame")) {
      stop("Grouped mapping function must return a data frame.")
    }
    if (!nrow(res)) {
      return(res)
    }
    key_rep <- key[rep(1, nrow(res)), , drop = FALSE]
    dplyr::bind_cols(key_rep, res)
  })
}

source(here::here("R", "data_ingest.R"))
source(here::here("R", "dose_response_models.R"))

default_ref_policies <- list(
  MDMA       = c("inactive_placebo", "active_non_psy_placebo", "active_placebo"),
  LSD        = c("inactive_placebo", "active_placebo", "active_non_psy_placebo"),
  PSILOCYBIN = c("inactive_placebo", "active_non_psy_placebo", "active_placebo"),
  AYAHUASCA  = c("inactive_placebo", "active_placebo", "active_non_psy_placebo"),
  .default   = c("inactive_placebo", "active_placebo", "active_non_psy_placebo")
)

fmt_p <- function(p){
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "<0.001",
                formatC(p, format = "f", digits = 3)))
}

sig_stars_vec <- function(p){
  out <- rep("", length(p))
  out[!is.na(p) & p < 0.001] <- "***"
  out[!is.na(p) & p >= 0.001 & p < 0.01] <- "**"
  out[!is.na(p) & p >= 0.01  & p < 0.05] <- "*"
  out
}

build_es_all <- function(base_xlsx, sheet, ref_policies){
  suppressPackageStartupMessages({ library(readxl) })
  df <- readxl::read_excel(base_xlsx, sheet = sheet)

  req <- c("study_id","author_year","arm_id","molecule",
           "n_participants_arm","ae_term","time_window","events","dose_mg")
  miss <- setdiff(req, names(df))
  if (length(miss)) stop("Missing expected Excel columns: ", paste(miss, collapse=", "))

  raw <- df %>%
    transmute(
      study_id,
      author_year,
      molecule         = as.character(molecule),
      ae_term          = as.character(ae_term),
      dose_mg          = suppressWarnings(as.numeric(dose_mg)),
      group            = as.character(arm_id),
      arm_type         = case_when(
        str_detect(tolower(arm_id), "placebo") ~ "placebo",
        TRUE                                   ~ "active"
      ),
      n_total          = suppressWarnings(as.numeric(n_participants_arm)),
      n_events         = suppressWarnings(as.numeric(events)),
      time_window      = {
        tw <- tolower(trimws(as.character(time_window)))
        tw <- ifelse(str_detect(tw, "follow|follow[-_ ]?up|^fu$"), "follow_up", tw)
        tw <- ifelse(str_detect(tw, "session|acute|day 0|seance|séance"), "session", tw)
        forcats::fct_drop(factor(tw, levels = c("session","follow_up")))
      }
    )

  need <- c("study_id","molecule","ae_term","time_window","group","arm_type",
            "dose_mg","n_events","n_total")
  miss2 <- setdiff(need, names(raw))
  if (length(miss2)) stop("Harmonized data missing: ", paste(miss2, collapse=", "))

  raw <- raw %>%
    filter(!is.na(molecule), !is.na(ae_term), !is.na(dose_mg),
           !is.na(n_total), !is.na(n_events), !is.na(time_window))

  contr <- build_pairwise_2x2(raw, ref_policies = ref_policies)
  es    <- build_escalc(contr)

  es %>%
    filter(!is.na(time_window), is.finite(yi), is.finite(vi), !is.na(dose_mg))
}

compare_global_dose_slope_by_window <- function(es, min_k = 4) {
  stopifnot(all(c("molecule","dose_mg","time_window","yi","vi") %in% names(es)))
  grouped <- es %>%
    group_by(molecule)

  .map_groups_dfr(grouped, function(dat, key) {
    dat <- dat %>% filter(is.finite(dose_mg), is.finite(yi), is.finite(vi))
    dat <- dat %>% mutate(time_window = factor(time_window, levels = c("session","follow_up")))

    empty <- tibble(
      beta_session   = double(),
      beta_interact  = double(),
      p_session      = double(),
      p_interaction  = double(),
      k              = integer()
    )

    if (length(unique(dat$time_window)) < 2 || nrow(dat) < min_k) return(empty)

    mod <- tryCatch(
      metafor::rma(yi ~ dose_mg * time_window, vi = vi, data = dat, method = "REML"),
      error = function(e) NULL
    )

    if (is.null(mod)) return(empty)

    tibble(
      beta_session   = as.numeric(coef(mod)["dose_mg"]),
      beta_interact  = as.numeric(coef(mod)[grep("^dose_mg:time_window", names(coef(mod)))]),
      p_session      = as.numeric(mod$pval["dose_mg"]),
      p_interaction  = as.numeric(mod$pval[grep("^dose_mg:time_window", names(mod$pval))]),
      k              = mod$k
    )
  }) %>%
    ungroup()
}

plot_global_comparison <- function(estimates, outfile){
  if (!nrow(estimates)) return(invisible(NULL))
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

  df_plot <- estimates %>%
    mutate(
      session_label = sprintf("β_session = %.3f (p=%s)", beta_session, fmt_p(p_session)),
      diff_label    = sprintf("Δβ = %.3f (p=%s)", beta_interact, fmt_p(p_interaction)),
      stars         = sig_stars_vec(p_interaction)
    )

  p <- ggplot(df_plot, aes(x = reorder(molecule, beta_session), y = beta_session)) +
    geom_col(fill = "#4E79A7") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_text(aes(label = stars), vjust = -0.6, color = "red3", size = 4) +
    coord_flip() +
    labs(
      title = "Dose–response slope during session (global AE)",
      x = "Molecule",
      y = "β (log-OR change per dose unit)"
    ) +
    theme_minimal(base_size = 13)

  ggsave(outfile, p, width = 7.5, height = 5.5, dpi = 300)
  invisible(p)
}

run_compare_global_session_followup <- function(
    data_xlsx = here::here("data", "Adverse-events-dose-v5.xlsx"),
    sheet = "Feuil1",
    out_dir = here::here("results", "compare", "global")
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  es <- build_es_all(data_xlsx, sheet, ref_policies = default_ref_policies)
  est <- compare_global_dose_slope_by_window(es, min_k = 4)
  write_csv(est, file.path(out_dir, "global_session_followup_dose_slopes.csv"))
  plot_global_comparison(est, file.path(out_dir, "global_session_vs_followup_slopes.pdf"))
  invisible(est)
}

if (identical(environment(), globalenv())) {
  run_compare_global_session_followup()
}
