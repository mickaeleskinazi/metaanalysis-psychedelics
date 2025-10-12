suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(metafor)
})

# ==== CONFIG ====
BASE_XLSX <- "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/data/Adverse-events-dose-v5.xlsx"
SHEET     <- "Feuil1"
OUT_DIR   <- "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/results_compare"

# Project functions
source("/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/scripts/utils_data.R")
source("/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/scripts/analysis_dose_response.R")

# Reference priority (same as your runs)
default_ref_policies <- list(
  MDMA       = c("inactive_placebo","active_non_psy_placebo","active_placebo"),
  LSD        = c("inactive_placebo","active_placebo","active_non_psy_placebo"),
  PSILOCYBIN = c("inactive_placebo","active_non_psy_placebo","active_placebo"),
  AYAHUASCA  = c("inactive_placebo","active_placebo","active_non_psy_placebo"),
  .default   = c("inactive_placebo","active_placebo","active_non_psy_placebo")
)

# ---- helpers ----
norm_win <- function(x){
  x |>
    as.character() |>
    tolower() |>
    gsub("[^a-z0-9]+","_", x = _) |>
    gsub("_+","_", x = _) |>
    sub("^_","", .) |>
    sub("_$","", .)
}

build_es_for_window <- function(window_value){
  raw_all <- load_data(BASE_XLSX, sheet = SHEET) %>%
    mutate(time_window = norm_win(time_window))
  if (!window_value %in% unique(raw_all$time_window)) {
    stop("Window '", window_value, "' not found. Available: ",
         paste(sort(unique(raw_all$time_window)), collapse=", "))
  }
  raw_win <- raw_all %>% filter(time_window == window_value)
  contr   <- build_pairwise_2x2(raw_win, ref_policies = default_ref_policies)
  es      <- build_escalc(contr)
  es$time_window <- window_value           # << keep the window label
  es
}

# ==== BUILD BOTH WINDOWS ====
message("== Building ES for session and follow_up …")
es_session  <- try(build_es_for_window("session"))
es_followup <- try(build_es_for_window("follow_up"))

es_all <- dplyr::bind_rows(
  if (!inherits(es_session,  "try-error")) es_session,
  if (!inherits(es_followup,"try-error")) es_followup
)

if (!nrow(es_all)) stop("No ES rows after building session/follow_up.")

message("Counts by window:")
print(es_all %>% count(time_window))

# ==== TABLE MAKERS (inline; same logic I gave you earlier) ====
.sig_stars <- function(p){
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

# Per-window slope within each molecule
dr_fit_per_window <- function(es, min_k_per_window = 2){
  stopifnot(all(c("molecule","dose_mg","time_window","yi","vi") %in% names(es)))
  es %>%
    filter(!is.na(dose_mg), is.finite(yi), is.finite(vi)) %>%
    group_by(molecule, time_window) %>%
    group_modify(~{
      dat <- .x
      if (nrow(dat) < min_k_per_window) return(NULL)
      m <- tryCatch(rma(yi ~ dose_mg, vi = vi, data = dat, method = "REML"), error = function(e) NULL)
      if (is.null(m)) return(NULL)
      tibble(
        molecule       = dat$molecule[[1]],
        time_window    = dat$time_window[[1]],
        k              = m$k,
        beta_dose      = as.numeric(coef(m)["dose_mg"]),
        se_dose        = as.numeric(m$se["dose_mg"]),
        z_dose         = as.numeric(m$zval["dose_mg"]),
        p_dose         = as.numeric(m$pval["dose_mg"]),
        ci_lb          = as.numeric(m$ci.lb[["dose_mg"]]),
        ci_ub          = as.numeric(m$ci.ub[["dose_mg"]]),
        I2             = suppressWarnings(tryCatch(m$I2,  error=function(e) NA_real_)),
        tau2           = suppressWarnings(tryCatch(m$tau2, error=function(e) NA_real_))
      )
    }) %>% ungroup() %>%
    mutate(stars = .sig_stars(p_dose))
}

# Interaction test within each molecule: yi ~ dose_mg * time_window
dr_test_session_vs_followup <- function(es, min_k_total = 4){
  stopifnot(all(c("molecule","dose_mg","time_window","yi","vi") %in% names(es)))
  es %>%
    filter(!is.na(dose_mg), is.finite(yi), is.finite(vi)) %>%
    group_by(molecule) %>%
    group_modify(~{
      dat <- .x
      if (length(unique(dat$time_window)) < 2 || nrow(dat) < min_k_total) return(NULL)
      m <- tryCatch(rma(yi ~ dose_mg * time_window, vi = vi, data = dat, method = "REML"), error=function(e) NULL)
      if (is.null(m)) return(NULL)
      tibble(
        molecule          = dat$molecule[[1]],
        beta_dose_main    = as.numeric(coef(m)["dose_mg"]),
        beta_interaction  = as.numeric(coef(m)[grep("^dose_mg:time_window", names(coef(m)))]),
        p_interaction     = as.numeric(m$pval[grep("^dose_mg:time_window", names(m$pval))])
      )
    }) %>% ungroup() %>%
    mutate(stars_interaction = .sig_stars(p_interaction))
}

# Wrapper to export all tables
make_dr_window_tables <- function(
    es,
    out_dir = OUT_DIR,
    min_k_per_window = 2,
    min_k_total      = 4
){
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  per_window <- dr_fit_per_window(es, min_k_per_window = min_k_per_window)
  write_csv(per_window, file.path(out_dir, "dr_per_window_by_molecule.csv"))
  
  per_window_wide <- per_window %>%
    select(molecule, time_window, beta_dose, ci_lb, ci_ub, p_dose, stars, k, I2, tau2) %>%
    mutate(time_window = case_when(
      time_window %in% c("session","Session") ~ "session",
      time_window %in% c("follow_up","follow-up","Follow_up") ~ "follow_up",
      TRUE ~ time_window
    )) %>%
    pivot_wider(
      id_cols = molecule,
      names_from = time_window,
      values_from = c(beta_dose, ci_lb, ci_ub, p_dose, stars, k, I2, tau2),
      names_sep = "."
    )
  write_csv(per_window_wide, file.path(out_dir, "dr_per_window_by_molecule_wide.csv"))
  
  interact <- dr_test_session_vs_followup(es, min_k_total = min_k_total)
  write_csv(interact, file.path(out_dir, "dr_session_vs_followup_interaction.csv"))
  
  pub_tab <- per_window_wide %>%
    left_join(interact %>% select(molecule, beta_interaction, p_interaction, stars_interaction), by="molecule") %>%
    mutate(
      session_effect   = ifelse(!is.na(beta_dose.session),
                                sprintf("%.3f [%.3f, %.3f]%s",
                                        beta_dose.session, ci_lb.session, ci_ub.session,
                                        ifelse(is.na(stars.session), "", paste0(" ", stars.session))), NA),
      followup_effect  = ifelse(!is.na(beta_dose.follow_up),
                                sprintf("%.3f [%.3f, %.3f]%s",
                                        beta_dose.follow_up, ci_lb.follow_up, ci_ub.follow_up,
                                        ifelse(is.na(stars.follow_up), "", paste0(" ", stars.follow_up))), NA),
      interaction_str  = ifelse(!is.na(p_interaction),
                                sprintf("β×window=%.3f; p=%.3g%s",
                                        beta_interaction, p_interaction,
                                        ifelse(is.na(stars_interaction), "", paste0(" ", stars_interaction))), NA)
    ) %>%
    transmute(
      Molecule = molecule,
      `Session slope (β, 95% CI)`   = session_effect,
      `Follow-up slope (β, 95% CI)` = followup_effect,
      `Window interaction`          = interaction_str
    )
  write_csv(pub_tab, file.path(out_dir, "dr_session_followup_publication_table.csv"))
  
  message("Saved tables in: ", out_dir)
  invisible(list(per_window = per_window,
                 per_window_wide = per_window_wide,
                 interaction = interact,
                 publication_table = pub_tab))
}

# ==== RUN ====
out <- make_dr_window_tables(es_all, out_dir = OUT_DIR, min_k_per_window = 2, min_k_total = 4)