suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(metafor)
  library(ggplot2)
})

# === CONFIG ===
BASE_XLSX <- "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/data/Adverse-events-dose-v5.xlsx"
SHEET     <- "Feuil1"

UTILS     <- "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/scripts/utils_data.R"
AN_DR     <- "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/scripts/analysis_dose_response.R"

OUT_DIR   <- "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/results_compare"
OUT_TBL   <- file.path(OUT_DIR, "tables")
OUT_FIG   <- file.path(OUT_DIR, "figures")

dir.create(OUT_TBL, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# Reference/contrast policy (same as your main runs)
default_ref_policies <- list(
  MDMA       = c("inactive_placebo","active_non_psy_placebo","active_placebo"),
  LSD        = c("inactive_placebo","active_placebo","active_non_psy_placebo"),
  PSILOCYBIN = c("inactive_placebo","active_non_psy_placebo","active_placebo"),
  AYAHUASCA  = c("inactive_placebo","active_placebo","active_non_psy_placebo"),
  .default   = c("inactive_placebo","active_placebo","active_non_psy_placebo")
)

# === SOURCE PROJECT HELPERS ===
source(UTILS)   # must define: load_data(), build_pairwise_2x2(), build_escalc()
source(AN_DR)   # ok if not used, but harmless

# --- Helper: normalise window labels to {session, follow_up} with "session" as baseline
.norm_window <- function(x){
  x <- tolower(trimws(as.character(x)))
  x <- ifelse(grepl("follow", x), "follow_up", x)
  x <- ifelse(grepl("^fu$", x), "follow_up", x)
  x <- ifelse(grepl("session|acute|day 0|seance|séance", x), "session", x)
  factor(x, levels = c("session","follow_up"))
}

# --- Safe extract from named vectors (coeffs/pvals)
.safe <- function(x, nm) {
  if (is.null(x)) return(NA_real_)
  nms <- names(x)
  if (is.null(nms) || !(nm %in% nms)) return(NA_real_)
  as.numeric(x[[nm]])
}

# --- Build one combined ES with both windows (all AE pooled)
build_es_all <- function(base_xlsx, sheet, ref_policies){
  suppressPackageStartupMessages({ library(readxl); library(dplyr); library(stringr); library(forcats) })
  
  # ---- 1) Read your v5 sheet exactly ----
  df <- readxl::read_excel(base_xlsx, sheet = sheet)
  
  # Sanity: required raw columns (your exact names)
  req <- c("study_id","author_year","arm_id","molecule",
           "n_participants_arm","ae_term","time_window","events","dose_mg")
  miss <- setdiff(req, names(df))
  if (length(miss)) stop("Missing expected Excel columns: ", paste(miss, collapse=", "))
  
  # ---- 2) Harmonize to what downstream needs ----
  raw <- df %>%
    transmute(
      study_id,
      author_year,
      molecule         = as.character(molecule),
      ae_term          = as.character(ae_term),
      dose_mg          = suppressWarnings(as.numeric(dose_mg)),
      # copy arm_id into group
      group            = as.character(arm_id),
      # arm_type from arm_id
      arm_type         = case_when(
        str_detect(tolower(arm_id), "placebo") ~ "placebo",
        TRUE                                   ~ "active"
      ),
      # counts
      n_total          = suppressWarnings(as.numeric(n_participants_arm)),
      n_events         = suppressWarnings(as.numeric(events)),
      # normalize window → {session, follow_up}
      time_window      = {
        tw <- tolower(trimws(as.character(time_window)))
        tw <- ifelse(str_detect(tw, "follow|follow[-_ ]?up|^fu$"), "follow_up", tw)
        tw <- ifelse(str_detect(tw, "session|acute|day 0|seance|séance"), "session", tw)
        forcats::fct_drop(factor(tw, levels = c("session","follow_up")))
      }
    )
  
  # Final check for downstream fields
  need <- c("study_id","molecule","ae_term","time_window","group","arm_type",
            "dose_mg","n_events","n_total")
  miss2 <- setdiff(need, names(raw))
  if (length(miss2)) stop("Harmonized data missing: ", paste(miss2, collapse=", "))
  
  # Basic cleaning
  raw <- raw %>%
    filter(!is.na(molecule), !is.na(ae_term), !is.na(dose_mg),
           !is.na(n_total), !is.na(n_events), !is.na(time_window))
  
  # ---- 3) Build pairwise contrasts + log-ORs using your existing helpers ----
  contr <- build_pairwise_2x2(raw, ref_policies = ref_policies)
  es    <- build_escalc(contr)
  
  # Keep usable rows
  es %>%
    filter(!is.na(time_window), is.finite(yi), is.finite(vi), !is.na(dose_mg))
}

# --- Core: compare global (all-AE pooled) dose slopes between session vs follow_up
compare_global_dose_slope_by_window <- function(es, min_k = 4) {
  stopifnot(all(c("molecule","dose_mg","time_window","yi","vi") %in% names(es)))
  es %>%
    group_by(molecule) %>%
    group_modify(~{
      dat <- .x
      # need both windows and enough rows
      if (length(unique(dat$time_window)) < 2 || nrow(dat) < min_k) {
        return(tibble(
          molecule = dat$molecule[[1]],
          beta_dose_session  = NA_real_,
          beta_dose_followup = NA_real_,
          beta_diff          = NA_real_,
          p_interaction      = NA_real_,
          stars              = ""
        ))
      }
      
      m <- tryCatch(
        metafor::rma(yi ~ dose_mg * time_window, vi = vi, data = dat, method = "REML"),
        error = function(e) NULL
      )
      if (is.null(m)) {
        return(tibble(
          molecule = dat$molecule[[1]],
          beta_dose_session  = NA_real_,
          beta_dose_followup = NA_real_,
          beta_diff          = NA_real_,
          p_interaction      = NA_real_,
          stars              = ""
        ))
      }
      
      b_dose <- .safe(m$b, "dose_mg")
      b_diff <- .safe(m$b, "dose_mg:time_windowfollow_up")
      p_int  <- .safe(m$pval, "dose_mg:time_windowfollow_up")
      
      tibble(
        molecule = dat$molecule[[1]],
        beta_dose_session  = b_dose,
        beta_dose_followup = ifelse(is.na(b_dose) | is.na(b_diff), NA_real_, b_dose + b_diff),
        beta_diff          = b_diff,
        p_interaction      = p_int,
        stars = case_when(
          is.na(p_int) ~ "",
          p_int < 0.001 ~ "***",
          p_int < 0.01  ~ "**",
          p_int < 0.05  ~ "*",
          TRUE ~ ""
        )
      )
    }) %>% ungroup()
}

# --- Optional: quick LaTeX builder (simple, journal-neutral)
to_latex <- function(df, caption, label){
  if (!nrow(df)) return("% empty table")
  cols <- names(df)
  align <- paste(rep("l", length(cols)), collapse = "")
  header <- paste(cols, collapse = " & ")
  body <- apply(df, 1, function(x) paste(x, collapse = " & ")) |> paste(collapse = " \\\\\n")
  sprintf("\\begin{table}[ht]\n\\centering\n\\caption{%s}\n\\label{%s}\n\\begin{tabular}{%s}\n\\toprule\n%s \\\\\n\\midrule\n%s \\\\\n\\bottomrule\n\\end{tabular}\n\\end{table}",
          caption, label, align, header, body)
}

# --- Optional: compact plot showing both slopes per molecule
plot_global_slopes <- function(tbl, outfile){
  df <- tbl %>%
    tidyr::pivot_longer(c(beta_dose_session, beta_dose_followup),
                        names_to = "window", values_to = "beta") %>%
    mutate(window = recode(window,
                           beta_dose_session = "session",
                           beta_dose_followup = "follow_up"))
  p <- ggplot(df, aes(x = molecule, y = beta, group = molecule)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(aes(color = molecule), linewidth = 0.6, show.legend = FALSE) +
    geom_point(aes(shape = window, color = molecule), size = 3) +
    scale_shape_manual(values = c(session = 16, follow_up = 17)) +
    labs(title = "Global dose–response slopes by molecule: session vs follow-up",
         y = "Slope (log-OR per mg)", x = NULL) +
    theme_bw() +
    theme(legend.position = "bottom")
  ggsave(outfile, p, width = 8, height = 5, dpi = 150)
}

# === RUN ===
message("== Building combined ES with session + follow_up …")
es_all <- build_es_all(BASE_XLSX, SHEET, default_ref_policies)

message("== Comparing global slopes by molecule …")
tbl_global_cmp <- compare_global_dose_slope_by_window(es_all, min_k = 4)

# Save CSV
csv_out <- file.path(OUT_TBL, "compare_global_dose_slope_by_molecule.csv")
readr::write_csv(tbl_global_cmp, csv_out)

# Save LaTeX
tex_out <- file.path(OUT_TBL, "compare_global_dose_slope_by_molecule.tex")
writeLines(
  to_latex(tbl_global_cmp,
           "Comparison of global dose--response slopes between session and follow-up for each molecule.",
           "tab:global_dose_slope"),
  tex_out
)

# Plot
fig_out <- file.path(OUT_FIG, "compare_global_dose_slope_by_molecule.png")
plot_global_slopes(tbl_global_cmp, fig_out)

message("✅ Done.")
message("Tables: ", csv_out, " and ", tex_out)
message("Figure: ", fig_out)