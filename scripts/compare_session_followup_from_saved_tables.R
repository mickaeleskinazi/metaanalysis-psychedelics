suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

# ================== CONFIG ==================
DIR_SESSION  <- "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/results_session/tables"
DIR_FOLLOWUP <- "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/results_followup/tables"
OUT_DIR      <- "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/results_compare/tables"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# helper étoiles
sig_stars <- function(p){
  case_when(
    is.na(p)        ~ "",
    p < 0.001       ~ "***",
    p < 0.01        ~ "**",
    p < 0.05        ~ "*",
    TRUE            ~ ""
  )
}

# ---------- 1) LECTURE SÛRE DES FICHIERS ----------
.read_csv_safe <- function(path){
  if (!file.exists(path)) {
    warning("Missing file: ", path)
    return(NULL)
  }
  tryCatch(readr::read_csv(path, show_col_types = FALSE),
           error = function(e) { warning("Read error on ", path, ": ", e$message); NULL })
}

# fichiers attendus
f_s_byMol   <- file.path(DIR_SESSION,  "significance_by_molecule_models.csv")
f_f_byMol   <- file.path(DIR_FOLLOWUP, "significance_by_molecule_models.csv")

f_s_byAE    <- file.path(DIR_SESSION,  "significance_by_ae_models.csv")
f_f_byAE    <- file.path(DIR_FOLLOWUP, "significance_by_ae_models.csv")

f_s_aggMol  <- file.path(DIR_SESSION,  "significance_agg_by_molecule.csv")
f_f_aggMol  <- file.path(DIR_FOLLOWUP, "significance_agg_by_molecule.csv")

f_s_aggAEM  <- file.path(DIR_SESSION,  "significance_agg_by_ae_molecule.csv")
f_f_aggAEM  <- file.path(DIR_FOLLOWUP, "significance_agg_by_ae_molecule.csv")

byMol_s  <- .read_csv_safe(f_s_byMol)
byMol_f  <- .read_csv_safe(f_f_byMol)
byAE_s   <- .read_csv_safe(f_s_byAE)
byAE_f   <- .read_csv_safe(f_f_byAE)
aggMol_s <- .read_csv_safe(f_s_aggMol)
aggMol_f <- .read_csv_safe(f_f_aggMol)
aggAEM_s <- .read_csv_safe(f_s_aggAEM)
aggAEM_f <- .read_csv_safe(f_f_aggAEM)

# ---------- 2) HARMONISATION MINIMALE ----------
# On se concentre sur le modèle spline_df3 si présent (ou on garde tout sinon)
.keep_or_all <- function(df){
  if (is.null(df)) return(NULL)
  if ("model" %in% names(df) && any(df$model == "spline_df3")) {
    df %>% filter(model == "spline_df3")
  } else df
}

byMol_s  <- .keep_or_all(byMol_s)
byMol_f  <- .keep_or_all(byMol_f)
byAE_s   <- .keep_or_all(byAE_s)
byAE_f   <- .keep_or_all(byAE_f)
aggMol_s <- .keep_or_all(aggMol_s)
aggMol_f <- .keep_or_all(aggMol_f)
aggAEM_s <- .keep_or_all(aggAEM_s)
aggAEM_f <- .keep_or_all(aggAEM_f)

# ---------- 3) COMPARAISON PAR MOLÉCULE (modèles détaillés) ----------
# On agrège par molécule les p-values pertinentes (hors intercept)
collapse_byMol <- function(df){
  if (is.null(df)) return(NULL)
  # colonnes possibles dans tes exports
  # molecule, term, pval, k, I2, tau2, QM, QMp, estimate/ci_low/ci_high, etc.
  df %>%
    mutate(term = ifelse(is.na(term), "", term)) %>%
    filter(!grepl("intrcpt", term, ignore.case = TRUE)) %>%
    group_by(molecule) %>%
    summarise(
      k_total = suppressWarnings(max(k, na.rm = TRUE)),
      I2      = suppressWarnings(max(I2, na.rm = TRUE)),
      tau2    = suppressWarnings(max(tau2, na.rm = TRUE)),
      p_overall = suppressWarnings(
        if ("QMp" %in% names(.)) min(QMp, na.rm = TRUE) else min(pval, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    mutate(stars = sig_stars(p_overall))
}

byMol_s_coll <- collapse_byMol(byMol_s) %>% mutate(window = "session")
byMol_f_coll <- collapse_byMol(byMol_f) %>% mutate(window = "follow_up")

cmp_byMol <- full_join(
  byMol_s_coll %>% rename(k_session = k_total, I2_session = I2, tau2_session = tau2,
                          p_overall_session = p_overall, stars_session = stars),
  byMol_f_coll %>% rename(k_follow = k_total, I2_follow = I2, tau2_follow = tau2,
                          p_overall_follow = p_overall, stars_follow = stars),
  by = "molecule"
)

write_csv(cmp_byMol, file.path(OUT_DIR, "compare_by_molecule_overall.csv"))

# ---------- 4) COMPARAISON PAR AE × MOLÉCULE (modèles détaillés) ----------
collapse_byAE <- function(df){
  if (is.null(df)) return(NULL)
  df %>%
    mutate(term = ifelse(is.na(term), "", term)) %>%
    filter(!grepl("intrcpt", term, ignore.case = TRUE)) %>%
    group_by(ae_term, molecule) %>%
    summarise(
      k_total = suppressWarnings(max(k, na.rm = TRUE)),
      p_overall = suppressWarnings(
        if ("QMp" %in% names(.)) min(QMp, na.rm = TRUE) else min(pval, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    mutate(stars = sig_stars(p_overall))
}

byAE_s_coll <- collapse_byAE(byAE_s) %>% mutate(window = "session")
byAE_f_coll <- collapse_byAE(byAE_f) %>% mutate(window = "follow_up")

cmp_byAE <- full_join(
  byAE_s_coll %>% rename(k_session = k_total, p_session = p_overall, stars_session = stars),
  byAE_f_coll %>% rename(k_follow  = k_total, p_follow  = p_overall, stars_follow  = stars),
  by = c("ae_term","molecule")
)

write_csv(cmp_byAE, file.path(OUT_DIR, "compare_by_ae_molecule.csv"))

# ---------- 5) COMPARAISON AGGRÉGÉE (déjà agrégée dans tes exports) ----------
# (même idée, mais on fait juste un side-by-side rapide)
cmp_agg_byMol <- full_join(
  (aggMol_s %>% rename(k_session = k_total, QM_session = QM, p_session = QMp, stars_session = stars)),
  (aggMol_f %>% rename(k_follow  = k_total, QM_follow  = QM, p_follow  = QMp, stars_follow  = stars)),
  by = "molecule"
)
write_csv(cmp_agg_byMol, file.path(OUT_DIR, "compare_agg_by_molecule.csv"))

cmp_agg_byAEM <- full_join(
  (aggAEM_s %>% rename(k_session = k_total, QM_session = QM, p_session = QMp, stars_session = stars)),
  (aggAEM_f %>% rename(k_follow  = k_total, QM_follow  = QM, p_follow  = QMp, stars_follow  = stars)),
  by = c("ae_term","molecule")
)
write_csv(cmp_agg_byAEM, file.path(OUT_DIR, "compare_agg_by_ae_molecule.csv"))

# ---------- 6) TABLEAU "TOPLINE" POUR LE PAPIER ----------
# Pour chaque molécule : p_overall session, p_overall follow_up, nb d'AE significatifs dans chaque fenêtre
count_sig_ae <- function(df_coll){
  if (is.null(df_coll)) return(tibble(molecule=character(), n_sig=integer()))
  df_coll %>%
    filter(!is.na(p_overall)) %>%
    mutate(sig = p_overall < 0.05) %>%
    group_by(molecule) %>%
    summarise(n_sig = sum(sig, na.rm = TRUE), .groups = "drop")
}

n_sig_s <- count_sig_ae(byAE_s_coll) %>% rename(n_sig_session = n_sig)
n_sig_f <- count_sig_ae(byAE_f_coll) %>% rename(n_sig_follow  = n_sig)

topline <- cmp_byMol %>%
  select(molecule, p_overall_session, stars_session, p_overall_follow, stars_follow) %>%
  left_join(n_sig_s, by = "molecule") %>%
  left_join(n_sig_f, by = "molecule") %>%
  mutate(
    session_str  = ifelse(is.na(p_overall_session), NA,
                          sprintf("p=%.3g %s", p_overall_session, ifelse(is.na(stars_session),"",stars_session))),
    followup_str = ifelse(is.na(p_overall_follow),  NA,
                          sprintf("p=%.3g %s", p_overall_follow,  ifelse(is.na(stars_follow),"",stars_follow)))
  ) %>%
  transmute(
    Molecule = molecule,
    `Overall (session)`  = session_str,
    `Overall (follow-up)`= followup_str,
    `# significant AE (session)`  = n_sig_session,
    `# significant AE (follow-up)`= n_sig_follow
  )

write_csv(topline, file.path(OUT_DIR, "compare_topline_publication.csv"))

# ---------- 7) EXPORT LATEX (tableaux compacts) ----------
to_latex <- function(df, caption, label){
  if (is.null(df) || !nrow(df)) return("")
  cols <- names(df)
  align <- paste(rep("l", length(cols)), collapse="")
  header <- paste(cols, collapse=" & ")
  body <- apply(df, 1, function(row) paste(ifelse(is.na(row), "", row), collapse=" & "))
  body <- paste(body, collapse=" \\\\\n")
  tex <- sprintf("\\begin{table}[ht]\n\\centering\n\\caption{%s}\n\\label{%s}\n\\begin{tabular}{%s}\n\\toprule\n%s \\\\\n\\midrule\n%s \\\\\n\\bottomrule\n\\end{tabular}\n\\end{table}\n",
                 caption, label, align, header, body)
  tex
}

tex1 <- to_latex(topline,
                 "Overall significance by molecule (session vs follow-up) and number of significant AE.",
                 "tab:cmp_topline")
writeLines(tex1, file.path(OUT_DIR, "compare_topline_publication.tex"))

tex2 <- to_latex(
  cmp_byAE %>%
    arrange(molecule, ae_term) %>%
    mutate(
      session = ifelse(is.na(p_session), "", sprintf("p=%.3g %s", p_session, ifelse(is.na(stars_session),"",stars_session))),
      follow  = ifelse(is.na(p_follow),  "", sprintf("p=%.3g %s", p_follow,  ifelse(is.na(stars_follow),"",stars_follow)))
    ) %>% select(molecule, ae_term, session, follow),
  "Per-AE comparison by molecule: session vs follow-up p-values.",
  "tab:cmp_ae_molecule"
)
writeLines(tex2, file.path(OUT_DIR, "compare_by_ae_molecule.tex"))

message("✅ Done. CSV + LaTeX saved in: ", OUT_DIR)