suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(janitor)
  library(readr)
})

# === Chemins ===
INFILE <- "/Users/mickaeleskinazi/Downloads/forest_by_molecule_ae_window_publication.xlsx"
OUTDIR <- "/Users/mickaeleskinazi/Downloads/tables_pub"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# === Utilitaires ===
.norm_window <- function(x){
  x0 <- tolower(trimws(gsub("[[:space:]]+", " ", x)))
  dplyr::case_when(
    grepl("session|acute|during", x0) ~ "session",
    grepl("follow[_ \\-]?up|post|fu|day|week|month", x0) ~ "follow_up",
    TRUE ~ x0
  )
}

# robust match of possibly accented/stat symbols
.find_col <- function(nms, candidates){
  ix <- which(tolower(nms) %in% tolower(candidates))
  if (length(ix)) return(nms[ix[1]])
  # try regex/contains
  ix2 <- which(Reduce(`|`, lapply(candidates, function(c) grepl(c, nms, ignore.case = TRUE))))
  if (length(ix2)) return(nms[ix2[1]])
  return(NA_character_)
}

# === Lecture & nettoyage ===
stopifnot(file.exists(INFILE))
raw <- readxl::read_excel(INFILE)
df <- janitor::clean_names(raw)

# Reconnaître les colonnes clés même si renommées
cn <- names(df)

col_mol   <- .find_col(cn, c("molecule","drug"))
col_ae    <- .find_col(cn, c("ae_term","ae","adverse_event","event"))
col_win   <- .find_col(cn, c("window","time_window","tw"))
col_k     <- .find_col(cn, c("k"))
col_or    <- .find_col(cn, c("^or$","odds_ratio"))
col_logor <- .find_col(cn, c("log\\(?or\\)?","log_or","yi"))
col_p     <- .find_col(cn, c("^p$","pval","p_value"))
col_stars <- .find_col(cn, c("^stars$","signif","sig"))

req <- c(col_mol, col_ae, col_win, col_k, col_or, col_logor, col_p)
if (anyNA(req)) {
  stop("Colonnes manquantes. Trouvées: ",
       "\n- molecule: ", col_mol,
       "\n- ae_term:  ", col_ae,
       "\n- window:   ", col_win,
       "\n- k:        ", col_k,
       "\n- OR:       ", col_or,
       "\n- log(OR):  ", col_logor,
       "\n- p:        ", col_p,
       "\nVérifie les en-têtes de ", INFILE)
}

df1 <- df %>%
  rename(
    molecule = !!col_mol,
    ae_term  = !!col_ae,
    window   = !!col_win,
    k        = !!col_k,
    OR       = !!col_or,
    logOR    = !!col_logor,
    p        = !!col_p,
    Stars    = dplyr::all_of(col_stars)
  ) %>%
  mutate(
    window = .norm_window(window),
    sig    = is.finite(p) & p < 0.05
  ) %>%
  filter(window %in% c("session","follow_up")) %>%
  mutate(
    molecule = as.character(molecule),
    ae_term  = as.character(ae_term)
  )

# === Table 1 : AE significatifs par molécule et par fenêtre (liste compacte) ===
ae_sig_long <- df1 %>%
  filter(sig) %>%
  group_by(molecule, window) %>%
  summarise(
    n_sig = n(),
    ae_list = paste(sort(unique(ae_term)), collapse = ", "),
    .groups = "drop"
  )

ae_sig_wide <- ae_sig_long %>%
  mutate(window = if_else(window=="session","session","follow_up")) %>%
  select(molecule, window, n_sig, ae_list) %>%
  pivot_wider(names_from = window,
              values_from = c(n_sig, ae_list),
              values_fill = list(n_sig = 0, ae_list = "")) %>%
  arrange(molecule)

# Sauvegardes
write_csv(ae_sig_wide, file.path(OUTDIR, "forest_sig_ae_by_molecule_window.csv"))

# LaTeX (liste)
tex1 <- ae_sig_wide %>%
  mutate(
    n_sig_session   = ifelse(is.na(n_sig_session), 0, n_sig_session),
    n_sig_follow_up = ifelse(is.na(n_sig_follow_up), 0, n_sig_follow_up)
  ) %>%
  transmute(
    Molecule = molecule,
    `# AE sig. (session)`   = n_sig_session,
    `AE (session)`          = ifelse(ae_list_session=="","—",ae_list_session),
    `# AE sig. (follow-up)` = n_sig_follow_up,
    `AE (follow-up)`        = ifelse(ae_list_follow_up=="","—",ae_list_follow_up)
  )

cat(
  "\\begin{table}[ht]\n\\centering\n\\small\n",
  "\\caption{Significant adverse events (forest meta-analysis, $p<0.05$) by molecule and time window.}\n",
  "\\label{tab:forest_ae_sig_by_window}\n",
  "\\begin{tabular}{l r L{5.8cm} r L{5.8cm}}\n\\toprule\n",
  "\\textbf{Molecule} & \\textbf{\\# sig. (Sess.)} & \\textbf{AE (Sess.)} & \\textbf{\\# sig. (FU)} & \\textbf{AE (FU)}\\\\\n\\midrule\n",
  sep = "", file = file.path(OUTDIR, "forest_ae_sig_by_molecule_window.tex")
)

apply(tex1, 1, function(r){
  line <- sprintf("%s & %s & %s & %s & %s \\\\",
                  r[["Molecule"]], r[[2]], r[[3]], r[[4]], r[[5]])
  cat(line, "\n", file = file.path(OUTDIR, "forest_ae_sig_by_molecule_window.tex"), append = TRUE)
})

cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n",
    file = file.path(OUTDIR, "forest_ae_sig_by_molecule_window.tex"),
    append = TRUE)

# === Table 2 : Classification transient / emergent / persistent ===
# pivot pour repérer la signif. par fenêtre
sig_mat <- df1 %>%
  select(molecule, ae_term, window, sig) %>%
  distinct() %>%
  pivot_wider(names_from = window, values_from = sig, values_fill = FALSE)

# statut par AE×molécule
status_ae <- sig_mat %>%
  mutate(
    status = case_when(
      session & follow_up ~ "Persistent",
      session & !follow_up ~ "Transient (session-only)",
      !session & follow_up ~ "Emergent (follow-up-only)",
      TRUE ~ "Non-significant"
    )
  )

# résumé par molécule
status_summary <- status_ae %>%
  count(molecule, status) %>%
  pivot_wider(names_from = status, values_from = n, values_fill = 0) %>%
  mutate(
    `Total AE tested` = `Persistent` + `Transient (session-only)` + `Emergent (follow-up-only)` + `Non-significant`
  ) %>%
  arrange(molecule)

write_csv(status_summary, file.path(OUTDIR, "forest_ae_persistence_summary.csv"))

# LaTeX (résumé persistance)
cols_keep <- c("molecule","Transient (session-only)","Emergent (follow-up-only)","Persistent","Total AE tested")
ss <- status_summary %>% select(all_of(cols_keep))

cat(
  "\\begin{table}[ht]\n\\centering\n\\small\n",
  "\\caption{Temporal classification of AE significance per molecule (forest meta-analysis).}\n",
  "\\label{tab:forest_persistence_summary}\n",
  "\\begin{tabular}{l r r r r}\n\\toprule\n",
  "\\textbf{Molecule} & \\textbf{Transient} & \\textbf{Emergent} & \\textbf{Persistent} & \\textbf{Total}\\\\\n\\midrule\n",
  sep = "", file = file.path(OUTDIR, "forest_persistence_summary.tex")
)

apply(ss, 1, function(r){
  line <- sprintf("%s & %d & %d & %d & %d \\\\",
                  r[["molecule"]], r[[2]], r[[3]], r[[4]], r[[5]])
  cat(line, "\n", file = file.path(OUTDIR, "forest_persistence_summary.tex"), append = TRUE)
})

cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n",
    file = file.path(OUTDIR, "forest_persistence_summary.tex"),
    append = TRUE)

# === (Optionnel) Export long propre pour annexes ===
df_export <- df1 %>%
  mutate(
    Stars = ifelse(is.na(Stars) & sig, "*", Stars),
    Stars = ifelse(is.na(Stars), "", Stars)
  ) %>%
  arrange(molecule, ae_term, window)

write_csv(df_export, file.path(OUTDIR, "forest_by_molecule_ae_window_clean.csv"))

message("✅ Tables générées dans: ", OUTDIR)