suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(janitor)
  library(tidyr)
  library(purrr)
  library(glue)
  library(readr)
})

# -------------------- Chemins --------------------
INPUT_XLS <- "/Users/mickaeleskinazi/Downloads/forest_by_molecule_ae_window_publication.xlsx"
OUT_DIR   <- "/Users/mickaeleskinazi/Downloads/tables_pub"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -------------------- Helpers --------------------
.norm_window <- function(x){
  x %>%
    as.character() %>%
    stringi::stri_trans_general("Any-Latin; Latin-ASCII") %>%
    str_squish() %>%
    tolower() %>%
    str_replace_all("[-_/]", " ") %>%
    { case_when(
      str_detect(., "\\bsession\\b|acute|during|day 0|seance|séance") ~ "session",
      str_detect(., "follow ?up|post|\\bfu\\b|suivi|week|month|jours") ~ "follow_up",
      TRUE ~ .
    )}
}

# Cherche une colonne parmi des candidats (noms déjà "clean_names")
.find_col <- function(nm, candidates, desc="column"){
  ix <- which(tolower(nm) %in% tolower(candidates))
  if (length(ix)) return(nm[ix[1]])
  stop(glue("Colonne '{desc}' introuvable. Colonnes disponibles : {paste(nm, collapse=', ')}"))
}

# Coalesce plusieurs colonnes potentielles vers une seule, puis supprime les originales
coalesce_to <- function(data, candidates, out_name){
  cols <- intersect(names(data), candidates)
  if (length(cols) == 0) return(data)
  if (length(cols) == 1) {
    data[[out_name]] <- data[[cols]]
  } else {
    data[[out_name]] <- dplyr::coalesce(!!!rlang::syms(cols))
  }
  data[, setdiff(names(data), cols), drop = FALSE]
}

# -------------------- Lecture & nettoyage --------------------
raw <- readxl::read_excel(INPUT_XLS, col_names = TRUE)

# Si aucune colonne "molecule"/"ae_term"/"window" détectable, tenter : 1ère ligne = en-têtes
needs_header_fix <- function(df){
  nm <- janitor::make_clean_names(names(df))
  !any(str_detect(nm, "molecule|drug")) ||
    !any(str_detect(nm, "ae_term|adverse|event|ae")) ||
    !any(str_detect(nm, "window|time_window|tw"))
}
if (needs_header_fix(raw)) {
  # tenter de promouvoir la première ligne en noms de colonnes
  hdr <- unlist(raw[1, ], use.names = FALSE)
  # si beaucoup de NA, ne pas forcer
  if (sum(is.na(hdr)) < length(hdr)) {
    names(raw) <- as.character(hdr)
    raw <- raw[-1,]
  }
}

df <- raw %>% janitor::clean_names()

# Unifier les fenêtres tôt si la colonne existe "en clair"
# (sinon ça sera fait après renommage)
if (any(str_detect(names(df), "^window$|^time_window$|^tw$"))) {
  wcol <- names(df)[str_detect(names(df), "^window$|^time_window$|^tw$")][1]
  df[[wcol]] <- .norm_window(df[[wcol]])
}

# -------------------- Gérer doublons et renommer proprement --------------------
# Many sheets doublent certaines colonnes : k/k_1, p/pval, etc.
df <- df %>%
  # coalescer 'k'
  coalesce_to(candidates = c("k","k_1","n","n_contrasts","studies"), out_name = "k_tmp") %>%
  # coalescer p
  coalesce_to(candidates = c("p","pval","p_value"), out_name = "p_tmp") %>%
  # coalescer OR
  coalesce_to(candidates = c("or","odds_ratio"), out_name = "or_tmp") %>%
  # coalescer log(OR)
  coalesce_to(candidates = c("log_or","logor","logor_1","yi","log_or_"), out_name = "logor_tmp") %>%
  # garder aussi les colonnes utiles si elles existent
  identity()

nm <- names(df)

molecule_col <- tryCatch(.find_col(nm, c("molecule","drug"), "molecule"), error = function(e) NA_character_)
ae_col       <- tryCatch(.find_col(nm, c("ae_term","ae","adverse_event","event"), "ae_term"), error = function(e) NA_character_)
win_col      <- tryCatch(.find_col(nm, c("window","time_window","tw"), "window"), error = function(e) NA_character_)

if (any(is.na(c(molecule_col, ae_col, win_col)))) {
  stop("Impossible d'identifier molecule/ae_term/window après nettoyage. Vérifie la première ligne de ton Excel.")
}

df1 <- df %>%
  rename(
    molecule = !!molecule_col,
    ae_term  = !!ae_col,
    window   = !!win_col,
    k        = k_tmp,
    OR       = or_tmp,
    logOR    = logor_tmp,
    p        = p_tmp
  ) %>%
  mutate(
    window   = .norm_window(window),
    molecule = as.character(molecule),
    ae_term  = as.character(ae_term),
    k        = suppressWarnings(as.numeric(k)),
    OR       = suppressWarnings(as.numeric(OR)),
    logOR    = suppressWarnings(as.numeric(logOR)),
    p        = suppressWarnings(as.numeric(p)),
    sig      = is.finite(p) & p < 0.05
  ) %>%
  filter(window %in% c("session","follow_up"))

# Sauvegarde intermédiaire (audit)
readr::write_csv(df1, file.path(OUT_DIR, "ae_sig_by_window_clean.csv"))

# -------------------- Statuts : transient / emergent / persistent --------------------
# Pour chaque (molecule, ae_term), regarder la significativité par fenêtre
sig_status <- df1 %>%
  select(molecule, ae_term, window, sig) %>%
  mutate(sig = ifelse(is.na(sig), FALSE, sig)) %>%
  distinct() %>%
  complete(molecule, ae_term, window = c("session","follow_up"), fill = list(sig = FALSE)) %>%
  group_by(molecule, ae_term) %>%
  summarise(
    sig_session  = any(sig[window == "session"], na.rm = TRUE),
    sig_followup = any(sig[window == "follow_up"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    status = case_when(
      sig_session & sig_followup ~ "Persistent",
      sig_session & !sig_followup ~ "Transient (session-only)",
      !sig_session & sig_followup ~ "Emergent (follow-up-only)",
      TRUE ~ "Non-significant"
    )
  )

# Comptes par molécule
counts_by_molecule <- sig_status %>%
  group_by(molecule) %>%
  summarise(
    `AE significant (session)`   = sum(sig_session,  na.rm = TRUE),
    `AE significant (follow-up)` = sum(sig_followup, na.rm = TRUE),
    `Transient (session-only)`   = sum(status == "Transient (session-only)"),
    `Emergent (follow-up-only)`  = sum(status == "Emergent (follow-up-only)"),
    `Persistent (both)`          = sum(status == "Persistent"),
    .groups = "drop"
  ) %>%
  arrange(molecule)

# Sauvegardes CSV
readr::write_csv(sig_status,         file.path(OUT_DIR, "ae_sig_transition_by_ae.csv"))
readr::write_csv(counts_by_molecule, file.path(OUT_DIR, "ae_sig_transition_counts.csv"))

# -------------------- Tableau LaTeX (comptes par molécule) --------------------
make_latex_table <- function(df_counts, outfile){
  # On fixe l’ordre LSD / MDMA / Psilocybin / Ayahuasca si présent
  mol_order <- c("LSD","MDMA","PSILOCYBIN","AYAHUASCA","Psilocybin","Ayahuasca")
  df_counts <- df_counts %>%
    mutate(molecule = factor(molecule, levels = unique(c(intersect(mol_order, molecule), setdiff(molecule, mol_order))))) %>%
    arrange(molecule)
  
  tot <- df_counts %>%
    summarise(
      molecule = "Total",
      `AE significant (session)`   = sum(`AE significant (session)`),
      `AE significant (follow-up)` = sum(`AE significant (follow-up)`),
      `Transient (session-only)`   = sum(`Transient (session-only)`),
      `Emergent (follow-up-only)`  = sum(`Emergent (follow-up-only)`),
      `Persistent (both)`          = sum(`Persistent (both)`)
    )
  
  df_tex <- bind_rows(df_counts, tot)
  
  # Génération LaTeX simple (booktabs)
  header <- "\\begin{table}[ht]\n\\centering\n\\caption{Adverse events (AEs) by time window and persistence status, per molecule. Counts reflect AE terms with p<0.05 in the forest models.}\n\\label{tab:ae_sig_persistence}\n\\begin{tabular}{lccccc}\n\\toprule\nMolecule & AE sig. (Session) & AE sig. (Follow-up) & Transient & Emergent & Persistent \\\\\n\\midrule\n"
  body <- paste0(
    apply(df_tex, 1, function(r){
      glue("{r[['molecule']]} & {r[['AE significant (session)']]} & {r[['AE significant (follow-up)']]} & {r[['Transient (session-only)']]} & {r[['Emergent (follow-up-only)']]} & {r[['Persistent (both)']]} \\\\")
    }),
    collapse = "\n"
  )
  footer <- "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n"
  cat(header, body, footer, file = outfile, sep = "")
}

make_latex_table(counts_by_molecule, file.path(OUT_DIR, "ae_sig_transition_counts.tex"))

message("✅ Fini. Fichiers écrits dans: ", OUT_DIR,
        "\n - ae_sig_by_window_clean.csv",
        "\n - ae_sig_transition_by_ae.csv",
        "\n - ae_sig_transition_counts.csv",
        "\n - ae_sig_transition_counts.tex")


suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(stringr); library(tidyr)
  library(janitor); library(knitr); library(kableExtra)
})

# ---- Parseur de p-values robuste (reprend celui qu’on a ajouté) ----
parse_p <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- gsub(",", ".", x)
  x <- trimws(x)
  x <- gsub("^p\\s*[<=>]\\s*", "", x, ignore.case = TRUE)
  x <- gsub("^<\\s*", "<", x)
  x <- gsub("^≤\\s*", "<=", x)
  x <- tolower(x)
  x[x %in% c("", "na", "n/a", "ns", "n.s.", "non sig", "non-significant")] <- NA
  is_lt  <- grepl("^<",  x)
  is_le  <- grepl("^<=", x)
  val <- suppressWarnings(as.numeric(gsub("^[^0-9.]*", "", x)))
  val[is_lt  & is.finite(val)] <- pmax(val[is_lt],  .Machine$double.xmin)
  val[is_le  & is.finite(val)] <- val[is_le]
  val
}

# ---- Fonction utilitaire: étoiles ---
p_to_stars <- function(p) {
  case_when(
    is.na(p)        ~ "",
    p < 0.001       ~ "***",
    p < 0.01        ~ "**",
    p < 0.05        ~ "*",
    TRUE            ~ ""
  )
}

# ---- GÉNÉRATION DES TABLEAUX PUBLICATION ----
make_sig_tables <- function(INPUT_XLS,
                            OUT_DIR,
                            p_thresh = 0.05,
                            top_n_per_cell = NULL  # optionnel: limite #AE listés par cellule
) {
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # 1) Lecture + nettoyage
  df_raw <- readxl::read_excel(INPUT_XLS) %>% clean_names()
  
  # On s’assure des colonnes attendues (noms du fichier publication)
  required <- c("molecule","ae_term","window","k","or","log_or","p","stars")
  miss <- setdiff(required, names(df_raw))
  if (length(miss)) stop("Colonnes manquantes dans l'input: ", paste(miss, collapse=", "))
  
  df <- df_raw %>%
    mutate(
      molecule = toupper(str_squish(molecule)),
      ae_term  = str_squish(ae_term),
      window   = tolower(str_replace_all(str_squish(window), "[-_/]", " ")),
      window   = case_when(
        grepl("\\bsession\\b|\\bacute\\b|day 0|during", window) ~ "session",
        grepl("follow ?up|post|week|month|after", window)        ~ "follow_up",
        TRUE ~ window
      ),
      p      = parse_p(p),
      stars  = ifelse(is.na(stars) | stars == "", p_to_stars(p), stars)
    ) %>%
    filter(window %in% c("session","follow_up"))
  
  # 2) Filtre significatif + consolidation AE (si doublons)
  sig_df <- df %>%
    filter(is.finite(p), p < p_thresh) %>%
    group_by(molecule, window, ae_term) %>%
    summarise(
      k       = sum(suppressWarnings(as.numeric(k)), na.rm = TRUE),
      or_med  = median(suppressWarnings(as.numeric(or)), na.rm = TRUE),
      p_min   = suppressWarnings(min(p, na.rm = TRUE)),
      stars   = p_to_stars(p_min)
      , .groups = "drop") %>%
    arrange(molecule, window, p_min, desc(abs(log(or_med))))
  
  # 3) TABLEAU 1 — AE significatifs par molécule × fenêtre (large)
  listify <- function(x, stars, or_med) {
    # AE ± étoiles ; on ajoute OR arrondi pour aider la lecture
    lab <- paste0(x, " ", stars, " (OR≈", formatC(or_med, format="f", digits=2), ")")
    lab
  }
  
  sig_wide <- sig_df %>%
    mutate(ae_label = listify(ae_term, stars, or_med)) %>%
    group_by(molecule, window) %>%
    summarise(ae_list = {
      v <- ae_label
      if (!is.null(top_n_per_cell) && length(v) > top_n_per_cell) {
        v <- v[seq_len(top_n_per_cell)]
        v <- c(v, "…")
      }
      paste(v, collapse = "; ")
    }, .groups = "drop") %>%
    tidyr::pivot_wider(names_from = window, values_from = ae_list) %>%
    arrange(molecule)
  
  # Sauvegardes
  write.csv(sig_df,   file.path(OUT_DIR, "significant_ae_long.csv"), row.names = FALSE)
  write.csv(sig_wide, file.path(OUT_DIR, "significant_ae_by_molecule_window.csv"), row.names = FALSE)
  
  # 4) TABLEAU 2 — Statut Transient / Persistent / Emergent
  status_df <- sig_df %>%
    mutate(is_sig = TRUE) %>%
    select(molecule, window, ae_term, is_sig) %>%
    tidyr::pivot_wider(names_from = window, values_from = is_sig, values_fill = FALSE) %>%
    mutate(
      status = case_when(
        session & follow_up ~ "Persistent",
        session & !follow_up ~ "Transient (session-only)",
        !session & follow_up ~ "Emergent (follow-up-only)",
        TRUE ~ "Non-significant"  # ne devrait pas arriver ici
      )
    )
  
  status_wide <- status_df %>%
    group_by(molecule, status) %>%
    summarise(ae_list = paste(sort(ae_term), collapse = "; "), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = status, values_from = ae_list) %>%
    arrange(molecule)
  
  write.csv(status_df,   file.path(OUT_DIR, "ae_transition_status_long.csv"), row.names = FALSE)
  write.csv(status_wide, file.path(OUT_DIR, "ae_transition_status_by_molecule.csv"), row.names = FALSE)
  
  # 5) Rendus LaTeX (via kable)
  # Table 1
  tab1_tex <- sig_wide %>%
    rename(`Molecule` = molecule, `Session (p<0.05)` = session, `Follow-up (p<0.05)` = follow_up) %>%
    kable(format = "latex", booktabs = TRUE, escape = TRUE,
          caption = "Significant adverse events by molecule and time window (session vs. follow-up). Entries list AE with significance stars and approximate pooled OR.",
          label = "tab:ae-signif-by-window",
          col.names = c("Molecule", "Session (p<0.05)", "Follow-up (p<0.05)")) %>%
    kable_styling(latex_options = c("hold_position"))
  
  cat(tab1_tex, file = file.path(OUT_DIR, "significant_ae_by_molecule_window.tex"))
  
  # Table 2
  # S’assure que toutes les colonnes existent :
  for (col in c("Persistent","Transient (session-only)","Emergent (follow-up-only)")) {
    if (!col %in% names(status_wide)) status_wide[[col]] <- ""
  }
  tab2_tex <- status_wide %>%
    rename(`Molecule` = molecule) %>%
    select(`Molecule`, `Transient (session-only)`, `Emergent (follow-up-only)`, `Persistent`) %>%
    kable(format = "latex", booktabs = TRUE, escape = TRUE,
          caption = "Temporal status of significant AEs by molecule: Transient (session-only), Emergent (follow-up-only), or Persistent (both windows).",
          label = "tab:ae-temporal-status",
          col.names = c("Molecule", "Transient (session-only)", "Emergent (follow-up-only)", "Persistent")) %>%
    kable_styling(latex_options = c("hold_position"))
  
  cat(tab2_tex, file = file.path(OUT_DIR, "ae_transition_status_by_molecule.tex"))
  
  message("✔ Tables saved in: ", OUT_DIR)
  invisible(list(sig_long = sig_df, sig_wide = sig_wide, status_long = status_df, status_wide = status_wide))
}
