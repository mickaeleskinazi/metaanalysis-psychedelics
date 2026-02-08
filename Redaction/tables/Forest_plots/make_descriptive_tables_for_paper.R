suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(kableExtra)
  library(readr)
})

# ========= 0) PARAMS =========
INFILE <- "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/data/Adverse-events-dose-v5.xlsx"
SHEET  <- "Feuil1"   # adapte si besoin
OUTDIR <- "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/results/tables_paper"

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ========= 1) HELPERS (vectorisés & robustes) =========

# a) Normalisation des fenêtres temporelles
norm_window <- function(x){
  x <- tolower(gsub("[-_/\\s]+", "_", trimws(as.character(x))))
  dplyr::recode(
    x,
    "followup"  = "follow_up",
    "follow-up" = "follow_up",
    .default    = x
  )
}

# b) Normalisation des molécules en UPPER (pour cohérence)
norm_molecule <- function(x) toupper(trimws(as.character(x)))

# c) Détection placebo (vectorisée)
is_placebo_vec <- function(arm_id, dose_mg, placebo_substance){
  arm_id <- tolower(as.character(arm_id))
  subs   <- tolower(as.character(placebo_substance))
  
  # zéro dose considéré placebo
  is_zero <- suppressWarnings(!is.na(as.numeric(dose_mg)) & as.numeric(dose_mg) == 0)
  
  # mots-clés explicites
  has_kw  <- str_detect(arm_id, "placebo") |
    str_detect(arm_id, "placebo_?actif|active_placebo|placebo_inactif")
  
  # substances indiquées
  has_subs <- !is.na(subs) & nzchar(subs)
  
  is_zero | has_kw | has_subs
}

# d) Typage placebo (tes libellés FR, + fallback)
#    Renvoie l'une des étiquettes suivantes:
#    "placebo_inactif" | "placebo_actif_non_psychedelique" | "placebo_actif_psychedelique" | "placebo_non_precise" | NA
placebo_type_vec <- function(arm_id, dose_mg, placebo_substance){
  arm_id <- tolower(as.character(arm_id))
  subs   <- tolower(as.character(placebo_substance))
  
  # listes de mots-clés côté substance (à étendre si besoin)
  inactive_kw <- c("inactive","lactose","mannitol","sugar","placebo")
  active_nonpsy_kw <- c("niacin")
  active_psy_kw <- c("lsd","mdma","psilocybin","psilocybine","ayahuasca","psilo")
  
  is_zero <- suppressWarnings(!is.na(as.numeric(dose_mg)) & as.numeric(dose_mg) == 0)
  has_kw  <- str_detect(arm_id, "placebo")
  
  out <- rep(NA_character_, length(arm_id))
  
  # cas où une substance est renseignée
  has_subs <- !is.na(subs) & nzchar(subs)
  out[has_subs & subs %in% inactive_kw]     <- "placebo_inactif"
  out[has_subs & subs %in% active_nonpsy_kw]<- "placebo_actif_non_psychedelique"
  out[has_subs & subs %in% active_psy_kw]   <- "placebo_actif_psychedelique"
  
  # fallback par arm_id (si pas de substance)
  out[is.na(out) & str_detect(arm_id, "placebo_inactif")] <- "placebo_inactif"
  out[is_na(out) & str_detect(arm_id, "placebo_?actif|active_placebo")] <- "placebo_actif_non_psychedelique"
  
  # si dose=0 & "placebo" dans le nom et rien déterminé → inactif par défaut
  out[is.na(out) & is_zero & has_kw] <- "placebo_inactif"
  
  # si identifié comme placebo mais sans type déterminé
  is_pl <- is_placebo_vec(arm_id, dose_mg, subs)
  out[is_pl & is.na(out)] <- "placebo_non_precise"
  
  out
}

# e) Helper pour concaténer des valeurs uniques (pour colonnes "substances")
collapse_unique <- function(x, nmax = 4){
  x <- unique(na.omit(x))
  if (!length(x)) return(NA_character_)
  x <- sort(x)
  if (length(x) <= nmax) paste(x, collapse=", ")
  else paste0(paste(x[1:nmax], collapse=", "), ", …")
}


# --- Robust time-window normalization -----------------------------------------
norm_window <- function(x){
  # 1) coerce to plain ASCII, lowercase, squish
  y <- as.character(x)
  y <- iconv(y, to = "ASCII//TRANSLIT")          # remove accents safely
  y <- tolower(y)
  y <- gsub("[^a-z0-9]+", "_", y)                # non-alnum -> underscore
  y <- gsub("^_+|_+$", "", y)                    # trim underscores
  
  # 2) fix common corruptions and synonyms
  y <- gsub("^_?e_ion$", "session", y)           # fixes "_e_ion" → "session"
  y <- gsub("^seance$", "session", y)            # FR → EN
  y <- gsub("^acute$", "session", y)
  
  # 3) final mapping
  out <- dplyr::case_when(
    grepl("\\b(session|day_?0)\\b", y)                  ~ "session",
    grepl("follow[_-]?up|post|week|month|\\bfu\\b", y)  ~ "follow_up",
    TRUE                                                ~ y
  )
  out
}

# Apply to your dataframe `df` right after reading/renaming:
df <- df %>%
  mutate(time_window = norm_window(time_window))

# Sanity check:
message("Distinct time_window values after normalization:")
print(sort(unique(df$time_window)))

# ========= 2) LECTURE & HARMONISATION =========
req <- c("study_id","author_year","arm_id","molecule",
         "n_participants_arm","ae_term","time_window",
         "events","dose_mg","placebo_substance")

raw <- readxl::read_excel(INFILE, sheet = SHEET)
missing <- setdiff(req, names(raw))
if (length(missing)) stop("Colonnes manquantes dans l’Excel: ", paste(missing, collapse=", "))

df <- raw %>%
  select(all_of(req)) %>%
  mutate(
    molecule     = norm_molecule(molecule),
    time_window  = norm_window(time_window),
    dose_mg      = suppressWarnings(as.numeric(dose_mg)),
    n_participants_arm = suppressWarnings(as.numeric(n_participants_arm)),
    # flags placebo/type
    is_placebo   = is_placebo_vec(arm_id, dose_mg, placebo_substance),
    placebo_type = placebo_type_vec(arm_id, dose_mg, placebo_substance),
    is_active    = !is_placebo
  )

# version "bras uniques" (anti double-comptage via AE)
arm_level <- df %>%
  distinct(study_id, author_year, molecule, arm_id,
           dose_mg, n_participants_arm, time_window,
           is_placebo, placebo_type, placebo_substance, is_active)

# ========= 3) TABLES DE BASE =========

# 3a) Participants par molécule (2 façons)
# - participants_arms_sum : somme brute des tailles de bras (unique arms)
# - participants_unique_study : somme des max(bras) par étude (meilleure approx sans double-compte)
participants_by_molecule <- arm_level %>%
  group_by(molecule, study_id) %>%
  summarise(max_n_study = suppressWarnings(max(n_participants_arm, na.rm = TRUE)), .groups="drop") %>%
  group_by(molecule) %>%
  summarise(
    participants_unique_study = sum(max_n_study, na.rm = TRUE),
    .groups="drop"
  ) %>%
  left_join(
    arm_level %>% group_by(molecule) %>%
      summarise(participants_arms_sum = sum(n_participants_arm, na.rm=TRUE), .groups="drop"),
    by="molecule"
  )

# 3b) Dose ranges (actives uniquement)
dose_by_molecule <- arm_level %>%
  filter(is_active, !is.na(dose_mg)) %>%
  group_by(molecule) %>%
  summarise(
    dose_min_active   = suppressWarnings(min(dose_mg, na.rm = TRUE)),
    dose_max_active   = suppressWarnings(max(dose_mg, na.rm = TRUE)),
    n_unique_doses    = n_distinct(dose_mg),
    .groups = "drop"
  ) %>%
  mutate(dose_range_active = ifelse(is.finite(dose_min_active) & is.finite(dose_max_active),
                                    paste0(format(dose_min_active, trim=TRUE), "–", format(dose_max_active, trim=TRUE)),
                                    "—"))

# 3c) Placebos par **étude** (pas par bras) — pour éviter de compter 2x session/follow_up
placebos_per_study <- arm_level %>%
  filter(is_placebo) %>%
  group_by(molecule, study_id) %>%
  summarise(
    placebo_types_in_study = collapse_unique(placebo_type, nmax = 10),
    placebo_subs_in_study  = collapse_unique(placebo_substance, nmax = 10),
    .groups="drop"
  )

# agrégation par molécule (par étude)
placebo_distribution_by_study <- placebos_per_study %>%
  group_by(molecule, placebo_types_in_study) %>%
  summarise(n_studies_with_type = n_distinct(study_id), .groups="drop") %>%
  arrange(molecule, placebo_types_in_study)

# 3d) Fenêtres disponibles & richesse AE / bras
by_window_arms <- arm_level %>%
  group_by(molecule, time_window) %>%
  summarise(
    n_studies      = n_distinct(study_id),
    n_arms         = n(),  # bras déjà distincts dans arm_level
    n_active_arms  = sum(is_active, na.rm = TRUE),
    n_active_doses = n_distinct(dose_mg[is_active & !is.na(dose_mg)]),
    .groups="drop"
  )

by_window_ae <- df %>%
  group_by(molecule, time_window) %>%
  summarise(n_unique_ae = n_distinct(ae_term), .groups = "drop")

by_window <- by_window_arms %>%
  left_join(by_window_ae, by = c("molecule","time_window"))

# ========= 4) GRAND TABLEAU AGRÉGÉ =========

agg_by_molecule <- arm_level %>%
  group_by(molecule) %>%
  summarise(
    n_studies = n_distinct(study_id),
    .groups="drop"
  ) %>%
  left_join(participants_by_molecule, by="molecule") %>%
  left_join(dose_by_molecule, by="molecule") %>%
  # repérage des fenêtres dispo
  left_join(
    by_window %>%
      mutate(has_win = TRUE) %>%
      select(molecule, time_window, has_win) %>%
      distinct() %>%
      pivot_wider(id_cols = molecule, names_from = time_window, values_from = has_win, values_fill = FALSE),
    by="molecule"
  ) %>%
  # Jolis labels de fenêtres
  mutate(windows = paste0(ifelse(session, "Session", "—"), " / ",
                          ifelse(follow_up, "Follow-up", "—"))) %>%
  # quelques colonnes de placebo (par étude)
  left_join(
    placebos_per_study %>%
      group_by(molecule) %>%
      summarise(
        studies_with_placebo = n_distinct(study_id),
        example_active_placebo_subs   = collapse_unique(placebo_subs_in_study[str_detect(placebo_types_in_study, "actif")], nmax=4),
        example_inactive_placebo_subs = collapse_unique(placebo_subs_in_study[placebo_types_in_study=="placebo_inactif"], nmax=4),
        .groups="drop"
      ),
    by="molecule"
  ) %>%
  # colonnes finales ordonnées
  transmute(
    molecule,
    n_studies,
    participants_unique_study,
    participants_arms_sum,
    dose_range_active,
    n_unique_doses,
    windows,
    example_inactive_placebo_subs,
    example_active_placebo_subs
  ) %>%
  arrange(molecule)

# ========= 5) TABLE WIDE par FENÊTRE =========

by_window_wide <- by_window %>%
  pivot_wider(
    id_cols = molecule,
    names_from  = time_window,
    values_from = c(n_studies, n_arms, n_unique_ae, n_active_arms, n_active_doses),
    names_sep = "_",
    values_fill = list(
      n_studies=0L, n_arms=0L, n_unique_ae=0L, n_active_arms=0L, n_active_doses=0L
    )
  )

# assurer toutes les colonnes
need_cols <- c("n_studies_session","n_studies_follow_up",
               "n_arms_session","n_arms_follow_up",
               "n_unique_ae_session","n_unique_ae_follow_up",
               "n_active_arms_session","n_active_arms_follow_up",
               "n_active_doses_session","n_active_doses_follow_up")
for (cc in need_cols) if (!cc %in% names(by_window_wide)) by_window_wide[[cc]] <- 0L

final_by_window <- by_window_wide %>%
  transmute(
    Molecule = molecule,
    `Studies (session)`           = n_studies_session,
    `Studies (follow-up)`         = n_studies_follow_up,
    `Arms (session)`              = n_arms_session,
    `Arms (follow-up)`            = n_arms_follow_up,
    `Unique AE terms (session)`   = n_unique_ae_session,
    `Unique AE terms (follow-up)` = n_unique_ae_follow_up,
    `Active arms (session)`       = n_active_arms_session,
    `Active arms (follow-up)`     = n_active_arms_follow_up,
    `Active doses (session)`      = n_active_doses_session,
    `Active doses (follow-up)`    = n_active_doses_follow_up
  ) %>%
  arrange(Molecule)

# ========= 6) EXPORT CSV (optionnel) =========
write_csv(agg_by_molecule, file.path(OUTDIR, "study_characteristics_aggregated.csv"))
write_csv(final_by_window, file.path(OUTDIR, "per_molecule_by_window.csv"))
write_csv(placebo_distribution_by_study, file.path(OUTDIR, "placebo_distribution_by_study.csv"))

# ========= 7) EXPORT LaTeX =========
# (a) Grand tableau agrégé
latex_main <- agg_by_molecule %>%
  select(
    Molecule = molecule,
    `Studies` = n_studies,
    `Participants (unique/study)` = participants_unique_study,
    `Participants (arms sum)`     = participants_arms_sum,
    `Dose range (mg, active)`     = dose_range_active,
    `# active doses`              = n_unique_doses,
    `Windows`                     = windows,
    `Inactive placebo subs.`      = example_inactive_placebo_subs,
    `Active placebo subs.`        = example_active_placebo_subs
  ) %>%
  kbl(format="latex", booktabs=TRUE, escape=TRUE,
      caption="Study-level characteristics by molecule: studies, participants, active dose range, and windows. Placebo substances are illustrative examples aggregated per study.",
      label="tab:study_characteristics_aggregated") %>%
  kable_styling(latex_options=c("hold_position","striped")) %>%
  column_spec(5, width="2.5cm") %>%
  column_spec(8:9, width="3.2cm")
writeLines(latex_main, file.path(OUTDIR, "study_characteristics_aggregated.tex"))

# (b) Dispo outcomes par fenêtre
latex_bywin <- final_by_window %>%
  kbl(format="latex", booktabs=TRUE, escape=TRUE,
      caption="Outcome availability by molecule and time window: counts of studies, arms, unique AE terms, active arms, and active doses.",
      label="tab:outcome_availability_by_window") %>%
  kable_styling(latex_options=c("hold_position","striped"))
writeLines(latex_bywin, file.path(OUTDIR, "outcome_availability_by_window.tex"))

# (c) Répartition des placebos PAR ÉTUDE
latex_placebo <- placebo_distribution_by_study %>%
  select(Molecule=molecule, `Placebo type(s) in study`=placebo_types_in_study,
         `# studies`=n_studies_with_type) %>%
  kbl(format="latex", booktabs=TRUE, escape=TRUE,
      caption="Distribution of placebo types by molecule (counted at the study level to avoid double-counting across windows).",
      label="tab:placebo_distribution_by_study") %>%
  kable_styling(latex_options=c("hold_position","striped"))
writeLines(latex_placebo, file.path(OUTDIR, "placebo_distribution_by_study.tex"))

# (d) Tableau “méga” (tout-en-un) si besoin d’un seul tableau compact
mega <- agg_by_molecule %>%
  left_join(final_by_window %>% rename(molecule = Molecule), by="molecule") %>%
  select(
    Molecule = molecule,
    Studies = n_studies,
    `Participants (unique/study)` = participants_unique_study,
    `Dose range (mg, active)` = dose_range_active,
    `# active doses` = n_unique_doses,
    Windows = windows,
    `Studies (session)`, `Studies (follow-up)`,
    `Arms (session)`, `Arms (follow-up)`,
    `Unique AE terms (session)`, `Unique AE terms (follow-up)`
  )
latex_mega <- mega %>%
  kbl(format="latex", booktabs=TRUE, escape=TRUE,
      caption="Compact “mega” table: core study characteristics plus per-window availability.",
      label="tab:mega_study_profile") %>%
  kable_styling(latex_options=c("hold_position","striped"))
writeLines(latex_mega, file.path(OUTDIR, "mega_study_profile.tex"))

message("✅ LaTeX tables saved in: ", OUTDIR)