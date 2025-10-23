suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(kableExtra)
})

# ==== 1) Paths ====
INFILE <- "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/data/Adverse-events-dose-v5.xlsx"
SHEET  <- "Feuil1"
OUTDIR <- "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/results/tables_paper"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ==== 2) Helpers ====
clean_chr <- function(x) tolower(trimws(as.character(x)))

collapse_unique <- function(x, nmax = 4) {
  x <- unique(na.omit(x))
  if (!length(x)) return(NA_character_)
  x <- sort(x)
  if (length(x) <= nmax) paste(x, collapse = ", ") else paste0(paste(x[1:nmax], collapse = ", "), ", …")
}

# Detect if an arm is placebo (using your French scheme + substance)
is_placebo <- function(arm_id, dose_mg, placebo_substance) {
  arm_id <- clean_chr(arm_id)
  subs   <- clean_chr(placebo_substance)
  is_zero <- suppressWarnings(!is.na(dose_mg) & as.numeric(dose_mg) == 0)
  
  # Keywords in arm_id
  has_placebo_kw <- str_detect(arm_id, "placebo|placebo_inactif|placebo_actif")
  
  # Substance-based hints
  inactive_subs   <- c("inactive","inactif","lactose","mannitol","sugar","placebo")
  active_nonpsy   <- c("niacin","vit_b3","vitamine_b3")
  active_psy_subs <- c("lsd","mdma","psilocybin","psilocybine","ayahuasca")
  
  is_inactive     <- !is.na(subs) & subs %in% inactive_subs
  is_act_nonpsy   <- !is.na(subs) & subs %in% active_nonpsy
  is_act_psy      <- !is.na(subs) & subs %in% active_psy_subs
  
  is_zero || has_placebo_kw || is_inactive || is_act_nonpsy || is_act_psy
}

# Classify placebo type (French labels you asked)
placebo_type <- function(arm_id, dose_mg, placebo_substance) {
  arm_id <- clean_chr(arm_id)
  subs   <- clean_chr(placebo_substance)
  is_zero <- suppressWarnings(!is.na(dose_mg) & as.numeric(dose_mg) == 0)
  
  inactive_subs   <- c("inactive","inactif","lactose","mannitol","sugar","placebo")
  active_nonpsy   <- c("niacin","vit_b3","vitamine_b3")
  active_psy_subs <- c("lsd","mdma","psilocybin","psilocybine","ayahuasca")
  
  case_when(
    # Substance overrides
    !is.na(subs) & subs %in% inactive_subs   ~ "placebo_inactif",
    !is.na(subs) & subs %in% active_nonpsy   ~ "placebo_actif_non_psychedelique",
    !is.na(subs) & subs %in% active_psy_subs ~ "placebo_actif_psychedelique",
    
    # Name-based fallbacks
    str_detect(arm_id, "placebo_inactif")    ~ "placebo_inactif",
    str_detect(arm_id, "placebo_?actif")     ~ "placebo_actif_non_psychedelique",
    is_zero & str_detect(arm_id, "placebo")  ~ "placebo_inactif",
    
    # If our general detector says placebo but no class matched
    is_placebo(arm_id, dose_mg, placebo_substance) ~ "placebo_non_precise",
    
    TRUE ~ NA_character_
  )
}

# ==== 3) Load & harmonise ====
stopifnot(file.exists(INFILE))
raw <- readxl::read_excel(INFILE, sheet = SHEET)

required <- c("study_id","author_year","arm_id","molecule","n_participants_arm",
              "ae_term","time_window","events","dose_mg","placebo_substance")
miss <- setdiff(required, names(raw))
if (length(miss)) stop("Colonnes manquantes dans l’Excel: ", paste(miss, collapse=", "))

df <- raw %>%
  select(all_of(required)) %>%
  mutate(
    molecule    = toupper(trimws(molecule)),
    time_window = tolower(gsub("[-_/\\s]+", "_", trimws(time_window))),
    dose_mg     = suppressWarnings(as.numeric(dose_mg)),
    is_placebo  = is_placebo(arm_id, dose_mg, placebo_substance),
    placebo_type = placebo_type(arm_id, dose_mg, placebo_substance),
    is_active   = !is_placebo
  )

# Arm-level de-dup (avoid double counting across AEs)
arm_level <- df %>%
  distinct(study_id, author_year, molecule, arm_id, dose_mg, n_participants_arm,
           time_window, is_placebo, placebo_type, placebo_substance, is_active)

# ==== 4) Aggregated table (per molecule) ====
agg_by_molecule <- arm_level %>%
  group_by(molecule) %>%
  summarise(
    n_studies       = n_distinct(study_id),
    participants_arms = sum(unique(n_participants_arm)),      # sum across unique arms
    dose_min_active = suppressWarnings(min(dose_mg[is_active], na.rm = TRUE)),
    dose_max_active = suppressWarnings(max(dose_mg[is_active], na.rm = TRUE)),
    n_unique_doses  = n_distinct(dose_mg[is_active]),
    n_placebo_inactif         = sum(placebo_type == "placebo_inactif", na.rm = TRUE),
    n_placebo_actif_non_psy   = sum(placebo_type == "placebo_actif_non_psychedelique", na.rm = TRUE),
    n_placebo_actif_psy       = sum(placebo_type == "placebo_actif_psychedelique", na.rm = TRUE),
    n_placebo_non_precise     = sum(placebo_type == "placebo_non_precise", na.rm = TRUE),
    subs_actif_non_psy = collapse_unique(placebo_substance[placebo_type == "placebo_actif_non_psychedelique"]),
    subs_inactif       = collapse_unique(placebo_substance[placebo_type == "placebo_inactif"]),
    subs_actif_psy     = collapse_unique(placebo_substance[placebo_type == "placebo_actif_psychedelique"]),
    has_session   = any(time_window == "session"),
    has_follow_up = any(time_window == "follow_up"),
    .groups = "drop"
  ) %>%
  mutate(
    dose_range_active = ifelse(
      is.finite(dose_min_active) & is.finite(dose_max_active),
      paste0(format(dose_min_active, trim=TRUE), "–", format(dose_max_active, trim=TRUE)),
      "—"
    ),
    windows = paste0(ifelse(has_session, "Session", "—"), " / ",
                     ifelse(has_follow_up, "Follow-up", "—"))
  ) %>%
  select(
    molecule, n_studies, participants_arms,
    dose_range_active, n_unique_doses,
    n_placebo_inactif, n_placebo_actif_non_psy, n_placebo_actif_psy, n_placebo_non_precise,
    subs_actif_non_psy, subs_inactif, subs_actif_psy, windows
  ) %>%
  arrange(molecule)

# ==== 5) Outcome availability (per molecule × window) ====
outcome_availability <- df %>%
  group_by(molecule, time_window) %>%
  summarise(
    n_studies   = n_distinct(study_id),
    n_arms      = n_distinct(arm_id),
    n_unique_ae = n_distinct(ae_term),
    .groups = "drop"
  ) %>%
  mutate(time_window = recode(time_window, session="Session", follow_up="Follow-up", .default = time_window)) %>%
  arrange(molecule, time_window)

# ==== 6) Placebo distribution tidy (counts + row %) ====
placebo_dist <- arm_level %>%
  filter(!is.na(placebo_type)) %>%
  count(molecule, placebo_type, name = "n_arms") %>%
  group_by(molecule) %>%
  mutate(pct = 100 * n_arms / sum(n_arms)) %>%
  ungroup() %>%
  arrange(molecule, desc(n_arms))

# ==== 7) Export CSV ====
write_csv(agg_by_molecule,       file.path(OUTDIR, "study_characteristics_aggregated.csv"))
write_csv(outcome_availability,  file.path(OUTDIR, "outcome_availability_by_window.csv"))
write_csv(placebo_dist,          file.path(OUTDIR, "placebo_distribution_by_molecule.csv"))

# ==== 8) Export LaTeX ====

# 8a) Main characteristics table (FR headers)
latex_main <- agg_by_molecule %>%
  rename(
    `Placebo (inactif)`              = n_placebo_inactif,
    `Placebo (actif non-psy)`        = n_placebo_actif_non_psy,
    `Placebo (actif psychédélique)`  = n_placebo_actif_psy,
    `Placebo (non précisé)`          = n_placebo_non_precise,
    `Substances actives (non-psy)`   = subs_actif_non_psy,
    `Substances inactives`           = subs_inactif,
    `Substances actives (psy)`       = subs_actif_psy
  ) %>%
  select(
    Molécule = molecule,
    Études = n_studies,
    `Participants (bras)` = participants_arms,
    `Plage de dose (mg, actifs)` = dose_range_active,
    `# doses actives` = n_unique_doses,
    `Placebo (inactif)`,
    `Placebo (actif non-psy)`,
    `Placebo (actif psychédélique)`,
    `Placebo (non précisé)`,
    `Substances actives (non-psy)`,
    `Substances inactives`,
    `Substances actives (psy)`,
    Fenêtres = windows
  ) %>%
  kbl(format = "latex", booktabs = TRUE, escape = TRUE,
      caption = "Caractéristiques par molécule : participants (bras uniques), plage de dose active, répartition des placebos et fenêtres disponibles.",
      label = "tab:caracteristiques_etudes_par_molecule",
      align = "lrrrrrrrrlll") %>%
  kable_styling(latex_options = c("hold_position", "striped")) %>%
  column_spec(4, width = "2.5cm") %>%
  column_spec(10:12, width = "3.2cm") %>%
  column_spec(13, width = "2.2cm")

writeLines(latex_main, file.path(OUTDIR, "study_characteristics_aggregated.tex"))

# 8b) Outcome availability table
latex_outcome <- outcome_availability %>%
  select(
    Molécule = molecule,
    Fenêtre  = time_window,
    Études   = n_studies,
    Bras     = n_arms,
    `Termes AE uniques` = n_unique_ae
  ) %>%
  kbl(format="latex", booktabs=TRUE, escape=TRUE,
      caption="Disponibilité des résultats par molécule et fenêtre (Session vs Follow-up) : nombre d’études, de bras et de termes AE uniques.",
      label="tab:disponibilite_outcomes_par_fenetre",
      align="llrrr") %>%
  kable_styling(latex_options=c("hold_position","striped")) %>%
  column_spec(5, width="2.8cm")

writeLines(latex_outcome, file.path(OUTDIR, "outcome_availability_by_window.tex"))

# 8c) Placebo distribution (counts + %)
latex_placebo <- placebo_dist %>%
  mutate(
    placebo_type = recode(placebo_type,
                          "placebo_inactif" ~ "Placebo (inactif)",
                          "placebo_actif_non_psychedelique" ~ "Placebo (actif non-psy)",
                          "placebo_actif_psychedelique" ~ "Placebo (actif psychédélique)",
                          "placebo_non_precise" ~ "Placebo (non précisé)",
                          .default = placebo_type
    ),
    pct = sprintf("%.1f%%", pct)
  ) %>%
  arrange(molecule, placebo_type) %>%
  kbl(format="latex", booktabs=TRUE, escape=TRUE,
      caption="Répartition des types de placebo par molécule (nombre de bras et pourcentages en ligne).",
      label="tab:repartition_placebos_par_molecule",
      col.names = c("Molécule","Type de placebo","Bras (n)","% (ligne)"),
      align="llrc") %>%
  kable_styling(latex_options=c("hold_position","striped"))

writeLines(latex_placebo, file.path(OUTDIR, "placebo_distribution_by_molecule.tex"))

message("✅ Descriptive tables written to: ", OUTDIR)