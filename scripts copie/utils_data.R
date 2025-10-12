suppressPackageStartupMessages({
<<<<<<< Updated upstream
  library(dplyr); library(janitor); library(rlang); library(stringr)
})

.detect_and_rename <- function(df, targets){
  cn <- names(df); lower <- tolower(cn)
  getcol <- function(cands){
    # exact first, then contains
    ix <- match(cands, lower, nomatch = 0L); ix <- ix[ix > 0]
    if (length(ix)) return(cn[ix[1]])
    for (c in cands){
      hit <- which(grepl(c, lower, fixed = TRUE))
      if (length(hit)) return(cn[hit[1]])
    }
    NA_character_
  }
  for (std in names(targets)){
    found <- getcol(targets[[std]])
    if (!is.na(found) && found != std && !(std %in% names(df))){
      df <- dplyr::rename(df, !!std := !!sym(found))
    }
  }
  df
}

load_data <- function(path, sheet){
  stopifnot(file.exists(path))
  raw <- readxl::read_excel(path, sheet = sheet)
  df  <- janitor::clean_names(raw)
  
  # Map typical aliases to our standard schema
  targets <- list(
    study_id    = c("study_id","study","studyid","study_code","study_name"),
    arm_id      = c("arm_id","arm","group_id","condition_id"),
    group       = c("group","arm_label","armname","armlabel","group_name","condition","treatment"),
    arm_type    = c("arm_type","placebo_type","control_type","group_type"),
    molecule    = c("molecule","drug","substance","compound"),
    ae_term     = c("ae_term","adverse_event","ae","outcome","event_name"),
    time_window = c("time_window","window","timepoint","time","visit"),
    dose_mg     = c("dose_mg","dose","dose_mg_numeric","dose_mg_num","dose_milligram","lsd_dose_mg"),
    events      = c("events","event","ae_n","cases","num_events","n_events"),
    n           = c("n","total","n_total","sample_size","denominator")
  )
  df <- .detect_and_rename(df, targets)
  
  # Friendly error if must-haves are missing (except group; we can synthesize it)
  must_have <- c("study_id","molecule","ae_term","time_window","dose_mg","events","n")
  missing <- setdiff(must_have, names(df))
  if (length(missing)){
    cat("Columns in sheet:\n"); print(names(df))
    stop("Missing required columns after auto-detect: ", paste(missing, collapse = ", "))
  }
  
  # Normalize basic types
  df <- df %>%
    mutate(
      study_id    = as.character(study_id),
      molecule    = toupper(as.character(molecule)),
      ae_term     = as.character(ae_term),
      time_window = as.character(coalesce(time_window, "session")),
      dose_mg     = suppressWarnings(as.numeric(gsub(",", ".", as.character(dose_mg)))),
      events      = suppressWarnings(as.numeric(events)),
      n           = suppressWarnings(as.numeric(n))
    )
  
  # --- NEW: derive group + arm_type from arm_id when needed ---
  # 1) start from any provided group; otherwise fall back to arm_id
  df <- df %>%
    mutate(
      arm_id = if ("arm_id" %in% names(df)) as.character(arm_id) else NA_character_,
      group  = dplyr::coalesce(as.character(if ("group" %in% names(df)) group else NA_character_), arm_id)
    )
  
  # 2) derive arm_type:
  #    - if we already have arm_type, keep but normalize FR/EN
  #    - else infer:
  #         * contains "placebo"  -> inactive if dose_mg==0, active if dose_mg>0
  #         * otherwise           -> "active"
  df <- df %>%
    mutate(
      arm_type = if ("arm_type" %in% names(df)) tolower(as.character(arm_type)) else NA_character_,
      arm_type = case_when(
        !is.na(arm_type) & arm_type %in% c("placebo_inactif","inactive_placebo","placebo","inactif") ~ "inactive_placebo",
        !is.na(arm_type) & arm_type %in% c("placebo_actif","active_placebo","placebo actif") ~ "active_placebo",
        !is.na(arm_type) & arm_type %in% c("active_non_psy_placebo","actif_non_psy","active_non_psy") ~ "active_non_psy_placebo",
        TRUE ~ arm_type
      ),
      # infer when missing
      arm_type = ifelse(
        is.na(arm_type) & !is.na(arm_id) & str_detect(tolower(arm_id), "placebo"),
        ifelse(is.finite(dose_mg) & dose_mg > 0, "active_placebo", "inactive_placebo"),
        arm_type
      ),
      arm_type = ifelse(is.na(arm_type), "active", arm_type)
    )
  
  # 3) make group label nice for dose tiers (mild/low/high) if present in arm_id
  df <- df %>%
    mutate(
      group = ifelse(
        !is.na(arm_id) &
          (str_detect(tolower(arm_id), "mild") |
             str_detect(tolower(arm_id), "low")  |
             str_detect(tolower(arm_id), "high")),
        arm_id, group
      )
    )
  
  # Final tidy filter
  df <- df %>%
    filter(
      !is.na(study_id), !is.na(group), !is.na(molecule), !is.na(ae_term),
      is.finite(dose_mg), is.finite(events), is.finite(n)
    )
  
  df
}

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

# Utilitaire interne: convertir events en comptes si ce sont des proportions
#. Règle: si max(events) <= 1 et n >= 1, alors events := round(events * n)
.norm_events_to_counts <- function(df){
  if (!("events" %in% names(df)) || !("n" %in% names(df))) return(df)
  has_prop <- suppressWarnings(max(df$events, na.rm = TRUE) <= 1)
  if (is.finite(has_prop) && isTRUE(has_prop)) {
    df <- df %>% mutate(events = round(as.numeric(events) * as.numeric(n)))
  }
  df
}

# Choix de la référence dans un sous-ensemble (une étude x molécule)
.pick_reference_row <- function(d_sub, ref_prefs){
  # priorité: types de placebo dans l'ordre
  for (typ in ref_prefs){
    ref <- d_sub %>% filter(tolower(arm_type) == typ) %>% slice_head(n=1)
    if (nrow(ref)) return(ref)
  }
  # fallback: plus faible dose
  d_sub %>% arrange(dose_mg) %>% slice_head(n = 1)
}

suppressPackageStartupMessages({library(dplyr)})

# Fonction pour construire des contrastes 2x2 (ref vs actif)
# ref_policies = liste de priorité de références par molécule
build_pairwise_2x2 <- function(raw, ref_policies){
  stopifnot(all(c("study_id","molecule","ae_term","time_window","group","arm_type","dose_mg","events","n") %in% names(raw)))
  
  raw <- .norm_events_to_counts(raw)
  
  out <- list(); ii <- 1L
  
  for (sid in unique(raw$study_id)){
    d1 <- raw %>% filter(study_id == sid)
    
    for (mol in unique(d1$molecule)){
      d <- d1 %>% filter(molecule == mol)
      if (!nrow(d)) next
      
      ref_pref <- if (!is.null(ref_policies[[mol]])) ref_policies[[mol]] else ref_policies$.default
      ref_row  <- .pick_reference_row(d, ref_pref)
      
      # on va travailler par AE x time_window
      splits <- d %>% group_by(ae_term, time_window) %>% group_split()
      
      for (g in splits){
        # référence correspondante (même AE/window si possible)
        ref_g <- ref_row %>%
          filter(ae_term %in% unique(g$ae_term), time_window %in% unique(g$time_window))
        if (!nrow(ref_g)) ref_g <- ref_row
        
        # bras actifs = tout sauf la ligne de référence exacte (même group & même dose)
        act <- g %>%
          filter(!(group == ref_g$group[[1]] & abs(dose_mg - ref_g$dose_mg[[1]]) < 1e-12))
        if (!nrow(act)) next
        
        # construire une ligne par bras actif
        for (k in seq_len(nrow(act))){
          a <- act[k,]
          # 2×2
          ai <- as.numeric(a$events)
          bi <- as.numeric(a$n) - ai
          ci <- as.numeric(ref_g$events[[1]])
          di <- as.numeric(ref_g$n[[1]]) - ci
          
          out[[ii]] <- tibble::tibble(
            study_id    = a$study_id,
            molecule    = a$molecule,
            ae_term     = a$ae_term,
            time_window = a$time_window,
            # info bras de référence
            ref_group    = ref_g$group[[1]],
            ref_arm_type = ref_g$arm_type[[1]],
            ref_dose_mg  = as.numeric(ref_g$dose_mg[[1]]),
            # info bras actif
            cmp_group   = a$group,
            dose_mg     = as.numeric(a$dose_mg),
            dose_diff   = as.numeric(a$dose_mg) - as.numeric(ref_g$dose_mg[[1]]),
            # cellules 2×2 pour metafor::escalc
            ai = ai, bi = bi, ci = ci, di = di
          )
          ii <- ii + 1L
        }
      }
    }
  }
  
  dplyr::bind_rows(out)
}
=======
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

# ---- 1) Load & validate ----
load_data <- function(path = "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/data/Adverse-events-dose-v4_dose.xlsx",
                      sheet = "Feuil2") {
  stopifnot(file.exists(path))
  df <- readxl::read_excel(path, sheet = sheet) |>
    janitor::clean_names()
  df <- standardize_columns(df)
  validate_columns(df)
  df <- normalize_time_window(df)
  df
}

standardize_columns <- function(df){
  # Try to normalize common column names to a standard schema:
  rename_candidates <- c(
    "studyid" = "study_id",
    "study" = "study_id",
    "molecule_name" = "molecule",
    "drug" = "molecule",
    "dose" = "dose_mg",
    "dose_mg_total" = "dose_mg",
    "arm" = "arm_id",
    "armid" = "arm_id",
    "group" = "group",
    "condition" = "group",
    "timewindow" = "time_window",
    "window" = "time_window",
    "ae" = "ae_term",
    "adverse_event" = "ae_term",
    "events" = "events",
    "n_events" = "events",
    "total" = "total",
    "n_total" = "total",
    "n_participants_arm" = "total"
  )
  
  for (from in names(rename_candidates)) {
    to <- rename_candidates[[from]]
    if (!to %in% names(df) && from %in% names(df)) {
      df <- dplyr::rename(df, !!to := dplyr::all_of(from))
    }
  }
  df
}

validate_columns <- function(df){
  required <- c("study_id","molecule","dose_mg","time_window",
                "ae_term","arm_id","events","total")
  miss <- setdiff(required, names(df))
  if (length(miss)) {
    stop("Missing required columns: ", paste(miss, collapse=", "))
  }
  # Basic types
  df <- df |>
    mutate(
      dose_mg = suppressWarnings(as.numeric(dose_mg)),
      events  = suppressWarnings(as.integer(events)),
      total   = suppressWarnings(as.integer(total))
    )
  if (anyNA(df$dose_mg)) message("⚠️ Some dose_mg are NA (will treat as 0 for placebo if group marks it).")
  if (any(df$events > df$total, na.rm=TRUE)) {
    stop("Found rows with events > total.")
  }
  df
}

normalize_time_window <- function(df){
  df |>
    mutate(time_window = tolower(trimws(as.character(time_window)))) |>
    mutate(time_window = dplyr::case_when(
      time_window %in% c("session","acute","during","on-session") ~ "session",
      time_window %in% c("followup","follow-up","post","after") ~ "followup",
      TRUE ~ time_window
    ))
}

# ---- 2) Build pairwise contrasts within study × AE × window ----
# We try to identify a reference arm:
#   - Prefer explicit placebo if present (group == "placebo" or dose_mg == 0)
#   - Otherwise use the lowest dose in that study block
#
# Returns tidy effect-size-ready table (ai, bi, ci, di)
build_pairwise_2x2 <- function(df){
  df <- df |>
    mutate(
      group = tolower(as.character(dplyr::coalesce(group, ""))),
      is_placebo = group %in% c("placebo","control","inactive") | dose_mg == 0 | is.na(dose_mg)
    )
  
  by_keys <- c("study_id","molecule","ae_term","time_window")
  
  blocks <- df |>
    arrange(study_id, ae_term, time_window, dose_mg, arm_id) |>
    group_split(across(all_of(by_keys)), .keep = TRUE)
  
  make_contrasts <- function(g){
    # choose reference row
    ref <- if (any(g$is_placebo, na.rm=TRUE)) {
      g |> dplyr::filter(is_placebo) |> slice(1)
    } else {
      g |> arrange(dplyr::coalesce(dose_mg, Inf)) |> slice(1)
    }
    
    act <- g |> anti_join(ref |> select(arm_id), by = "arm_id")
    if (nrow(act) == 0L) return(tibble())
    
    tibble(
      study_id   = act$study_id,
      molecule   = act$molecule,
      ae_term    = act$ae_term,
      time_window= act$time_window,
      ref_arm_id = ref$arm_id[[1]],
      cmp_arm_id = act$arm_id,
      ref_dose_mg= ref$dose_mg[[1]],
      dose_mg    = act$dose_mg,
      ai = act$events,               # events in active
      bi = act$total - act$events,   # non-events in active
      ci = ref$events[[1]],          # events in control
      di = ref$total[[1]] - ref$events[[1]]
    )
  }
  
  purrr::map_dfr(blocks, make_contrasts)
}
>>>>>>> Stashed changes
