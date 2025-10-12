suppressPackageStartupMessages({
  library(dplyr); library(janitor); library(rlang); library(stringr)
})

.map_groups_dfr <- function(.data, .f, ..., .drop_keys = TRUE, .template = NULL) {
  if (!dplyr::is_grouped_df(.data)) {
    stop("`.map_groups_dfr()` requires a grouped data frame.")
  }

  fn <- rlang::as_function(.f)

  pieces <- dplyr::group_map(
    .data,
    function(.x, .y) {
      out <- fn(.x, .y, ...)

      if (is.null(out)) {
        out <- dplyr::tibble()
      } else {
        out <- dplyr::as_tibble(out)
      }

      if (!is.null(.template) && !nrow(out)) {
        out <- .template[0, , drop = FALSE]
      }

      if (.drop_keys) {
        keep_cols <- setdiff(names(out), names(.y))
        out <- out[, keep_cols, drop = FALSE]

        if (!nrow(out)) {
          return(out)
        }

        key_df <- .y[rep(1, nrow(out)), , drop = FALSE]
        return(dplyr::bind_cols(key_df, out))
      }

      if (!nrow(out)) {
        return(out)
      }

      missing_keys <- setdiff(names(.y), names(out))
      if (length(missing_keys)) {
        key_df <- .y[rep(1, nrow(out)), missing_keys, drop = FALSE]
        out <- dplyr::bind_cols(key_df, out)
      }

      out
    },
    .keep = TRUE
  )

  if (!length(pieces)) {
    if (!is.null(.template)) {
      return(.template[0, , drop = FALSE])
    }
    return(dplyr::tibble())
  }

  res <- dplyr::bind_rows(pieces)
  if (!nrow(res) && !is.null(.template)) {
    return(.template[0, , drop = FALSE])
  }

  res
}

.norm_key <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x[is.na(x)] <- ""
  x <- tolower(x)
  gsub("[^a-z0-9]+", "", x)
}

.detect_and_rename <- function(df, targets){
  cn <- names(df)
  cn_norm <- .norm_key(cn)
  getcol <- function(cands){
    cands_norm <- .norm_key(cands)
    # exact first, then contains (using normalized tokens for accent-insensitive match)
    for (cand in cands_norm){
      ix <- which(cn_norm == cand)
      if (length(ix)) return(cn[ix[1]])
    }
    for (cand in cands_norm){
      hit <- which(grepl(cand, cn_norm, fixed = TRUE))
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
    ae_term     = c(
      "ae_term","ae_terms","adverse_event","ae","outcome","event_name","event_term",
      "event_label","ae_name","ae_label","ae_description","adverse_event_term",
      "adverse_event_name","adverse_event_label","adverse_effect","side_effect",
      "symptom","reaction","effet_indesirable","effets_indesirables",
      "effet_secondaire","effets_secondaires","se_adverse_event"
    ),
    time_window = c("time_window","window","timepoint","time","visit"),
    dose_mg     = c("dose_mg","dose","dose_mg_numeric","dose_mg_num","dose_milligram","lsd_dose_mg"),
    events      = c(
      "events","event","ae_n","cases","num_events","n_events",
      "events_arm","ae_count","number_of_events","count_events"
    ),
    n           = c(
      "n","total","n_total","sample_size","denominator",
      "participants","participants_total","n_participants",
      "n_participants_arm","arm_n","n_arm","participants_arm"
    )
  )
  df <- .detect_and_rename(df, targets)

  # Heuristic fallback for AE labels: when none of the aliases matched above we
  # try to infer a reasonable candidate based on keyword combinations while
  # avoiding the "events" column. This mirrors the intuition analysts use when
  # manually inspecting unfamiliar spreadsheets.
  if (!("ae_term" %in% names(df))) {
    nm_lower <- tolower(iconv(names(df), to = "ASCII//TRANSLIT"))
    keyword_sets <- list(
      c("ae", "term"),
      c("ae", "label"),
      c("adverse", "event"),
      c("event", "term"),
      c("event", "label"),
      c("event", "desc"),
      c("event", "name"),
      c("side", "effect"),
      c("effet", "ind"),
      c("effet", "second"),
      "symptom"
    )
    pick_keywords <- function(keywords){
      if (length(keywords) == 1L) {
        hits <- which(grepl(keywords, nm_lower, fixed = TRUE))
      } else {
        hits <- which(vapply(nm_lower, function(nm){
          all(vapply(keywords, function(kw) grepl(kw, nm, fixed = TRUE), logical(1)))
        }, logical(1)))
      }
      # filter out the events column which is frequently just "events"
      hits[names(df)[hits] != "events"]
    }
    idx <- NULL
    for (kw in keyword_sets){
      hits <- pick_keywords(kw)
      if (length(hits)) {
        idx <- hits[[1]]
        break
      }
    }
    if (!is.null(idx)) {
      picked <- names(df)[[idx]]
      rlang::inform(paste0("Auto-detected AE label column '", picked, "'"))
      df <- dplyr::rename(df, ae_term = !!rlang::sym(picked))
    }
  }

  if (!("ae_term" %in% names(df))) {
    exclude <- c(
      "study_id", "author_year", "arm_id", "group", "arm_type", "molecule",
      "time_window", "dose_mg", "events", "n", "n_total", "n_events"
    )
    candidates <- setdiff(names(df), exclude)
    if (length(candidates)) {
      score_col <- function(col){
        nm <- .norm_key(col)
        score <- 0
        pattern_hits <- c(
          "aeterm", "aelabel", "adverseevent", "sideeffect", "symptom",
          "effetindesirable", "effetsindesirables", "effetsecondaire",
          "effetssecondaires", "reaction", "safety", "teae"
        )
        for (pat in pattern_hits) {
          if (grepl(pat, nm, fixed = TRUE)) score <- score + 2
        }
        vals <- df[[col]]
        if (is.character(vals) || is.factor(vals)) {
          val_norm <- unique(.norm_key(trimws(as.character(vals))))
          val_norm <- val_norm[nchar(val_norm) > 0]
          if (length(val_norm)) {
            if (any(grepl("event", val_norm, fixed = TRUE))) score <- score + 1
            if (any(grepl("effet",  val_norm, fixed = TRUE))) score <- score + 1
            if (any(grepl("symptom", val_norm, fixed = TRUE))) score <- score + 1
            if (length(val_norm) > 1) score <- score + 1
          }
        }
        score
      }
      scores <- vapply(candidates, score_col, numeric(1))
      if (length(scores) && max(scores) >= 2) {
        picked <- candidates[[which.max(scores)]]
        rlang::inform(
          paste0(
            "Heuristically picked column '", picked,
            "' to serve as adverse event labels based on its contents."
          )
        )
        df <- dplyr::rename(df, ae_term = !!rlang::sym(picked))
      }
    }
  }

  # Some legacy ingestion scripts used ``n_total`` / ``n_events`` for counts.
  # Harmonise them here if they slipped through the automatic detection above.
  if ("n_total" %in% names(df) && !("n" %in% names(df))) {
    df <- dplyr::rename(df, n = n_total)
  }
  if ("n_events" %in% names(df) && !("events" %in% names(df))) {
    df <- dplyr::rename(df, events = n_events)
  }
  
  # Friendly error if must-haves are missing (except group; we can synthesize it)
  must_have <- c("study_id","molecule","ae_term","time_window","dose_mg","events","n")
  missing <- setdiff(must_have, names(df))
  if (length(missing)){
    rlang::inform(paste0("Columns detected in sheet: ", paste(names(df), collapse = ", ")))
    rlang::abort(
      message = paste0(
        "Missing required columns after auto-detect: ",
        paste(missing, collapse = ", ")
      )
    )
  }
  
  # Normalize basic types
  df <- df %>%
    mutate(
      study_id    = as.character(study_id),
      molecule    = toupper(as.character(molecule)),
      ae_term     = dplyr::na_if(trimws(as.character(ae_term)), ""),
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
  needed <- c("study_id","molecule","ae_term","time_window","group",
              "arm_type","dose_mg","events","n")
  missing_cols <- setdiff(needed, names(raw))
  if (length(missing_cols)) {
    rlang::abort(
      message = paste0(
        "build_pairwise_2x2() is missing required columns: ",
        paste(missing_cols, collapse = ", "),
        ". Available columns: ",
        paste(names(raw), collapse = ", ")
      )
    )
  }

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