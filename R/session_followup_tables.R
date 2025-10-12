suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(metafor)
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

.pluck_term_or_tail <- function(x, term) {
  if (is.null(x)) return(NA_real_)
  nm <- names(x)
  if (!is.null(nm) && length(nm)) {
    idx <- which(nm == term)
    if (length(idx) == 0) {
      idx <- which(tolower(nm) == tolower(term))
    }
    if (length(idx) >= 1) {
      return(as.numeric(x[[idx[[1]]]]))
    }
  }
  len <- length(x)
  if (len == 0) return(NA_real_)
  if (len == 1) return(as.numeric(x[[1]]))
  as.numeric(x[[len]])
}

norm_window_label <- function(x) {
  x |>
    as.character() |>
    tolower() |>
    gsub("[^a-z0-9]+", "_", x = _) |>
    gsub("_+", "_", x = _) |>
    sub("^_", "", x = _) |>
    sub("_$", "", x = _)
}

build_es_for_window <- function(raw_all, window_value, ref_policies) {
  stopifnot("time_window" %in% names(raw_all))
  raw_all <- raw_all %>% mutate(time_window = norm_window_label(time_window))
  if (!window_value %in% unique(raw_all$time_window)) {
    stop(
      "Window '", window_value, "' not found. Available: ",
      paste(sort(unique(raw_all$time_window)), collapse = ", ")
    )
  }
  raw_win <- raw_all %>% filter(time_window == window_value)
  contr   <- build_pairwise_2x2(raw_win, ref_policies = ref_policies)
  es      <- build_escalc(contr)
  es$time_window <- window_value
  es
}

sig_stars <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

dr_fit_per_window <- function(es, min_k_per_window = 2) {
  stopifnot(all(c("molecule", "dose_mg", "time_window", "yi", "vi") %in% names(es)))
  grouped <- es %>%
    filter(!is.na(dose_mg), is.finite(yi), is.finite(vi)) %>%
    group_by(molecule, time_window)

  .map_groups_dfr(grouped, function(dat, key) {
    empty <- tibble(
      k              = integer(),
      beta_dose      = double(),
      se_dose        = double(),
      z_dose         = double(),
      p_dose         = double(),
      ci_lb          = double(),
      ci_ub          = double(),
      I2             = double(),
      tau2           = double()
    )
    if (nrow(dat) < min_k_per_window) return(empty)
    m <- tryCatch(rma(yi ~ dose_mg, vi = vi, data = dat, method = "REML"), error = function(e) NULL)
    if (is.null(m)) return(empty)
    tibble(
      k              = m$k,
      beta_dose      = as.numeric(coef(m)["dose_mg"]),
      se_dose        = as.numeric(m$se["dose_mg"]),
      z_dose         = as.numeric(m$zval["dose_mg"]),
      p_dose         = as.numeric(m$pval["dose_mg"]),
      ci_lb          = .pluck_term_or_tail(m$ci.lb, "dose_mg"),
      ci_ub          = .pluck_term_or_tail(m$ci.ub, "dose_mg"),
      I2             = suppressWarnings(tryCatch(m$I2,  error = function(e) NA_real_)),
      tau2           = suppressWarnings(tryCatch(m$tau2, error = function(e) NA_real_))
    )
  }) %>%
    ungroup() %>%
    mutate(stars = sig_stars(p_dose))
}

dr_test_session_vs_followup <- function(es, min_k_total = 4) {
  stopifnot(all(c("molecule", "dose_mg", "time_window", "yi", "vi") %in% names(es)))
  grouped <- es %>%
    filter(!is.na(dose_mg), is.finite(yi), is.finite(vi)) %>%
    group_by(molecule)

  .map_groups_dfr(grouped, function(dat, key) {
    empty <- tibble(
      beta_dose_main    = double(),
      beta_interaction  = double(),
      p_interaction     = double()
    )
    if (length(unique(dat$time_window)) < 2 || nrow(dat) < min_k_total) return(empty)
    m <- tryCatch(rma(yi ~ dose_mg * time_window, vi = vi, data = dat, method = "REML"), error = function(e) NULL)
    if (is.null(m)) return(empty)
    tibble(
      beta_dose_main    = as.numeric(coef(m)["dose_mg"]),
      beta_interaction  = as.numeric(coef(m)[grep("^dose_mg:time_window", names(coef(m)))]),
      p_interaction     = as.numeric(m$pval[grep("^dose_mg:time_window", names(m$pval))])
    )
  }) %>%
    ungroup() %>%
    mutate(stars_interaction = sig_stars(p_interaction))
}

make_dr_window_tables <- function(
    es,
    out_dir,
    min_k_per_window = 2,
    min_k_total = 4
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  per_window <- dr_fit_per_window(es, min_k_per_window = min_k_per_window)
  write_csv(per_window, file.path(out_dir, "dr_per_window_by_molecule.csv"))

  per_window_wide <- per_window %>%
    select(molecule, time_window, beta_dose, ci_lb, ci_ub, p_dose, stars, k, I2, tau2) %>%
    mutate(time_window = case_when(
      time_window %in% c("session", "Session") ~ "session",
      time_window %in% c("follow_up", "follow-up", "Follow_up") ~ "follow_up",
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
    left_join(interact %>% select(molecule, beta_interaction, p_interaction, stars_interaction), by = "molecule") %>%
    mutate(
      session_effect = ifelse(
        !is.na(beta_dose.session),
        sprintf("%.3f [%.3f, %.3f]%s",
                beta_dose.session, ci_lb.session, ci_ub.session,
                ifelse(is.na(stars.session), "", paste0(" ", stars.session))),
        NA_character_
      ),
      followup_effect = ifelse(
        !is.na(beta_dose.follow_up),
        sprintf("%.3f [%.3f, %.3f]%s",
                beta_dose.follow_up, ci_lb.follow_up, ci_ub.follow_up,
                ifelse(is.na(stars.follow_up), "", paste0(" ", stars.follow_up))),
        NA_character_
      ),
      interaction = ifelse(
        !is.na(beta_interaction),
        sprintf("%.3f (p=%s)%s",
                beta_interaction,
                formatC(p_interaction, format = "f", digits = 3),
                ifelse(is.na(stars_interaction), "", paste0(" ", stars_interaction))),
        NA_character_
      )
    ) %>%
    select(
      molecule,
      session_effect,
      followup_effect,
      interaction
    )
  write_csv(pub_tab, file.path(out_dir, "dr_session_followup_publication_table.csv"))

  invisible(list(
    per_window = per_window,
    per_window_wide = per_window_wide,
    interaction = interact,
    publication = pub_tab
  ))
}
