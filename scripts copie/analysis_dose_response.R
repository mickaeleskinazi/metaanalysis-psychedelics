suppressPackageStartupMessages({
<<<<<<< Updated upstream
  library(dplyr); library(purrr); library(tibble); library(metafor); library(splines)
})

# --- Escalc en conservant dose_mg, ref_dose_mg, dose_diff ---
build_escalc <- function(contrasts, add = 0.5){
  if (is.null(contrasts) || !nrow(contrasts)) return(tibble())
=======
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(metafor)
  library(splines)   # for ns()
})

# ---- 1) Convert 2x2 to log-OR yi/vi via metafor::escalc ----
build_escalc <- function(contrasts, add = 0.5){
  if (!nrow(contrasts)) return(tibble())
>>>>>>> Stashed changes
  esc <- metafor::escalc(
    measure = "OR",
    ai = ai, bi = bi, ci = ci, di = di,
    add = add, to = "only0",
    data = contrasts
  )
<<<<<<< Updated upstream
  esc %>%
    transmute(
      study_id, molecule, ae_term, time_window,
      dose_mg = as.numeric(dose_mg),
      ref_dose_mg = as.numeric(ref_dose_mg),
      dose_diff = as.numeric(dose_diff),
=======
  esc |>
    transmute(
      study_id, molecule, ae_term, time_window,
      dose_mg,
>>>>>>> Stashed changes
      yi, vi
    )
}

<<<<<<< Updated upstream
# --- Helpers ---
make_dose_grid <- function(x, n = 120){
  x <- x[is.finite(x)]
  if (!length(x)) return(numeric(0))
  r <- range(x)
  if (r[1] == r[2]) return(sort(unique(x)))
  seq(r[1], r[2], length.out = n)
}

fit_dr_linear <- function(dat){
  dat <- dat %>% filter(is.finite(yi), is.finite(vi), is.finite(dose_diff))
  if (nrow(dat) < 2 || dplyr::n_distinct(dat$dose_diff) < 2) return(NULL)
  tryCatch(metafor::rma(yi, vi, mods = ~ dose_diff, data = dat, method = "REML"),
           error = function(e) NULL)
}

fit_dr_spline <- function(dat, df_spline = 3){
  dat <- dat %>% filter(is.finite(yi), is.finite(vi), is.finite(dose_diff))
  if (dplyr::n_distinct(dat$dose_diff) < 3 || nrow(dat) < df_spline + 1) return(NULL)
  tryCatch(metafor::rma(yi, vi, mods = ~ ns(dose_diff, df = df_spline),
                        data = dat, method = "REML"),
           error = function(e) NULL)
=======
# ---- 2) Fit dose–response models (linear + optional spline) ----
# For stability, we fit per (molecule × time_window) and optionally per AE.
fit_dr_linear <- function(dat){
  # dat columns: yi, vi, dose_mg, study_id
  if (nrow(dat) < 2) return(NULL)
  # Ensure finite & variation in dose
  dat <- dat |> filter(is.finite(yi), is.finite(vi), is.finite(dose_mg))
  if (n_distinct(dat$dose_mg) < 2) return(NULL)
  
  tryCatch({
    m <- rma(yi, vi, mods = ~ dose_mg, data = dat, method = "REML")
    m
  }, error = function(e) NULL)
}

fit_dr_spline <- function(dat, df_spline = 3){
  dat <- dat |> filter(is.finite(yi), is.finite(vi), is.finite(dose_mg))
  if (nrow(dat) < df_spline + 1 || n_distinct(dat$dose_mg) < df_spline) return(NULL)
  
  tryCatch({
    m <- rma(yi, vi, mods = ~ ns(dose_mg, df = df_spline), data = dat, method = "REML")
    m
  }, error = function(e) NULL)
>>>>>>> Stashed changes
}

tidy_rma <- function(m, model_label){
  if (is.null(m)) return(tibble())
<<<<<<< Updated upstream
  coefs <- tryCatch(coef(m), error = function(e) numeric(0))
  if (length(coefs) == 0) return(tibble())
  V   <- tryCatch(vcov(m), error = function(e) NULL)
  ses <- if (!is.null(V)) sqrt(diag(V)) else rep(NA_real_, length(coefs))
  z   <- coefs / ses
  p   <- 2*pnorm(abs(z), lower.tail = FALSE)
  ci  <- tryCatch(confint(m), error = function(e) NULL)
  ci_low  <- if (!is.null(ci) && "ci.lb" %in% names(ci)) as.numeric(ci$ci.lb) else rep(NA_real_, length(coefs))
  ci_high <- if (!is.null(ci) && "ci.ub" %in% names(ci)) as.numeric(ci$ci.ub) else rep(NA_real_, length(coefs))
  QM  <- suppressWarnings(as.numeric(m$QM)); QMp <- suppressWarnings(as.numeric(m$QMp))
=======
  coefs <- coef(m)
  ses   <- sqrt(diag(vcov(m)))
  z     <- coefs / ses
  p     <- 2*pnorm(abs(z), lower.tail = FALSE)
  ci    <- confint(m)
  
>>>>>>> Stashed changes
  tibble(
    model   = model_label,
    term    = names(coefs),
    estimate= as.numeric(coefs),
    se      = as.numeric(ses),
    z       = as.numeric(z),
    pval    = as.numeric(p),
<<<<<<< Updated upstream
    ci_low  = ci_low,
    ci_high = ci_high,
    tau2    = suppressWarnings(as.numeric(m$tau2)),
    QE      = suppressWarnings(as.numeric(m$QE)),
    k       = suppressWarnings(as.integer(m$k)),
    I2      = suppressWarnings(pmax(0, (m$QE - (m$k - 1))/m$QE) * 100),
    QM      = QM, QMp = QMp
  )
}

predict_spline <- function(m, doses_diff){
  if (is.null(m) || !length(doses_diff)) return(tibble())
  X_fit <- tryCatch(m$X, error = function(e) NULL); if (is.null(X_fit)) return(tibble())
  p <- ncol(X_fit) - 1L; if (p <= 0) return(tibble())
  B <- splines::ns(doses_diff, df = p) %>% as.matrix()
  colnames(B) <- colnames(X_fit)[-1]
  pr <- predict(m, newmods = B)
  tibble(dose_diff = doses_diff, fit = as.numeric(pr$pred), lwr = as.numeric(pr$ci.lb), upr = as.numeric(pr$ci.ub))
}

predict_linear <- function(m, doses_diff){
  if (is.null(m) || !length(doses_diff)) return(tibble())
  pr <- predict(m, newmods = as.numeric(doses_diff))
  tibble(dose_diff = doses_diff, fit = as.numeric(pr$pred), lwr = as.numeric(pr$ci.lb), upr = as.numeric(pr$ci.ub))
}

# --- DR globale par molécule ---
run_dr_by_molecule <- function(es, min_k = 4, fit_spline = TRUE){
  if (!nrow(es)) return(list(models = tibble(), preds = tibble()))
  es2 <- es %>%
    group_by(molecule) %>%
    filter(n() >= min_k, n_distinct(dose_diff[is.finite(dose_diff)]) >= 2) %>%
    ungroup()
  groups <- es2 %>% group_by(molecule) %>% group_split()
  
  models_out <- purrr::map_dfr(groups, function(g){
    m_lin <- fit_dr_linear(g); t_lin <- tidy_rma(m_lin, "linear")
    m_spl <- if (fit_spline) fit_dr_spline(g) else NULL
    t_spl <- tidy_rma(m_spl, "spline_df3")
    bind_rows(t_lin, t_spl) %>% mutate(molecule = g$molecule[[1]], .before = 1)
  })
  
  preds_out <- purrr::map_dfr(groups, function(g){
    doses <- make_dose_grid(g$dose_diff, n = 120)
    if (!length(doses)) return(tibble())
    m_lin <- fit_dr_linear(g); m_spl <- if (fit_spline) fit_dr_spline(g) else NULL
    p <- if (!is.null(m_spl)) predict_spline(m_spl, doses) %>% mutate(model = "spline_df3")
    else if (!is.null(m_lin)) predict_linear(m_lin, doses) %>% mutate(model = "linear")
    else tibble()
    if (!nrow(p)) return(p)
    ref_mean <- mean(g$ref_dose_mg, na.rm = TRUE)
    p %>% mutate(molecule = g$molecule[[1]], dose_mg = dose_diff + ref_mean, .before = 1)
  })
  
  list(models = models_out, preds = preds_out)
}

# --- DR par AE x molécule ---
run_dr_by_ae <- function(es, min_k = 4, fit_spline = TRUE){
  if (!nrow(es)) return(list(models = tibble(), preds = tibble()))
  es2 <- es %>%
    group_by(ae_term, molecule) %>%
    filter(n() >= min_k, n_distinct(dose_diff[is.finite(dose_diff)]) >= 2) %>%
    ungroup()
  groups <- es2 %>% group_by(ae_term, molecule) %>% group_split()
  
  models_out <- purrr::map_dfr(groups, function(g){
    m_lin <- fit_dr_linear(g); t_lin <- tidy_rma(m_lin, "linear")
    m_spl <- if (fit_spline) fit_dr_spline(g) else NULL
    t_spl <- tidy_rma(m_spl, "spline_df3")
    bind_rows(t_lin, t_spl) %>%
      mutate(ae_term = g$ae_term[[1]], molecule = g$molecule[[1]], .before = 1)
  })
  
  preds_out <- purrr::map_dfr(groups, function(g){
    doses <- make_dose_grid(g$dose_diff, n = 120)
    if (!length(doses)) return(tibble())
    m_lin <- fit_dr_linear(g); m_spl <- if (fit_spline) fit_dr_spline(g) else NULL
    p <- if (!is.null(m_spl)) predict_spline(m_spl, doses) %>% mutate(model = "spline_df3")
    else if (!is.null(m_lin)) predict_linear(m_lin, doses) %>% mutate(model = "linear")
    else tibble()
    if (!nrow(p)) return(p)
    ref_mean <- mean(g$ref_dose_mg, na.rm = TRUE)
    p %>% mutate(ae_term = g$ae_term[[1]], molecule = g$molecule[[1]], dose_mg = dose_diff + ref_mean, .before = 1)
  })
  
  list(models = models_out, preds = preds_out)
}

# ---- LSD dose-response: linear vs spline (robust to dropped spline cols) ----
run_dr_lsd_compare <- function(es){
  library(metafor); library(dplyr)
  dat <- es %>% filter(molecule == "LSD",
                       is.finite(yi), is.finite(vi),
                       is.finite(dose_diff))
  if (!nrow(dat)) return(list(linear=NULL, spline=NULL))
  
  mod_lin <- tryCatch(
    rma(yi, vi, mods = ~ dose_diff, data = dat, method="REML"),
    error = function(e) NULL
  )
  
  mod_spl <- tryCatch(
    rma(yi, vi, mods = ~ splines::ns(dose_diff, df=3), data = dat, method="REML"),
    error = function(e) NULL
  )
  
  list(linear = mod_lin, spline = mod_spl, data = dat)
}

.predict_with_same_basis <- function(model, new_dose_diff){
  # Build a newmods matrix whose column names exactly match model$X (minus intercept)
  X <- tryCatch(model$X, error = function(e) NULL)
  if (is.null(X)) return(NULL)
  col_target <- colnames(X)
  # first column is intercept "(Intrcpt)" in metafor; drop it
  col_target <- col_target[setdiff(seq_along(col_target), 1L)]
  p <- length(col_target)
  if (p == 0) {
    # no moderators; predict with NULL newmods
    return(predict(model, newmods = NULL))
  }
  # Build a ns() with df = p to match kept columns, then force names to match
  B <- as.matrix(splines::ns(new_dose_diff, df = p))
  # Ensure column count matches; if needed, truncate or recycle safely
  if (ncol(B) != p) {
    # Align shape
    if (ncol(B) > p) B <- B[, seq_len(p), drop = FALSE]
    if (ncol(B) < p) {
      # pad with zeros (conservative)
      B <- cbind(B, matrix(0, nrow = nrow(B), ncol = p - ncol(B)))
    }
  }
  colnames(B) <- col_target
  predict(model, newmods = B)
}

plot_dr_lsd_compare <- function(models){
  library(ggplot2); library(dplyr)
  dat <- models$data
  if (is.null(dat) || !nrow(dat)) return(NULL)
  
  # grid over the observed range
  grid <- seq(min(dat$dose_diff, na.rm=TRUE),
              max(dat$dose_diff, na.rm=TRUE),
              length.out = 120)
  
  preds <- list()
  
  if (!is.null(models$linear)) {
    p <- predict(models$linear, newmods = grid)
    preds$linear <- tibble(
      model = "linear",
      dose_diff = grid,
      fit = as.numeric(p$pred),
      lwr = as.numeric(p$ci.lb),
      upr = as.numeric(p$ci.ub)
    )
  }
  
  if (!is.null(models$spline)) {
    p <- .predict_with_same_basis(models$spline, grid)
    if (!inherits(p, "try-error") && !is.null(p)) {
      preds$spline <- tibble(
        model = "spline",
        dose_diff = grid,
        fit = as.numeric(p$pred),
        lwr = as.numeric(p$ci.lb),
        upr = as.numeric(p$ci.ub)
      )
    }
  }
  
  dfp <- bind_rows(preds)
  if (!nrow(dfp)) return(NULL)
  
  ggplot(dfp, aes(dose_diff, fit, color = model, fill = model)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.18, color = NA) +
    labs(
      x = "LSD dose difference (mg vs reference)",
      y = "log(OR) vs reference",
      title = "LSD dose–response: linear vs spline (robust prediction)"
    ) +
    theme_bw()
}

# (optional) print linear coef/p-value for quick interpretation
print_lsd_linear_summary <- function(models){
  if (is.null(models$linear)) { message("No linear model fitted."); return(invisible(NULL)) }
  s <- summary(models$linear)
  cat("\nLSD linear meta-regression (logOR ~ dose_diff)\n",
      "Intercept:", round(s$beta[1,1], 3), 
      " [", round(s$ci.lb[1],3), ", ", round(s$ci.ub[1],3), "]\n",
      "Slope (per mg):", round(s$beta[2,1], 3),
      " [", round(s$ci.lb[2],3), ", ", round(s$ci.ub[2],3), "]",
      "  p=", signif(s$pval[2], 3), "\n", sep = "")
  invisible(s)
=======
    ci_low  = as.numeric(ci$ci.lb),
    ci_high = as.numeric(ci$ci.ub),
    tau2    = as.numeric(m$tau2),
    QE      = as.numeric(m$QE),
    k       = as.integer(m$k),
    I2      = pmax(0, (m$QE - (m$k - 1))/m$QE) * 100
  )
}

# ---- 3) High-level runners ----
run_dose_response <- function(es, by_ae = FALSE, fit_spline = TRUE){
  if (!nrow(es)) return(tibble())
  
  group_keys <- c("molecule","time_window")
  if (by_ae) group_keys <- c(group_keys, "ae_term")
  
  blocks <- es |> group_split(across(all_of(group_keys)), .keep = TRUE)
  
  results <- purrr::map_dfr(blocks, function(g){
    meta <- g |> select(any_of(group_keys)) |> slice(1)
    # Linear
    m_lin <- fit_dr_linear(g |> select(yi, vi, dose_mg, study_id))
    t_lin <- tidy_rma(m_lin, "linear")
    # Spline
    t_spl <- tibble()
    if (fit_spline) {
      m_spl <- fit_dr_spline(g |> select(yi, vi, dose_mg, study_id))
      t_spl <- tidy_rma(m_spl, "spline_df3")
    }
    bind_rows(t_lin, t_spl) |>
      mutate(
        molecule   = meta$molecule[[1]],
        time_window= meta$time_window[[1]],
        ae_term    = if ("ae_term" %in% names(meta)) meta$ae_term[[1]] else NA_character_
      ) |>
      relocate(molecule, time_window, ae_term, .before = model)
  })
  
  results
>>>>>>> Stashed changes
}
