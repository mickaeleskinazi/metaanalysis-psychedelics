# =========================
# FILE: R/dose_response_models.R
# =========================
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(metafor)
  library(splines)
})

# ============================================================
# Effect sizes (OR) with HA correction only for zero cells
# ============================================================
build_escalc <- function(contrasts, add = 0.5){
  empty_tpl <- tibble(
    study_id    = character(),
    molecule    = character(),
    ae_term     = character(),
    time_window = character(),
    dose_mg     = double(),
    ref_dose_mg = double(),
    dose_diff   = double(),
    yi          = double(),
    vi          = double()
  )
  if (is.null(contrasts) || !nrow(contrasts)) return(empty_tpl)
  
  needed <- c("ai","bi","ci","di","study_id","molecule","ae_term",
              "time_window","dose_mg","ref_dose_mg","dose_diff")
  miss <- setdiff(needed, names(contrasts))
  if (length(miss)) stop("build_escalc(): missing columns: ", paste(miss, collapse=", "))
  
  esc <- metafor::escalc(
    measure = "OR",
    ai = ai, bi = bi, ci = ci, di = di,
    data = contrasts,
    add = add, to = "only0"
  )
  
  as_tibble(esc) %>%
    transmute(
      study_id    = as.character(study_id),
      molecule    = toupper(as.character(molecule)),
      ae_term     = as.character(ae_term),
      time_window = as.character(time_window),
      dose_mg     = as.numeric(dose_mg),
      ref_dose_mg = as.numeric(ref_dose_mg),
      dose_diff   = as.numeric(dose_diff),
      yi          = as.numeric(yi),
      vi          = as.numeric(vi)
    ) %>%
    filter(is.finite(yi), is.finite(vi), is.finite(dose_mg))
}

# ============================================================
# Grid helper
# ============================================================
make_grid <- function(doses_obs, grid = c("observed","truncated","continuous"), n_grid = 120){
  grid <- match.arg(grid)
  doses_obs <- sort(unique(as.numeric(doses_obs)))
  doses_obs <- doses_obs[is.finite(doses_obs)]
  if (!length(doses_obs)) return(double())
  
  if (grid == "observed") return(doses_obs)
  
  if (grid == "truncated") {
    if (length(doses_obs) < 3) return(doses_obs)
    q <- stats::quantile(doses_obs, c(0.05, 0.95), na.rm = TRUE, names = FALSE)
    return(doses_obs[doses_obs >= q[1] & doses_obs <= q[2]])
  }
  
  rng <- range(doses_obs, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) return(doses_obs)
  seq(rng[1], rng[2], length.out = n_grid)
}

# ============================================================
# Fit DR models
# - linear: yi ~ dose_mg
# - spline: yi ~ ns(dose_mg, df)
#
# IMPORTANT FIX:
# metafor::predict.rma() matches newmods by *names*.
# We store the exact names of fitted mod columns and reuse them at prediction.
# ============================================================
fit_dr <- function(dat, model = c("linear","spline"), df_spline = 3){
  model <- match.arg(model)
  dat <- dat %>% filter(is.finite(yi), is.finite(vi), is.finite(dose_mg))
  if (nrow(dat) < 2 || n_distinct(dat$dose_mg) < 2) return(NULL)
  
  if (model == "linear") {
    m <- tryCatch(
      metafor::rma(yi, vi, mods = ~ dose_mg, data = dat, method = "REML"),
      error = function(e) NULL
    )
    return(m)
  }
  
  # spline needs at least 3 distinct doses to be meaningful/stable
  if (nrow(dat) < 3 || n_distinct(dat$dose_mg) < 3) return(NULL)
  
  X <- splines::ns(dat$dose_mg, df = df_spline)
  m <- tryCatch(
    metafor::rma(yi, vi, mods = X, data = dat, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)
  
  # Capture exact mod column names used internally by metafor
  mod_names <- tryCatch(colnames(m$X)[-1], error = function(e) NULL)  # drop intercept
  
  list(
    model     = m,
    df_spline = df_spline,
    knots     = attr(X, "knots"),
    bknots    = attr(X, "Boundary.knots"),
    mod_names = mod_names
  )
}

predict_dr <- function(m, doses, model = c("linear","spline")){
  model <- match.arg(model)
  if (is.null(m) || !length(doses)) return(tibble())
  
  if (model == "linear") {
    pr <- predict(m, newmods = doses)
    return(tibble(dose_mg = doses, fit = pr$pred, lwr = pr$ci.lb, upr = pr$ci.ub))
  }
  
  # spline: m is a list created by fit_dr()
  Xnew <- splines::ns(
    doses,
    df = m$df_spline,
    knots = m$knots,
    Boundary.knots = m$bknots
  )
  
  # Force names to match fitted model (prevents "Could not find variable 's3' in the model")
  if (!is.null(m$mod_names) && ncol(Xnew) == length(m$mod_names)) {
    colnames(Xnew) <- m$mod_names
  } else {
    colnames(Xnew) <- NULL  # fallback: positional matching
  }
  
  pr <- predict(m$model, newmods = Xnew)
  tibble(dose_mg = doses, fit = pr$pred, lwr = pr$ci.lb, upr = pr$ci.ub)
}

# ============================================================
# Global DR by molecule (linear + spline supported)
# ============================================================
run_dr_by_molecule <- function(es, min_k = 2,
                               model = c("linear","spline"),
                               df_spline = 3,
                               grid = c("observed","truncated","continuous"),
                               n_grid = 120){
  model <- match.arg(model)
  grid  <- match.arg(grid)
  
  es2 <- es %>%
    filter(is.finite(yi), is.finite(vi), is.finite(dose_mg)) %>%
    group_by(molecule) %>%
    filter(n() >= min_k, n_distinct(dose_mg) >= 2) %>%
    ungroup()
  
  preds <- map_dfr(split(es2, es2$molecule), function(g){
    doses <- make_grid(g$dose_mg, grid = grid, n_grid = n_grid)
    m <- fit_dr(g, model = model, df_spline = df_spline)
    if (is.null(m)) return(tibble())
    
    predict_dr(m, doses, model = model) %>%
      mutate(molecule = g$molecule[1], model = model)
  })
  
  models <- es2 %>%
    split(.$molecule) %>%
    map_dfr(function(g){
      m <- fit_dr(g, model = model, df_spline = df_spline)
      if (is.null(m)) return(tibble())
      
      mm <- if (model == "spline") m$model else m
      
      tibble(
        molecule  = g$molecule[1],
        model     = model,
        df_spline = ifelse(model == "spline", df_spline, NA_integer_),
        k         = mm$k,
        QM        = mm$QM,
        QMp       = mm$QMp,
        I2        = mm$I2,
        # linear: slope + p; spline: omnibus p = QMp
        beta      = if (model == "linear") as.numeric(coef(mm)[2]) else NA_real_,
        pval      = if (model == "linear") as.numeric(mm$pval[2]) else as.numeric(mm$QMp)
      )
    })
  
  list(preds = preds, models = models)
}

# ============================================================
# DR by AE × molecule
# (You can still run spline here; if you don't want spline for AEs,
# just call run_dr_by_ae(..., model="linear") in run_main_analysis)
# ============================================================
run_dr_by_ae <- function(es, min_k = 2,
                         model = c("linear","spline"),
                         df_spline = 3,
                         grid = c("observed","truncated","continuous"),
                         n_grid = 120){
  model <- match.arg(model)
  grid  <- match.arg(grid)
  
  es2 <- es %>%
    filter(is.finite(yi), is.finite(vi), is.finite(dose_mg)) %>%
    group_by(molecule, ae_term) %>%
    filter(n() >= min_k, n_distinct(dose_mg) >= 2) %>%
    ungroup()
  
  preds <- map_dfr(split(es2, interaction(es2$molecule, es2$ae_term, drop = TRUE)), function(g){
    doses <- make_grid(g$dose_mg, grid = grid, n_grid = n_grid)
    m <- fit_dr(g, model = model, df_spline = df_spline)
    if (is.null(m)) return(tibble())
    
    predict_dr(m, doses, model = model) %>%
      mutate(molecule = g$molecule[1], ae_term = g$ae_term[1], model = model)
  })
  
  models <- es2 %>%
    split(interaction(es2$molecule, es2$ae_term, drop = TRUE)) %>%
    map_dfr(function(g){
      m <- fit_dr(g, model = model, df_spline = df_spline)
      if (is.null(m)) return(tibble())
      
      mm <- if (model == "spline") m$model else m
      
      tibble(
        molecule  = g$molecule[1],
        ae_term   = g$ae_term[1],
        model     = model,
        df_spline = ifelse(model == "spline", df_spline, NA_integer_),
        k         = mm$k,
        beta      = if (model == "linear") as.numeric(coef(mm)[2]) else NA_real_,
        pval      = if (model == "linear") as.numeric(mm$pval[2]) else as.numeric(mm$QMp)
      )
    })
  
  list(preds = preds, models = models)
}