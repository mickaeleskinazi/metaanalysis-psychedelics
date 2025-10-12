suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(purrr); library(stringr); library(tidyr)
})

# try to unnest a common list-col if present
.unlist_preds <- function(p) {
  if (!is.data.frame(p)) stop("preds is not a data.frame")
  # Unnest typical list columns if they exist
  for (lc in c("preds","data","df","tbl")) {
    if (lc %in% names(p) && is.list(p[[lc]])) {
      p <- tidyr::unnest(p, tidyselect::all_of(lc))
    }
  }
  p
}

# normalize any predictions df to: molecule, x, pred, lwr, upr
.normalize_pred_cols <- function(p){
  p <- .unlist_preds(p)
  nm <- names(p)
  
  # find molecule col
  molcol <- if ("molecule" %in% nm) "molecule" else NULL
  if (is.null(molcol)) {
    stop("Missing 'molecule' column. Columns are: ", paste(nm, collapse=", "))
  }
  
  # find dose axis
  x_candidates <- c("xdose","dose","dose_mg","x","dose_std","dose_norm")
  xcol <- intersect(x_candidates, nm)
  if (length(xcol) == 0) {
    stop("Could not find a dose column among: ", paste(x_candidates, collapse=", "),
         ". Columns are: ", paste(nm, collapse=", "))
  }
  xcol <- xcol[1]
  
  # find pred & CI columns
  pred_sets <- list(
    c("pred","lwr","upr"),
    c("fit","ci_low","ci_high"),
    c("estimate","ci_low","ci_high"),
    c("fit","ci_lb","ci_ub"),
    c("estimate","ci_lb","ci_ub"),
    c("y","ylwr","yupr"),
    c("mu","lwr","upr"),
    c("fit","lwr","upr")    # ← add this line
  )
  found <- NULL
  for (cs in pred_sets) if (all(cs %in% nm)) { found <- cs; break }
  if (is.null(found)) {
    stop(paste0(
      "Could not find prediction columns. Need one of:\n",
      "  - pred,lwr,upr\n",
      "  - fit,ci_low,ci_high\n",
      "  - estimate,ci_low,ci_high\n",
      "  - fit,ci_lb,ci_ub\n",
      "  - estimate,ci_lb,ci_ub\n",
      "  - y,ylwr,yupr\n",
      "  - mu,lwr,upr\n",
      "Columns are: ", paste(nm, collapse=", ")
    ))
  }
  
  p %>%
    dplyr::rename(
      x    = !!xcol,
      pred = !!found[1],
      lwr  = !!found[2],
      upr  = !!found[3]
    ) %>%
    dplyr::transmute(
      molecule = .data[[molcol]],
      x        = suppressWarnings(as.numeric(.data$x)),
      pred     = suppressWarnings(as.numeric(.data$pred)),
      lwr      = suppressWarnings(as.numeric(.data$lwr)),
      upr      = suppressWarnings(as.numeric(.data$upr))
    ) %>%
    dplyr::filter(is.finite(x), is.finite(pred), is.finite(lwr), is.finite(upr))
}

# Compare session vs follow_up dose-response per molecule (overlay ribbons/lines)
dr_compare_all_molecules <- function(preds_session, preds_followup, outfile,
                                     width = 11, height = 8.5, dpi = 150) {
  if (is.null(preds_session) || !nrow(preds_session)) stop("Empty preds_session")
  if (is.null(preds_followup) || !nrow(preds_followup)) stop("Empty preds_followup")
  
  # normalize both tables
  ps <- .normalize_pred_cols(preds_session)  %>% mutate(time_window = "session")
  pf <- .normalize_pred_cols(preds_followup) %>% mutate(time_window = "follow_up")
  
  # keep only molecules present in BOTH windows
  keep_mol <- intersect(unique(ps$molecule), unique(pf$molecule))
  if (!length(keep_mol)) stop("No overlapping molecules between session and follow_up.")
  pdat <- bind_rows(ps, pf) %>% filter(molecule %in% keep_mol)
  
  # plot
  p <- ggplot(pdat, aes(x = x, y = pred, color = time_window, fill = time_window)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, color = NA) +
    geom_line(size = 1) +
    facet_wrap(~ molecule, scales = "free_x") +
    labs(
      title = "Dose–response: session vs follow-up (by molecule)",
      x = "Dose (mg or study scale)",
      y = "log(OR) vs reference"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  ggsave(outfile, p, width = width, height = height, dpi = dpi)
  invisible(p)
}