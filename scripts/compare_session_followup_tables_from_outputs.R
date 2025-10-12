#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
})

source(here::here("scripts", "compare_session_followup_from_saved_tables.R"))

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(idx, default) {
  if (length(args) >= idx && nzchar(args[[idx]])) {
    args[[idx]]
  } else {
    default
  }
}

# Allow running against the legacy directory structure created by run_main_analysis.
# Users can override the defaults by passing explicit paths for the session tables,
# follow-up tables, and output directory (in that order).
dir_session <- get_arg(1, here::here("results_session", "tables"))
dir_followup <- get_arg(2, here::here("results_followup", "tables"))
out_dir <- get_arg(3, here::here("results_compare", "tables"))

message("Reading session tables from: ", dir_session)
message("Reading follow-up tables from: ", dir_followup)
message("Writing comparison tables to: ", out_dir)

compare_session_followup_from_saved_tables(
  dir_session = dir_session,
  dir_followup = dir_followup,
  out_dir = out_dir
)
