#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
})

# Delegate to the canonical pipeline in R/run_main_analysis.R so the CLI entry
# point always runs the full analysis workflow (significance, robustness,
# publication tables, session/follow-up comparison outputs, etc.).
runner_env <- new.env(parent = globalenv())
sys.source(here::here("R", "run_main_analysis.R"), envir = runner_env)

runner_env$run_main_analysis()
