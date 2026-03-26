#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
})

# Lightweight CLI wrapper:
# - runs the canonical pipeline in R/run_main_analysis_core.R
# - keeps defaults (including run_reporting_tables = TRUE)
runner <- new.env(parent = globalenv())
sys.source(here::here("R", "run_main_analysis_core.R"), envir = runner)

runner$run_main_analysis()
