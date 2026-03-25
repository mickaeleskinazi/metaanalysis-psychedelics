#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
})

# Canonical implementation lives in R/run_main_analysis_core.R
source(here::here("R", "run_main_analysis_core.R"))

if (identical(environment(), globalenv())) {
  run_main_analysis()
}
