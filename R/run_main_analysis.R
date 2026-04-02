# Backward-compatible entrypoint for code that previously sourced
# R/run_main_analysis.R directly.
#
# The canonical implementation now lives in scripts/run_main_analysis.R.
source(here::here("scripts", "run_main_analysis.R"))
