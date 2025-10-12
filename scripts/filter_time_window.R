#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(writexl)
})

# ---- CONFIG ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript scripts/filter_time_window.R <infile.xlsx> <sheet>\n",
       "Example: Rscript scripts/filter_time_window.R data/Adverse-events-dose-v5.xlsx Feuil1")
}

infile <- args[1]
sheet  <- args[2]

# ---- LOAD DATA ----
message("Reading Excel: ", infile)
df <- readxl::read_excel(infile, sheet = sheet)

if (!"time_window" %in% names(df)) {
  stop("Column 'time_window' not found in Excel sheet.")
}

# ---- NORMALIZE ----
df <- df %>%
  mutate(
    time_window = str_trim(tolower(time_window))
  )

# ---- SPLIT BY WINDOW ----
for (tw in c("session", "follow_up")) {
  df_tw <- df %>% filter(time_window == tw)
  if (nrow(df_tw) == 0) {
    message("⚠️ No rows found for ", tw, ". Skipping.")
    next
  }
  
  outfile <- sub("\\.xlsx$", paste0("_", tw, ".xlsx"), infile)
  writexl::write_xlsx(df_tw, outfile)
  message("✅ Saved: ", outfile, " (", nrow(df_tw), " rows)")
}

message("Done. Files generated for session and follow_up.")