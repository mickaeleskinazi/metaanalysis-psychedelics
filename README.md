# Meta-analysis of Psychedelics – Adverse Events  

This repository contains data and analysis scripts for a **meta-analysis of adverse events (AEs) in psychedelic studies**.  
The main goal is to explore **dose–response relationships** for different molecules and adverse effects, using the [`metafor`](https://cran.r-project.org/package=metafor) package in R.  

---

## 📊 Overview  
- **Objective**: Assess how psychedelic molecules (MDMA, LSD, psilocybin, ayahuasca, etc.) influence adverse events.  
- **Approach**:  
  - Dose–response analysis by molecule across all studies and time windows (session vs. follow-up).  
  - Dose–response analysis by molecule × adverse effect.  
  - Forest plots, dose–response curves, and summary tables with model estimates and significance.  

---

## 📂 Project structure

```
metaanalysis-psychedelics/
├── data/                     # Raw Excel inputs
├── R/                        # Reusable analysis modules (data prep, modelling, plotting)
├── scripts/                  # Entry-point scripts orchestrating analyses
├── results*/                 # Generated outputs (tables, plots, publication artefacts)
├── README.md
└── ...
```

> **Why two folders with R code?**
>
> * `R/` behaves like an internal package: it stores **functions** (helpers for loading
>   data, fitting models, making plots, exporting tables, etc.). These files are not
>   meant to be sourced one by one by hand.
> * `scripts/` contains the **entry points** that orchestrate the end-to-end workflows.
>   Each script starts by sourcing the modules under `R/` and then runs a complete
>   analysis (main follow-up, session vs follow-up comparison, global slope check,
>   overlay plots, …).

Key entry-point scripts:

- `scripts/run_main_analysis.R` – end-to-end pipeline for the main dataset (follow-up window).
- `scripts/run_session_followup_analysis.R` – compares session vs follow-up time windows and saves per-window artefacts.
- `scripts/compare_global_session_followup.R` – global slope comparison between windows.
- `scripts/compare_session_followup_from_saved_tables.R` – reconciles previously exported tables for reporting.
- `scripts/analysis_plots_by_ae_overlay.R` – generates AE overlays when you need molecule comparisons.


---

## ⚙️ Methods  

### Analyses  
1. **Global dose–response**  
   - By molecule, across all studies/time windows.  
2. **Specific adverse effects**  
   - Dose–response by molecule × AE.  
   - Test significance of dose effect on each AE.  
3. **Outputs**  
   - Forest plots (odds ratios).  
   - Dose–response curves.  
   - Tables of estimates, 95% CI, p-values, significance stars (`*`, `**`, `***`).  

---

## 🛠️ Requirements

- R (≥ 4.2.0)
- Suggested packages (install with `install.packages()`):
  - `metafor`, `dplyr`, `tidyr`, `readxl`, `ggplot2`, `purrr`, `janitor`, `stringr`, `here`,
    `patchwork`, `gt`, `kableExtra`, `flextable`, `officer`

## 🚀 How to run

### 1. Install R packages once

Open R (or RStudio) in the project root and install the required packages:

```r
install.packages(c(
  "metafor", "dplyr", "tidyr", "readxl", "ggplot2", "purrr",
  "janitor", "stringr", "here", "patchwork", "gt", "kableExtra",
  "flextable", "officer"
))
```

### 2. Choose how you want to launch the workflows

You can run the entry-point scripts from **within R** (interactive) or from the
**command line** with `Rscript`. The scripts automatically source the helpers in
`R/`, so you only need to call one command per workflow.

| Goal | Run from an R console | Run from a terminal |
| --- | --- | --- |
| Main follow-up analysis | `source("scripts/run_main_analysis.R")` | `Rscript scripts/run_main_analysis.R` |
| Session vs follow-up comparison | `source("scripts/run_session_followup_analysis.R")` | `Rscript scripts/run_session_followup_analysis.R` |
| Optional global slope check | `source("scripts/compare_global_session_followup.R")` | `Rscript scripts/compare_global_session_followup.R` |

Each script exposes arguments (see the top of the file) in case you need to point
to a different Excel file or results directory, but the defaults reproduce the
current study.

### 3. Confirm everything finished

1. Watch the console: the scripts print progress messages such as
   `"1) Load & harmonize …"`, `"Comparative dose–response overlays …"`, etc.
2. Inspect the outputs:
   - `results/main/…` – follow-up analysis artefacts.
   - `results_session/…` and `results_followup/…` – window-specific artefacts.
   - `results_compare/…` – side-by-side plots and tables for session vs follow-up.
   - `results/paper_tables/…` – publication-ready tables (when enabled).
3. Open key CSV/PDF files (e.g. `results_compare/tables/dr_session_followup_publication_table.csv`
   or `results_compare/dose_response/dr_session_vs_followup.pdf`) to verify that
   estimates, confidence intervals, and plots have been generated.

If a script stops with an error, the message will indicate the missing package,
file, or column so you can address it and rerun the command.

## 📈 Outputs

`results/main/` – artefacts from the primary follow-up analysis (tables, forest plots, dose–response curves, master figures).

`results_session/` & `results_followup/` – window-specific contrasts, models and plots.

`results_compare/` – side-by-side plots, comparative tables, and publication-ready summaries for session vs follow-up.

`results/paper_tables/` – formatted tables for manuscripts or slide decks.

## 🔮 Suggested enhancements

- Add study-level moderators (e.g., dosing paradigm, psychotherapy support, participant diagnosis) to explore heterogeneity.
- Incorporate risk-of-bias assessments and conduct sensitivity analyses excluding high-risk studies.
- Harmonise adverse-event terminology using controlled vocabularies (MedDRA) to ease cross-study comparisons.
- Extend the database with emerging molecules (e.g., DMT, mescaline) and longitudinal outcomes beyond acute/follow-up windows.
- Export machine-readable metadata (JSON) for downstream dashboards or reproducible manuscripts (e.g., Quarto bookdown).

📄 License

MIT License
