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
├── results/                  # Generated outputs (tables, plots, publication artefacts)
├── README.md
└── ...
```

Key scripts:

- `scripts/run_main_analysis.R` – end-to-end pipeline for the main dataset (follow-up window).
- `scripts/run_session_followup_analysis.R` – compares session vs follow-up time windows and saves per-window artefacts.
- `scripts/compare_global_session_followup.R` – global slope comparison between windows.
- `scripts/compare_session_followup_from_saved_tables.R` – reconciles previously exported tables for reporting.
- `scripts/compare_session_followup_tables_from_outputs.R` – builds the session vs follow-up comparison tables from existing
  CSV exports (no model re-fitting required).
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

From the project root you can execute the entry-point scripts directly with `Rscript`:

```bash
# Main analysis (fits the models and generates plots/tables)
Rscript scripts/run_main_analysis.R

# Session vs follow-up analysis with fresh model fits
Rscript scripts/run_session_followup_analysis.R

# Session vs follow-up tables only (reuses saved CSVs)
# Optional arguments: <session_dir> <followup_dir> <output_dir>
Rscript scripts/compare_session_followup_tables_from_outputs.R

# Session vs follow-up tables only (reuse saved CSVs)
source("scripts/compare_session_followup_tables_from_outputs.R")

# Global slope comparison (optional)
Rscript scripts/compare_global_session_followup.R
```

Each script accepts parameters (see function definitions) so you can point to alternative files or output folders. The
session/follow-up table helper defaults to the legacy `results_session/`, `results_followup/`, and `results_compare/`
directories, but you can override them by supplying explicit paths when calling the script.

### 🔍 Verifying the pipeline end-to-end

1. **Install packages** listed in the requirements section (`install.packages(c("metafor", "here", ...))`).
2. **Run the scripts** above from an interactive R session or with `Rscript scripts/run_main_analysis.R` (etc.).
3. **Inspect the console log** – each major step prints a numbered progress message so you can see where the pipeline is.
4. **Check the output folders**:
   - `results/main/` (follow-up analysis) – should contain `tables/`, `forest_plots/`, `dose_response/`, and `master/` sub-folders.
   - `results/session/`, `results/follow_up/`, and `results/compare/` – confirm comparable artefacts exist for the session vs. follow-up workflow.
   - `results/paper_tables/` – verify CSV/Word/HTML tables for manuscript use when `make_paper_tables = TRUE`.
5. **Review generated CSVs** (e.g., `results/compare/tables/dr_session_followup_publication_table.csv`) to confirm model estimates are populated and that `stars` columns mark significant effects.

If any of the scripts stop with an error, the message will call out the missing file, package, or column that needs attention.

## 📈 Outputs

`results/main/` – artefacts from the primary follow-up analysis (tables, forest plots, dose–response curves, master figures).

`results/session/` & `results/follow_up/` – window-specific contrasts, models and plots.

`results/compare/` – side-by-side plots, comparative tables, and publication-ready summaries for session vs follow-up.

`results/paper_tables/` – formatted tables for manuscripts or slide decks.

## 🔮 Suggested enhancements

- Add study-level moderators (e.g., dosing paradigm, psychotherapy support, participant diagnosis) to explore heterogeneity.
- Incorporate risk-of-bias assessments and conduct sensitivity analyses excluding high-risk studies.
- Harmonise adverse-event terminology using controlled vocabularies (MedDRA) to ease cross-study comparisons.
- Extend the database with emerging molecules (e.g., DMT, mescaline) and longitudinal outcomes beyond acute/follow-up windows.
- Export machine-readable metadata (JSON) for downstream dashboards or reproducible manuscripts (e.g., Quarto bookdown).

📄 License

MIT License
