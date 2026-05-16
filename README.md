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


### Outcome definition (for manuscript reporting)
Use the following **complete Outcome text** in Methods (before Statistical analyses and Results):

> **Outcome.** The primary outcome was the odds of adverse-event (AE) occurrence as a function of psychedelic dose, estimated from arm-level AE counts across randomized controlled trials. To avoid combining clinically distinct temporal profiles, outcomes were pre-specified and analysed separately according to assessment window: acute/session outcomes (AEs occurring during administration or in the immediate post-session period) and follow-up outcomes (AEs reported at later post-acute assessments). The same dose–response meta-analytic framework was applied independently within each time window at both the global molecule level (pooling all extractable AEs within each substance) and the molecule × AE level (modelling individual harmonized AE terms). AE-level models were fit separately for each molecule (LSD, MDMA, psilocybin) and were not pooled into a single pan-substance effect. Multiplicity-adjusted significance reporting was applied. Session-versus-follow-up comparisons were treated as secondary analyses and are reported in dedicated tables and figures.


Suggested wording when placebo composition is partially missing:
- **Placebo condition**: Placebo composition was not reported in several trials. Among studies with known composition, inactive placebos (lactose capsules or saline) were used in both ayahuasca trials (2/2, 100%), 6 MDMA studies, 3 LSD studies, and 1 psilocybin study.

Repository mapping for this outcome structure:
- Session outputs: `results/main/session/` and `results/paper_tables/session/`.
- Follow-up outputs: `results/main/follow_up/` and `results/paper_tables/follow_up/`.
- Cross-window comparisons: `results/main/compare/` (especially `tables/topline_session_followup_summary.csv` and `tables/dr_session_followup_publication_table.csv`).

---


### Sensitivity analyses (for manuscript Methods)
Use the following text in Methods to match the Results statements:

> **Sensitivity analyses.** We evaluated model-specification robustness by refitting all dose–response analyses under two prespecified functional forms: (i) a linear dose term and (ii) a restricted cubic spline specification (natural spline, 3 degrees of freedom). This was done separately within each assessment window (session and follow-up), both at the global substance level and at the molecule × AE level, using the same eligibility thresholds as the primary analyses. For each analysis unit, we compared significance patterns across specifications and classified results as concordant (both significant), linear-only, spline-only, or non-significant. Robustness summaries were reported in Supplementary Table S1 (global substance level) and Supplementary Table S2 (AE level counts by molecule/window).

Operational mapping in this repository:
- Robustness inputs are generated in `scripts/run_main_analysis.R` and saved per window under `results/main/<window>/robustness/`.
- Main files: `robustness_linear_vs_spline_by_molecule.csv` and `robustness_linear_vs_spline_by_ae_molecule.csv`.
- Supplement-ready summary tables are produced by `R/robustness_report_tables.R` (S1/S2 exports).

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

# Regenerate manuscript-ready LaTeX tables and narrative (no model refit)
Rscript scripts/build_results_tables.R
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
