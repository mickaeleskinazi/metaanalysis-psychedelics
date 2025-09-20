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


metaanalysis-psychedelics/
│
├── data/ # Raw and processed datasets
│ ├── psychedelic_adverse_events.xlsx
│ └── README.md
│
├── scripts/ # R scripts for analysis
│ ├── run_analysis.R
│ ├── utils_data.R
│ ├── analysis_dose_response.R
│ ├── analysis_forest_plots.R
│ └── analysis_tables.R
│
├── results/ # Outputs
│ ├── tables/
│ │ └── significance_results.csv
│ ├── forest_plots/
│ ├── dose_response/
│ └── global_summary/
│
├── README.md # Project description
├── LICENSE
└── .gitignore


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
- R packages:  
  ```r
  install.packages(c("metafor", "dplyr", "tidyr", "readxl", "ggplot2"))

🚀 How to run

Clone the repository:git clone https://github.com/your-username/metaanalysis-psychedelics.git
cd metaanalysis-psychedelics
Run the main script in R:source("scripts/run_analysis.R")
Results (tables + plots) will be saved in results/.

📈 Outputs

results/tables/ – Statistical tables (estimates, CI, p-values, significance).

results/forest_plots/ – Forest plots per AE.

results/dose_response/ – Dose–response plots by molecule.

results/global_summary/ – Overview graphics of molecules × AEs.

📄 License

MIT License
