suppressPackageStartupMessages({ library(dplyr) })

# -- chemins --
DATA_XLSX <- "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/data/Adverse-events-dose-v5.xlsx"
SHEET     <- "Feuil1"
ROOT      <- "/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/results"


default_ref_policies <- list(
  MDMA       = c("inactive_placebo","active_non_psy_placebo","active_placebo"),
  LSD        = c("inactive_placebo","active_placebo","active_non_psy_placebo"),
  PSILOCYBIN = c("inactive_placebo","active_non_psy_placebo","active_placebo"),
  AYAHUASCA  = c("inactive_placebo","active_placebo","active_non_psy_placebo"),
  .default   = c("inactive_placebo","active_placebo","active_non_psy_placebo")
)

# -- params --
MIN_K      <- 3
FIT_SPLINE <- TRUE

# -- sources --
source("/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/scripts/utils_data.R")
source("/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/scripts/analysis_dose_response.R")
source("/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/scripts/analysis_forest_plots.R")
source("/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/scripts/analysis_plots_dr.R")
source("/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/scripts/analysis_tables.R")
source("/Users/mickaeleskinazi/Documents/GitHub/metaanalysis-psychedelics/scripts/analysis_plots_master.R")

# -- pipeline --
message("1) Load & harmonize …")
raw <- load_data(DATA_XLSX, sheet = SHEET)

message("2) Build pairwise contrasts (reference-aware) …")
contr <- build_pairwise_2x2(raw, ref_policies = default_ref_policies)

# (optionnel) audit LSD
# audit_doses(raw, contr, "LSD")

message("3) Escalc …")
es <- build_escalc(contr)

message("4) DR by molecule (all AE, all windows) …")
dr_mol <- run_dr_by_molecule(es, min_k = MIN_K, fit_spline = FIT_SPLINE)

message("5) DR by AE (curves per molecule) …")
dr_ae  <- run_dr_by_ae(es, min_k = MIN_K, fit_spline = FIT_SPLINE)

dir.create(file.path(ROOT, "tables"), recursive = TRUE, showWarnings = FALSE)

message("6) Save significance tables …")
save_significance_table(dr_mol$models, file.path(ROOT,"tables/significance_by_molecule_models.csv"))
save_significance_table(dr_ae$models,  file.path(ROOT,"tables/significance_by_ae_models.csv"))
save_agg_significance_by_molecule(dr_mol$models, file.path(ROOT,"tables/significance_agg_by_molecule.csv"))
save_agg_significance_by_ae_molecule(dr_ae$models,  file.path(ROOT,"tables/significance_agg_by_ae_molecule.csv"))

message("7) Forest plots …")
make_forest_plots(es, file.path(ROOT,"forest_plots"))
make_forest_plots_per_molecule_pdf(es, file.path(ROOT,"forest_plots_by_molecule"))
make_forest_summary_per_molecule(es, file.path(ROOT,"forest_plots_summary"))

message("8) Dose–response plots …")
plot_dr_by_molecule_split(dr_mol$preds, dr_mol$models, file.path(ROOT,"dose_response/by_molecule_split"))
plot_dr_per_molecule_across_ae_facets(dr_ae$preds, dr_ae$models, file.path(ROOT,"dose_response/by_molecule_facets"))
plot_dr_per_ae_normalized_dose(dr_ae$preds, file.path(ROOT,"dose_response/by_ae_normalized"))

message("MASTER PLOTS …")
# A) DR global par molécule
plot_master_dr_by_molecule(dr_mol$preds, file.path(ROOT,"master/master_dr_by_molecule.pdf"))
# B) Forest globaux (pooled par AE, facets par molécule)
plot_master_forest_by_molecule(es, file.path(ROOT,"master/master_forest_by_molecule.pdf"), min_k = 2)
# C) DR par AE (grille)
plot_master_dr_by_ae(
  preds               = dr_ae$preds,
  models              = dr_ae$models,
  outfile             = file.path(ROOT, "master/master_dr_by_ae.pdf"),
  max_ae_per_molecule = 20,
  significant_only    = TRUE
)


message("MASTER PLOTS … done ✅")
message("Done ✅")
