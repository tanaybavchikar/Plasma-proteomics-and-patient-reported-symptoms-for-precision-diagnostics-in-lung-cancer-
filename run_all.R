# =============================================================================
# run_all.R  —  the numbered pipeline (00–07), in order.
# -----------------------------------------------------------------------------
# Working directory MUST be LC_Data.  Run:  source("../code/run_all.R")
# All paths are relative to LC_Data (not to the script), so the layout works
# as long as wd = LC_Data and helpers stay in code/helpers/.
#
# Full pipeline takes hours (the tuning steps dominate). Comment out what you
# don't need. Dependencies: 01 -> 02 -> {03,04} -> 05 -> {06,07}. 00 is
# descriptive and needs 01+02's outputs. 08 and 09 are not built yet.
# =============================================================================

run <- function(path) {
  message("\n=====================  ", path, "  =====================")
  source(path, chdir = FALSE)          # keep wd = LC_Data throughout
}

run("../code/01_clean_data.R")             # raw -> cleaned data (symptoms, LCP1, LCP2)
run("../code/02_create_ML_data.R")         # cleaned -> ML tables + CV folds
run("../code/00_univariate_table1.R")      # Task 1 (descriptive; needs 01+02 outputs)
run("../code/03_symptombg_411_model.R")    # EN+RF on 411 symptoms+bg (for selection)
run("../code/04_LCP2_proteomics_model.R")  # EN+RF on LCP2 proteins (for selection)
run("../code/05_select_features.R")        # stability selection -> top proteins + vars
run("../code/06_symptombg_87_model.R")     # symptoms+bg-only model on the 87
run("../code/07_proteomics_87_model.R")    # proteomics-only model on the 87
run("../code/08_early_integrated_87_model.R")  # combined (proteins + symptoms+bg) on the 87
run("../code/09_compare_performance.R")        # compare the 3 approaches (ROC + DeLong)
run("../code/10_importance_cv_vs_all.R")       # importance: CV folds vs all-87 (Venn + final lists)

message("\nPipeline 00-10 complete.")
