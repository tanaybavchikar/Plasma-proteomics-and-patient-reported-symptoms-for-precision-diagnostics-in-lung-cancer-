# ============================================================
# 09_compare_performance.R  — compare the 3 modelling approaches on the 87 cohort
# ------------------------------------------------------------
# Reads the outputs of models_87.R (== 06 / 07 / 08) for both selection criteria
# and keeps ONE criterion per model family (MODEL_CRIT below): EN from impfrac,
# RF from selfreq. So each model family is shown at its chosen criterion, giving
# one EN + one RF summary. Change MODEL_CRIT to swap which criterion each uses.
#
#   06_symptombg_87        -> symptoms+background
#   07_proteomics_87       -> proteomics
#   08_early_integrated_87 -> integrated
#
# All models trained EN + RF on the SAME 87 patients with the SAME inline folds
# (seed 900 / v=5 / repeats=10 / patient order), so out-of-fold predictions align
# by .row across approaches AND criteria -> averaging and paired DeLong are valid.
#
# Writes (results/09_compare_performance/):
#   cv_performance_all.csv   CV metrics (averaged over criteria) per model x approach
#   auc_pooled_oof.csv       pooled out-of-fold AUC per model_set
#   roc_EN.png / roc_RF.png  ROC curves, one figure per family (3 approaches each)
#   auc_barplot.png          CV ROC-AUC bar chart, EN vs RF, WITH numbers
#   metrics_EN.png / metrics_RF.png   grouped metric bars, coloured by modality
#   delong_pairwise.csv      paired DeLong tests between approaches (per family)
#
# NOTE: impfrac and selfreq select different feature sets, so EN-impfrac and
# RF-selfreq are simply the two chosen models; nothing is averaged. Every plot
# (ROC per family, metric bars, AUC bar, DeLong) is built from exactly these.
#
# Run with working directory = LC_Data, after models_87.R has produced both criteria.
# ============================================================

required <- c("tidymodels")
to_install <- required[!required %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)
library(tidymodels)
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")

OUT      <- "results/09_compare_performance"
CRITERIA <- c("_impfrac", "_selfreq")   # both read in, then each family keeps ONE (below)
# each model family keeps ONE criterion (NOT averaged): EN from impfrac, RF from selfreq
MODEL_CRIT <- c(EN = "_impfrac", RF = "_selfreq")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

approaches <- c(
  "symptoms+background" = "results/06_symptombg_87",
  "proteomics"          = "results/07_proteomics_87",
  "integrated"          = "results/08_early_integrated_87")

# split "EN · proteomics" -> model = "EN", modality = "proteomics"
split_ms <- function(df) tidyr::separate(df, model_set, into = c("model", "modality"),
                                         sep = " · ", remove = FALSE)

## ---- 0) read a given file base across ALL criteria x approaches ---------------
read_across <- function(fname_base) {
  purrr::map_dfr(CRITERIA, function(crit) {
    purrr::imap_dfr(approaches, function(dir, appr) {
      f <- file.path(dir, paste0(fname_base, crit, ".csv"))
      if (!file.exists(f)) { message("SKIP (missing): ", f); return(NULL) }
      read.csv(f) %>% dplyr::mutate(approach = appr, criterion = crit)
    })
  })
}

cv_raw  <- read_across("cv_performance")
oof_raw <- read_across("oof_predictions")
if (is.null(cv_raw) || !nrow(cv_raw))
  stop("No cv_performance", paste(CRITERIA, collapse = "/"), ".csv found under ",
       "results/06_*, 07_*, 08_*. Run models_87.R first.")

## ---- 1) keep ONE criterion per model family (EN<-impfrac, RF<-selfreq) --------
# no averaging: each model_set comes from exactly one criterion
pick <- function(df) df %>% split_ms() %>%
  dplyr::filter(criterion == MODEL_CRIT[model]) %>%
  dplyr::select(-model, -modality)               # rebuilt later by split_ms()

cv_all  <- pick(cv_raw)
write.csv(cv_all, file.path(OUT, "cv_performance_all.csv"), row.names = FALSE)

oof_all <- pick(oof_raw) %>%
  dplyr::mutate(Lungcancer = factor(Lungcancer, levels = c("Yes", "No")))

cat("\nCV ROC-AUC by model x approach (EN: impfrac, RF: selfreq):\n")
print(as.data.frame(cv_all %>% dplyr::filter(.metric == "roc_auc") %>%
                      dplyr::select(model_set, approach, criterion, mean, std_err)), row.names = FALSE)

## ---- 2) pooled OOF AUC per model_set (drives the ROC legends) -----------------
aucs <- oof_all %>% dplyr::group_by(model_set) %>%
  roc_auc(truth = Lungcancer, .pred_Yes) %>%
  dplyr::ungroup() %>% split_ms() %>%
  dplyr::mutate(legend = sprintf("%s  (AUC %.3f)", modality, .estimate)) %>%
  dplyr::select(model_set, auc = .estimate, legend)
write.csv(aucs, file.path(OUT, "auc_pooled_oof.csv"), row.names = FALSE)

## ---- 3) ROC curves — one figure per family (EN, RF) --------------------------
roc_dat <- oof_all %>% dplyr::group_by(model_set) %>%
  roc_curve(truth = Lungcancer, .pred_Yes) %>%
  dplyr::left_join(dplyr::select(aucs, model_set, legend), by = "model_set") %>%
  split_ms()

plot_roc_family <- function(fam) {
  d <- dplyr::filter(roc_dat, model == fam)
  if (!nrow(d)) { message("no ROC data for ", fam); return(invisible()) }
  p <- ggplot(d, aes(1 - specificity, sensitivity, colour = legend)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_path(linewidth = 0.9) + coord_equal() +
    labs(title = sprintf("LCP1 87 cohort — %s: symptoms+bg vs proteomics vs integrated", fam),
         subtitle = "Pooled out-of-fold predictions (EN from impfrac, RF from selfreq)",
         x = "1 - specificity", y = "sensitivity", colour = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
  ggsave(file.path(OUT, paste0("roc_", fam, ".png")), p, width = 8, height = 7.2, dpi = 150)
}
plot_roc_family("EN")
plot_roc_family("RF")

## ---- colours shared by the AUC bar chart and the metric bars -----------------
mod_cols <- c("symptoms+background" = "#4CA9A6",   # teal
              "proteomics"          = "#E7A6AE",   # pink
              "integrated"          = "#B0A4CE")   # purple

## ---- 4) CV ROC-AUC bar chart: clustered by EN / RF, coloured by approach ------
bar_dat <- cv_all %>% dplyr::filter(.metric == "roc_auc") %>% split_ms() %>%
  dplyr::mutate(model    = factor(model, levels = c("EN", "RF")),
                approach = factor(approach,
                           levels = c("symptoms+background", "proteomics", "integrated")))
p_bar <- ggplot(bar_dat, aes(model, mean, fill = approach)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean - std_err, ymax = mean + std_err),
                position = position_dodge(0.8), width = 0.2) +
  geom_text(aes(label = sprintf("%.2f", mean), y = mean + std_err),
            position = position_dodge(0.8), vjust = -0.5, size = 3.4) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
  scale_fill_manual(values = mod_cols, name = "Approach") +
  coord_cartesian(ylim = c(0.4, 1)) +
  labs(title = "CV ROC-AUC: EN vs RF by approach (87 cohort)", x = NULL,
       y = "CV ROC-AUC (mean +/- SE)") +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(OUT, "auc_barplot.png"), p_bar, width = 7.5, height = 5.2, dpi = 150)

## ---- 5) grouped metric bars per family, coloured by modality -----------------
# like the reference figure: metric on x, one coloured bar per modality, error bars
metric_labels <- c(roc_auc = "ROC-AUC", sens = "Sensitivity", spec = "Specificity",
                   bal_accuracy = "Balanced acc.", j_index = "Youden's J")

metr <- cv_all %>% split_ms() %>%
  dplyr::filter(.metric %in% names(metric_labels)) %>%
  dplyr::mutate(metric   = factor(metric_labels[.metric], levels = unname(metric_labels)),
                modality = factor(modality,
                           levels = c("symptoms+background", "proteomics", "integrated")))

plot_metrics_family <- function(fam) {
  d <- dplyr::filter(metr, model == fam)
  if (!nrow(d)) { message("no metrics for ", fam); return(invisible()) }
  p <- ggplot(d, aes(metric, mean, fill = modality)) +
    geom_col(position = position_dodge(0.8), width = 0.72) +
    geom_errorbar(aes(ymin = pmax(0, mean - std_err), ymax = pmin(1, mean + std_err)),
                  position = position_dodge(0.8), width = 0.25) +
    scale_fill_manual(values = mod_cols, name = "Modality") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = sprintf("%s — CV metrics by modality (87 cohort)", fam),
         subtitle = "mean +/- SE  (EN: impfrac, RF: selfreq)", x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 20, hjust = 1))
  ggsave(file.path(OUT, paste0("metrics_", fam, ".png")), p, width = 8.5, height = 5.4, dpi = 150)
}
plot_metrics_family("EN")
plot_metrics_family("RF")

## ---- 6) paired DeLong tests between approaches (within each family) -----------
if (requireNamespace("pROC", quietly = TRUE)) {
  wide <- oof_all %>% dplyr::select(.row, model_set, .pred_Yes) %>%
    tidyr::pivot_wider(names_from = model_set, values_from = .pred_Yes) %>%
    dplyr::arrange(.row)
  truth <- oof_all %>% dplyr::distinct(.row, Lungcancer) %>%
    dplyr::arrange(.row) %>% dplyr::pull(Lungcancer)

  delong <- function(a, b) {
    if (!all(c(a, b) %in% names(wide))) { message("skip DeLong (missing): ", a, " / ", b); return(NULL) }
    r1 <- pROC::roc(truth, wide[[a]], quiet = TRUE, levels = c("No", "Yes"), direction = "<")
    r2 <- pROC::roc(truth, wide[[b]], quiet = TRUE, levels = c("No", "Yes"), direction = "<")
    tt <- pROC::roc.test(r1, r2, method = "delong", paired = TRUE)
    tibble(comparison = paste(a, "vs", b),
           auc_1 = as.numeric(pROC::auc(r1)), auc_2 = as.numeric(pROC::auc(r2)),
           Z = unname(tt$statistic), p_value = tt$p.value)
  }
  pairs <- list(
    c("EN · proteomics",  "EN · symptoms+background"),
    c("EN · integrated",  "EN · proteomics"),
    c("EN · integrated",  "EN · symptoms+background"),
    c("RF · proteomics",  "RF · symptoms+background"),
    c("RF · integrated",  "RF · proteomics"),
    c("RF · integrated",  "RF · symptoms+background"))
  dl <- purrr::map_dfr(pairs, ~ delong(.x[1], .x[2]))
  if (nrow(dl)) {
    write.csv(dl, file.path(OUT, "delong_pairwise.csv"), row.names = FALSE)
    cat("\nPaired DeLong tests (same 87 patients; EN from impfrac, RF from selfreq):\n")
    print(as.data.frame(dl), row.names = FALSE)
  }
} else {
  message("pROC not installed - skipping DeLong (install.packages('pROC'))")
}

message("\nDone. Comparison (EN from ", MODEL_CRIT["EN"], ", RF from ", MODEL_CRIT["RF"],
        ") saved in ", OUT, "/")
