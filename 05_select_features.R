# ============================================================
# 05_select_features.R  — feature selection for the 87-cohort models
# ------------------------------------------------------------
# Reduces each modality (proteins from LCP2, symptoms+background from the 411) to
# a small feature set to carry into the 87-patient models. Selecting on the LARGE
# independent cohorts and evaluating on the 87 avoids selection leakage.
#
# Two things happen here:
#
# (1) BARPLOTS of three per-feature quantities, for each of the 4 model runs, so
#     you can see the distributions and sanity-check the cutoffs:
#       - mean importance          (mean_imp)
#       - scaled mean importance   (mean_imp / max mean_imp  -> 0..1)
#       - selection metric         (EN: selection_freq ; RF: freq_top10)
#
# (2) TWO selection criteria, each producing its own feature set (union of EN & RF):
#       - "impfrac": scaled mean importance >= IMP_FRAC_CUTOFF (default 0.10 = 10% of max)
#       - "selfreq": selection metric      >= SELFREQ_CUTOFF   (default 0.10)
#                    (EN uses selection_freq; RF uses freq_top10, because RF's raw
#                     selection_freq is ~1 for everything and useless.)
#     Written to top_variables_<modality>_impfrac.csv and _selfreq.csv, so
#     models_87.R can run all three approaches under each criterion via SUFFIX.
#
# Run with working dir = LC_Data, AFTER 03 (en/rf_symptbg) and 04 (en/rf_LCP2).
# ============================================================

library(tidyverse)

OUT <- "results/task3_variable_selection"
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

## ---- tunables (PER-MODALITY: loosen proteins without touching symptoms) -------
# criterion "impfrac": keep if scaled mean importance (mean_imp / max) >= cutoff
IMP_FRAC_CUTOFF_PROT <- 0.05   # proteins  (lower = more proteins kept)
IMP_FRAC_CUTOFF_SYMP <- 0.10   # symptoms + background
# criterion "selfreq": keep if selection metric (EN selection_freq / RF freq_top10) >= cutoff
SELFREQ_CUTOFF_PROT  <- 0.05   # proteins
SELFREQ_CUTOFF_SYMP  <- 0.10   # symptoms + background
COMBINE_RULE         <- "union"  # combine EN & RF: "union" (either) or "intersect" (both)

# the four model varimp files, with each model's selection-frequency column
models <- tibble::tribble(
  ~modality,             ~model, ~path,                                                 ~selfreq_col,
  "symptoms_background", "EN",   "results/en_symptbg/varimp_across_folds.csv",          "selection_freq",
  "symptoms_background", "RF",   "results/rf_symptbg/varimp_across_folds.csv",          "freq_top10",
  "proteins",            "EN",   "results/en_LCP2_proteomics/varimp_across_folds.csv",  "selection_freq",
  "proteins",            "RF",   "results/rf_LCP2_proteomics/varimp_across_folds.csv",  "freq_top10")

## ---- (1) distribution barplots -----------------------------------------------
bar_ranked <- function(d, valcol, title, out_png) {
  dd <- d %>% dplyr::filter(.data[[valcol]] > 0) %>%
    dplyr::arrange(dplyr::desc(.data[[valcol]])) %>%
    dplyr::mutate(rank = dplyr::row_number())
  p <- ggplot(dd, aes(rank, .data[[valcol]])) +
    geom_col(fill = "steelblue", width = 1) +
    labs(title = title, x = "feature rank (most -> least)", y = valcol) +
    theme_minimal()
  ggsave(out_png, p, width = 9, height = 5, dpi = 150)
}

for (i in seq_len(nrow(models))) {
  m <- models[i, ]
  if (!file.exists(m$path)) { message("SKIP barplots (missing): ", m$path); next }
  d <- read.csv(m$path)
  d$scaled_mean_imp <- d$mean_imp / max(d$mean_imp, na.rm = TRUE)
  tag <- sprintf("%s_%s", m$modality, m$model)
  bar_ranked(d, "mean_imp",        sprintf("%s: mean importance", tag),        file.path(OUT, sprintf("bar_meanimp_%s.png", tag)))
  bar_ranked(d, "scaled_mean_imp", sprintf("%s: scaled mean importance", tag), file.path(OUT, sprintf("bar_scaledimp_%s.png", tag)))
  bar_ranked(d, m$selfreq_col,     sprintf("%s: %s", tag, m$selfreq_col),      file.path(OUT, sprintf("bar_selfreq_%s.png", tag)))
}

## ---- (2) selection under each criterion --------------------------------------
# returns the kept rows for ONE model under ONE criterion, carrying both metrics
select_one <- function(m, criterion) {
  if (!file.exists(m$path)) { message("SKIP (missing): ", m$path); return(NULL) }
  d <- read.csv(m$path)
  if (!"label" %in% names(d)) d$label <- d$Variable
  d$scaled_mean_imp <- d$mean_imp / max(d$mean_imp, na.rm = TRUE)

  is_prot <- m$modality == "proteins"
  if (criterion == "impfrac") {
    d$score <- d$scaled_mean_imp
    cutoff  <- if (is_prot) IMP_FRAC_CUTOFF_PROT else IMP_FRAC_CUTOFF_SYMP
  } else {                              # "selfreq"
    d$score <- d[[m$selfreq_col]]
    cutoff  <- if (is_prot) SELFREQ_CUTOFF_PROT else SELFREQ_CUTOFF_SYMP
  }
  keep <- d$score >= cutoff
  cat(sprintf("  [%s] %-20s %-2s: %d of %d kept (cutoff %.2f)\n",
              criterion, m$modality, m$model, sum(keep, na.rm = TRUE), nrow(d), cutoff))
  d %>% dplyr::filter(keep) %>%
    dplyr::transmute(Variable, label, model = m$model,
                     score, mean_imp, scaled_mean_imp, selection_freq, freq_top10)
}

# for one criterion: union EN & RF per modality, write a wide file with both metrics
run_selection <- function(criterion, crit_suffix) {
  cat(sprintf("\n===== criterion: %s =====\n", criterion))
  for (mod in unique(models$modality)) {
    parts <- dplyr::bind_rows(lapply(which(models$modality == mod),
                                     function(i) select_one(models[i, ], criterion)))
    if (is.null(parts) || !nrow(parts)) { message("nothing selected for ", mod); next }

    by_model <- split(parts$Variable, parts$model)
    vars <- if (COMBINE_RULE == "intersect" && length(by_model) > 1)
              Reduce(intersect, by_model) else unique(unlist(by_model))
    sel  <- parts %>% dplyr::filter(Variable %in% vars)
    meta <- sel %>% dplyr::group_by(Variable) %>%
      dplyr::summarise(models = paste(sort(unique(model)), collapse = "+"), .groups = "drop")

    out <- sel %>%
      dplyr::select(Variable, label, model, score, mean_imp, scaled_mean_imp,
                    selection_freq, freq_top10) %>%
      tidyr::pivot_wider(names_from = model,
                         values_from = c(score, mean_imp, scaled_mean_imp,
                                         selection_freq, freq_top10)) %>%
      dplyr::left_join(meta, by = "Variable")

    f <- file.path(OUT, sprintf("top_variables_%s%s.csv", mod, crit_suffix))
    write.csv(out, f, row.names = FALSE)
    cat(sprintf("  -> %s : %d variables (%s of EN/RF)\n", basename(f), nrow(out), COMBINE_RULE))
  }
}

run_selection("impfrac", "_impfrac")   # scaled importance >= 10% of max
run_selection("selfreq", "_selfreq")   # selection metric  >= 10%

message("\nDone. Barplots + two selected sets in ", OUT, "/  ",
        "(top_variables_*_impfrac.csv and *_selfreq.csv).")
message("Feed into models_87.R (it loops both criteria via SUFFIX). ",
        "Adjust IMP_FRAC_CUTOFF / SELFREQ_CUTOFF if a set is too large/small.")
