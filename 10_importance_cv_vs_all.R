# ============================================================
# 10_importance_cv_vs_all.R  — compare variable importance: CV folds vs all-87 fit
# ------------------------------------------------------------
# models_87.R writes, per (model x modality x criterion), into each results folder:
#   varimp_final_<tag>.csv  importance from the model trained on ALL 87 (best HPs)
#   varimp_cv_<tag>.csv     importance AVERAGED across the CV folds (mean_imp, ...)
# where tag = "<modality>_<model><criterion>", e.g. "proteomics_EN_impfrac".
#
# This script joins the two, per feature, and asks: is a feature "important" the
# SAME way whether you read it off the all-87 fit or the fold-averaged CV fit?
# A feature is called IMPORTANT if its scaled importance (imp / max imp) >= CUTOFF.
#   - all-87 : scaled_all = importance / max(importance)
#   - CV     : scaled_cv  = mean_imp   / max(mean_imp)
# For each run it draws a 2-set Venn (all-87-important vs CV-important) and counts
# both / all-only / cv-only. Features important BOTH ways are the robust ones.
#
# Final lists (the deliverable): proteins from the proteomics runs (07), symptoms
# from the symptoms+background runs (06). A feature makes the list if it is
# important BOTH ways in >= 1 run; we also report in how many runs (EN/RF x
# impfrac/selfreq) it was robust, so you can sort by reproducibility.
#
# Run with working directory = LC_Data, AFTER models_87.R (with EXTRACT_CV_IMP=TRUE).
# ============================================================

library(tidyverse)

OUT <- "results/10_importance_cv_vs_all"
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

IMP_CUTOFF <- 0.10     # scaled importance (imp / max) to call a feature "important"

# every run models_87.R could have produced
runs <- tidyr::expand_grid(
  modality = c("symptoms+background", "proteomics", "integrated"),
  model    = c("EN", "RF"),
  crit     = c("_impfrac", "_selfreq")) %>%
  dplyr::mutate(
    dir = dplyr::recode(modality,
      "symptoms+background" = "results/06_symptombg_87",
      "proteomics"          = "results/07_proteomics_87",
      "integrated"          = "results/08_early_integrated_87"),
    tag = paste0(gsub("[^A-Za-z0-9]+", "_", modality), "_", model, crit))

## ---- 2-set Venn (no extra packages) -----------------------------------------
venn2 <- function(left_only, both, right_only, left_lab, right_lab, title, file) {
  circ <- function(x0, r = 1.5, n = 120) {
    t <- seq(0, 2 * pi, length.out = n); data.frame(x = x0 + r * cos(t), y = r * sin(t)) }
  L <- circ(-0.75); R <- circ(0.75)
  p <- ggplot() +
    geom_polygon(data = L, aes(x, y), fill = "#2E8C8A", alpha = 0.35, colour = "#2E8C8A", linewidth = 0.8) +
    geom_polygon(data = R, aes(x, y), fill = "#6D4C91", alpha = 0.35, colour = "#6D4C91", linewidth = 0.8) +
    annotate("text", x = -1.55, y = 0, label = left_only,  size = 7) +
    annotate("text", x =  0.00, y = 0, label = both,       size = 7, fontface = "bold") +
    annotate("text", x =  1.55, y = 0, label = right_only, size = 7) +
    annotate("text", x = -0.9, y = 1.85, label = left_lab,  size = 4.2, colour = "#1f6664") +
    annotate("text", x =  0.9, y = 1.85, label = right_lab, size = 4.2, colour = "#523670") +
    coord_equal() + theme_void() +
    labs(title = title) + theme(plot.title = element_text(hjust = 0.5, size = 12))
  ggsave(file, p, width = 6, height = 5, dpi = 150)
}

## ---- compare ONE run ---------------------------------------------------------
compare_one <- function(r) {
  ff <- file.path(r$dir, paste0("varimp_final_", r$tag, ".csv"))
  cf <- file.path(r$dir, paste0("varimp_cv_",    r$tag, ".csv"))
  if (!file.exists(ff) || !file.exists(cf)) { message("SKIP (missing): ", r$tag); return(NULL) }

  fin <- read.csv(ff) %>%
    dplyr::transmute(Variable, imp_all = importance,
                     scaled_all = imp_all / max(imp_all, na.rm = TRUE))
  cvv <- read.csv(cf) %>%
    dplyr::transmute(Variable, imp_cv = mean_imp,
                     scaled_cv = imp_cv / max(mean_imp, na.rm = TRUE))

  m <- dplyr::full_join(fin, cvv, by = "Variable") %>%
    dplyr::mutate(dplyr::across(c(imp_all, scaled_all, imp_cv, scaled_cv),
                                ~ tidyr::replace_na(., 0)),
                  important_all = scaled_all >= IMP_CUTOFF,
                  important_cv  = scaled_cv  >= IMP_CUTOFF,
                  status = dplyr::case_when(
                    important_all &  important_cv ~ "both",
                    important_all & !important_cv ~ "all87_only",
                   !important_all &  important_cv ~ "cv_only",
                    TRUE                          ~ "neither"),
                  modality = r$modality, model = r$model,
                  criterion = sub("_", "", r$crit), tag = r$tag)

  # per-run merged table + Venn
  write.csv(m %>% dplyr::arrange(dplyr::desc(scaled_all)),
            file.path(OUT, paste0("compare_", r$tag, ".csv")), row.names = FALSE)
  n_both <- sum(m$status == "both"); n_all <- sum(m$status == "all87_only"); n_cv <- sum(m$status == "cv_only")
  venn2(n_all, n_both, n_cv, "all-87 fit", "CV folds",
        sprintf("%s · %s · %s  (cutoff %.2f)", r$model, r$modality, sub("_", "", r$crit), IMP_CUTOFF),
        file.path(OUT, paste0("venn_", r$tag, ".png")))
  m
}

all_cmp <- purrr::pmap_dfr(runs, function(...) compare_one(tibble::tibble(...)))
if (is.null(all_cmp) || !nrow(all_cmp))
  stop("No varimp_final_/varimp_cv_ files found. Run models_87.R with EXTRACT_CV_IMP = TRUE first.")

## ---- overlap summary across every run ----------------------------------------
overlap <- all_cmp %>%
  dplyr::group_by(modality, model, criterion, tag) %>%
  dplyr::summarise(both = sum(status == "both"),
                   all87_only = sum(status == "all87_only"),
                   cv_only = sum(status == "cv_only"),
                   jaccard = both / pmax(1, both + all87_only + cv_only),
                   .groups = "drop") %>%
  dplyr::arrange(modality, model, criterion)
write.csv(overlap, file.path(OUT, "overlap_summary.csv"), row.names = FALSE)
cat("\nCV-vs-all87 agreement (features important both ways / union), per run:\n")
print(as.data.frame(overlap), row.names = FALSE)

## ---- FINAL lists: robust features (important BOTH ways) ----------------------
# proteins  <- from the proteomics runs (07); symptoms <- from the symptoms runs (06)
final_list <- function(mod_name, out_csv, what) {
  d <- all_cmp %>% dplyr::filter(modality == mod_name, status == "both")
  if (!nrow(d)) { message("no robust ", what, " found"); return(invisible()) }
  out <- d %>% dplyr::group_by(Variable) %>%
    dplyr::summarise(n_runs_robust = dplyr::n(),                 # out of EN/RF x 2 criteria = up to 4
                     runs = paste(sort(paste0(model, ":", criterion)), collapse = ", "),
                     mean_scaled_all = mean(scaled_all),
                     mean_scaled_cv  = mean(scaled_cv), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_runs_robust), dplyr::desc(mean_scaled_all))
  write.csv(out, file.path(OUT, out_csv), row.names = FALSE)
  cat(sprintf("\n#### FINAL %s important in BOTH CV and all-87 (%d features) ####\n", what, nrow(out)))
  print(as.data.frame(utils::head(out, 25)), row.names = FALSE)
}
final_list("proteomics",          "final_important_proteins.csv", "PROTEINS")
final_list("symptoms+background", "final_important_symptoms.csv", "SYMPTOMS/BACKGROUND")

message("\nDone. Per-run Venns + compare tables, overlap_summary.csv, and the two ",
        "final_important_*.csv lists are in ", OUT, "/")
message("Note: 'both' = scaled importance >= ", IMP_CUTOFF,
        " in BOTH the all-87 fit and the fold-averaged CV fit. ",
        "Integrated (08) gets Venns/compare tables too, but the final lists are ",
        "sourced from the single-modality models (07 proteins, 06 symptoms) to keep ",
        "protein-vs-symptom labelling unambiguous.")
