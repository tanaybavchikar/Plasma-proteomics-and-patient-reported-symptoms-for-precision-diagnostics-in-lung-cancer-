# ============================================================
# models_87.R  — one function for all three 87-cohort modelling approaches
# ------------------------------------------------------------
# A reusable alternative to running 06 + 07 + 08 one by one. `run_modality()`
# takes a vector of selected proteins and a vector of selected symptom/background
# variables, builds the corresponding 87-patient table, and trains EN + RF on it:
#   proteins only            -> pass symp_feats = character(0)   (== 07)
#   symptoms+background only  -> pass prot_feats = character(0)   (== 06)
#   both (early integration)  -> pass both                        (== 08)
#
# It writes to the SAME results folders as 06/07/08, so 09_compare_performance.R
# reads its output unchanged. Folds are built inline (v=5, repeats=10, seed 900),
# identical across the three approaches, so the comparison stays fair.
#
# Run with working directory = LC_Data, after 02 and 05.
# ============================================================

required <- c("tidymodels", "glmnet", "ranger", "future", "vip")
to_install <- required[!required %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)
library(tidymodels); library(glmnet); library(ranger); library(future); library(vip)

source("../code/helpers/varimp_across_folds.R")   # fold_importance / summarise_importance

plan(multisession, workers = 18)
SUFFIX <- ""              # "_loose" for the sensitivity variable set
EXTRACT_CV_IMP <- TRUE    # also compute fold-averaged importance (for CV-vs-all-data comparison;
                          # set FALSE to skip the extra ~100 fits/model)

sanitize <- function(x) make.names(gsub("[^A-Za-z0-9_]+", "_", x), unique = TRUE)

## ---- the reusable function ---------------------------------------------------
# prot_feats, symp_feats : character vectors of feature names (either may be empty)
# modality_label         : e.g. "proteomics" / "symptoms+background" / "integrated"
# out_dir                : where to write cv_performance + oof (e.g. results/07_proteomics_87)
run_modality <- function(prot_feats, symp_feats, modality_label, out_dir, suffix = SUFFIX) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  ## -- build the 87 table from the requested features (aligned to canon order) --
  prot87 <- read.csv("processed_data/ml87_proteomics.csv", check.names = FALSE)
  canon  <- prot87$Patient
  tab    <- prot87 %>% dplyr::select(Patient, Lungcancer)     # id + label, canon order

  if (length(prot_feats)) {
    keep_p <- intersect(prot_feats, names(prot87))
    if (length(setdiff(prot_feats, names(prot87))))
      warning("proteins not found: ", paste(setdiff(prot_feats, names(prot87)), collapse = ", "))
    tab <- tab %>% dplyr::left_join(prot87 %>% dplyr::select(Patient, dplyr::all_of(keep_p)), by = "Patient")
  }
  if (length(symp_feats)) {
    symp87 <- read.csv("processed_data/symptoms_bg_test87.csv", check.names = FALSE)
    names(symp87) <- sanitize(names(symp87))
    keep_s <- intersect(symp_feats, names(symp87))
    if (length(setdiff(symp_feats, names(symp87))))
      warning("symptom vars not found: ", paste(setdiff(symp_feats, names(symp87)), collapse = ", "))
    tab <- tab %>% dplyr::left_join(symp87 %>% dplyr::select(Patient, dplyr::all_of(keep_s)), by = "Patient")
  }
  stopifnot(identical(as.character(tab$Patient), as.character(canon)),
            !anyNA(tab %>% dplyr::select(-Patient, -Lungcancer)))

  key <- gsub("[^A-Za-z0-9]+", "_", modality_label)
  write.csv(tab, file.path("processed_data", paste0("ml87_sel_", key, suffix, ".csv")), row.names = FALSE)
  cat(sprintf("\n[%s] %d patients x %d predictors\n", modality_label, nrow(tab), ncol(tab) - 2))

  ## -- CV folds built ONCE, shared by EN and RF (same seed/order as every modality) --
  train <- tab %>% dplyr::select(-Patient) %>%
    mutate(Lungcancer = factor(Lungcancer, levels = c("Yes", "No")))
  set.seed(900)
  folds <- vfold_cv(train, v = 5, repeats = 10, strata = Lungcancer)
  en_metrics <- metric_set(roc_auc, sens, yardstick::spec, bal_accuracy, j_index)

  fit_model <- function(model) {
    p <- ncol(train) - 1
    if (model == "EN") {
      rec  <- recipe(Lungcancer ~ ., data = train) %>%
        step_zv(all_predictors()) %>% step_normalize(any_of("Age"))   # no-op if no Age
      spec <- logistic_reg(penalty = tune(), mixture = tune()) %>%
        set_engine("glmnet") %>% set_mode("classification")
      grid <- expand_grid(mixture = seq(0, 1, by = 0.05),
                          penalty = 10^seq(-4, 0.5, length.out = 40))
    } else {
      rec  <- recipe(Lungcancer ~ ., data = train) %>% step_zv(all_predictors())
      spec <- rand_forest(mtry = tune(), min_n = tune(), trees = tune()) %>%
        set_engine("ranger", importance = "permutation", num.threads = 1) %>%
        set_mode("classification")
      mtry_hi <- max(2L, min(p, as.integer(ceiling(sqrt(p) * 3))))
      grid <- grid_regular(mtry(range = c(2L, mtry_hi)), min_n(range = c(2L, 20L)),
                           trees(range = c(500L, 2000L)), levels = 5)
    }
    wf  <- workflow() %>% add_recipe(rec) %>% add_model(spec)
    res <- tune_grid(wf, resamples = folds, grid = grid, metrics = en_metrics,
                     control = control_grid(save_pred = TRUE, parallel_over = "everything"))
    best  <- select_best(res, metric = "roc_auc")
    label <- paste(model, "·", modality_label)   # e.g. "EN · proteomics"

    ## --- FINAL model on ALL samples (best hyperparameters) + its importance ----
    # NOTE: no honest performance from this fit (resubstitution); the honest
    # performance is the CV out-of-fold AUC below. This is for the IMPORTANCE
    # comparison: all-data fit (varimp_final_*) vs fold-averaged (varimp_cv_*).
    wf_final  <- finalize_workflow(wf, best)
    final_fit <- wf_final %>% fit(data = train)
    if (model == "EN") {
      extract_fn <- function(fit, split) {
        g <- extract_fit_engine(fit); co <- coef(g, s = best$penalty)
        tibble(Variable = rownames(co), importance = abs(as.numeric(co))) %>%
          dplyr::filter(Variable != "(Intercept)")
      }
    } else {
      extract_fn <- function(fit, split)
        extract_fit_parsnip(fit) %>% vip::vi() %>% dplyr::transmute(Variable, importance = Importance)
    }
    tag <- paste0(gsub("[^A-Za-z0-9]+", "_", modality_label), "_", model, suffix)
    extract_fn(final_fit, NULL) %>% dplyr::arrange(dplyr::desc(importance)) %>%
      write.csv(file.path(out_dir, paste0("varimp_final_", tag, ".csv")), row.names = FALSE)
    saveRDS(final_fit, file.path(out_dir, paste0("final_model_", tag, ".rds")))
    if (EXTRACT_CV_IMP) {
      fold_importance(wf_final, folds, extract_fn) %>% summarise_importance() %>%
        write.csv(file.path(out_dir, paste0("varimp_cv_", tag, ".csv")), row.names = FALSE)
    }

    cv <- collect_metrics(res) %>%
      semi_join(best, by = intersect(names(best), names(collect_metrics(res)))) %>%
      dplyr::select(.metric, mean, std_err, n) %>% dplyr::mutate(model_set = label)
    oof <- collect_predictions(res, parameters = best) %>%
      dplyr::group_by(.row) %>%
      dplyr::summarise(.pred_Yes = mean(.pred_Yes),
                       Lungcancer = dplyr::first(Lungcancer), .groups = "drop") %>%
      dplyr::mutate(model_set = label)
    cat(sprintf("  %s | CV ROC-AUC = %.3f\n", label, cv$mean[cv$.metric == "roc_auc"]))
    list(cv = cv, oof = oof)
  }

  r_en <- fit_model("EN"); r_rf <- fit_model("RF")
  cv_all  <- dplyr::bind_rows(r_en$cv,  r_rf$cv)
  oof_all <- dplyr::bind_rows(r_en$oof, r_rf$oof)
  write.csv(cv_all,  file.path(out_dir, paste0("cv_performance",  suffix, ".csv")), row.names = FALSE)
  write.csv(oof_all, file.path(out_dir, paste0("oof_predictions", suffix, ".csv")), row.names = FALSE)
  invisible(list(cv = cv_all, oof = oof_all))
}

## ---- driver: run BOTH selection criteria x the three approaches --------------
# Each criterion's feature lists come from 05 (top_variables_*_impfrac.csv /
# _selfreq.csv). The criterion tag is used as the file SUFFIX throughout, so each
# writes cv_performance<crit>.csv / oof_predictions<crit>.csv into the same
# results/06_*, 07_*, 08_* folders -> run 09 once per criterion (SUFFIX = crit).
for (crit in c("_impfrac", "_selfreq")) {
  message("\n########## selection criterion: ", crit, " ##########")
  top_p <- read.csv(file.path("results/task3_variable_selection",
                              paste0("top_variables_proteins", crit, ".csv")))$Variable
  top_s <- read.csv(file.path("results/task3_variable_selection",
                              paste0("top_variables_symptoms_background", crit, ".csv")))$Variable

  run_modality(character(0), top_s,        "symptoms+background", "results/06_symptombg_87",        suffix = crit)  # == 06
  run_modality(top_p,        character(0), "proteomics",          "results/07_proteomics_87",       suffix = crit)  # == 07
  run_modality(top_p,        top_s,        "integrated",          "results/08_early_integrated_87", suffix = crit)  # == 08
}

plan(sequential)
message("\nDone. Both criteria x three approaches trained; results in results/06_*, 07_*, 08_* ",
        "(files tagged _impfrac / _selfreq). Run 09 with SUFFIX = '_impfrac' then '_selfreq'.")
