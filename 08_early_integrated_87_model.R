# ============================================================
# 08_early_integrated_87_model.R  — early-integration model on the 87 cohort
# ------------------------------------------------------------
# EARLY fusion: concatenate the selected proteins (from 05, chosen on LCP2) and
# the selected symptom/background variables (from 05, chosen on the 411) into ONE
# feature table for the 87 patients, and train EN + RF on it. The third of the
# three modelling approaches; 09_compare_performance compares this with 06
# (symptoms+bg only) and 07 (proteomics only).
#
# Feature selection happened on the independent cohorts (411, LCP2), never on the
# 87, so the 87 stays a clean testbed. Folds are built inline with the SAME seed /
# v / repeats / patient order as 06 and 07, so all three approaches are scored on
# identical CV splits.
#
# Run with working directory = LC_Data, after 02, 03, 04, 05.
# ============================================================

required <- c("tidymodels", "glmnet", "ranger", "future")
to_install <- required[!required %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)
library(tidymodels); library(glmnet); library(ranger); library(future)

# CV folds are built inline below (no folds87_utils.R, no .rds file)
plan(multisession, workers = 18)
set.seed(900)

OUT    <- "results/08_early_integrated_87"
SUFFIX <- ""                    # "_loose" for the sensitivity variable set
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

sanitize <- function(x) make.names(gsub("[^A-Za-z0-9_]+", "_", x), unique = TRUE)

## ---- build the COMBINED reduced table (proteins + symptoms+bg), aligned ------
# ml87_proteomics.csv sets the canonical Patient order that the folds use; the
# symptom table is re-ordered to match, then the two are joined on Patient.
prot87 <- read.csv("processed_data/ml87_proteomics.csv", check.names = FALSE)
canon  <- prot87$Patient

top_p  <- read.csv(file.path("results/task3_variable_selection",
                             paste0("top_variables_proteins", SUFFIX, ".csv")))$Variable
keep_p <- intersect(top_p, names(prot87))
if (length(setdiff(top_p, names(prot87))))
  warning("proteins not found in 87 panel: ", paste(setdiff(top_p, names(prot87)), collapse = ", "))
prot_part <- prot87 %>% dplyr::select(Patient, Lungcancer, dplyr::all_of(keep_p))   # already canon order

symp87 <- read.csv("processed_data/symptoms_bg_test87.csv", check.names = FALSE)
names(symp87) <- sanitize(names(symp87))
top_s  <- read.csv(file.path("results/task3_variable_selection",
                             paste0("top_variables_symptoms_background", SUFFIX, ".csv")))$Variable
keep_s <- intersect(top_s, names(symp87))
if (length(setdiff(top_s, names(symp87))))
  warning("symptom vars not found: ", paste(setdiff(top_s, names(symp87)), collapse = ", "))
symp_part <- symp87 %>%
  dplyr::select(Patient, dplyr::all_of(keep_s)) %>%
  dplyr::slice(match(canon, Patient))

# proteins are prot_*, symptoms are Q*/Br_*/... -> no name clashes to worry about
combined <- prot_part %>% dplyr::left_join(symp_part, by = "Patient")
stopifnot(nrow(combined) == length(canon),
          identical(as.character(combined$Patient), as.character(canon)),
          !anyNA(combined %>% dplyr::select(-Patient, -Lungcancer)))
write.csv(combined, file.path("processed_data", paste0("ml87_top_combined", SUFFIX, ".csv")),
          row.names = FALSE)
cat(sprintf("combined 87 table: %d patients x %d predictors (%d proteins + %d symptom/bg)\n",
            nrow(combined), length(keep_p) + length(keep_s), length(keep_p), length(keep_s)))

## ---- fit EN + RF on shared inline folds --------------------------------------
en_metrics <- metric_set(roc_auc, sens, yardstick::spec, bal_accuracy, j_index)

train <- combined %>% dplyr::select(-Patient) %>%
  mutate(Lungcancer = factor(Lungcancer, levels = c("Yes", "No")))   # positive first

# CV folds: built ONCE here, shared by EN and RF (no .rds, no retarget).
# Same 87 patients/order/label as 06/07 + same seed => identical folds, so 09 can
# compare the three approaches fairly.
set.seed(900)
folds <- vfold_cv(train, v = 5, repeats = 10, strata = Lungcancer)

run_one <- function(train, model, label) {
  p <- ncol(train) - 1
  # `folds` is the shared object built once above (identical for EN and RF)
  if (model == "EN") {
    rec  <- recipe(Lungcancer ~ ., data = train) %>%
      step_zv(all_predictors()) %>% step_normalize(any_of("Age"))
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
  best <- select_best(res, metric = "roc_auc")
  cv   <- collect_metrics(res) %>%
    semi_join(best, by = intersect(names(best), names(collect_metrics(res)))) %>%
    dplyr::select(.metric, mean, std_err, n) %>% dplyr::mutate(model_set = label)
  oof <- collect_predictions(res, parameters = best) %>%
    dplyr::group_by(.row) %>%
    dplyr::summarise(.pred_Yes = mean(.pred_Yes),
                     Lungcancer = dplyr::first(Lungcancer), .groups = "drop") %>%
    dplyr::mutate(model_set = label)
  cat(sprintf("%s done | CV ROC-AUC = %.3f\n", label, cv$mean[cv$.metric == "roc_auc"]))
  list(cv = cv, oof = oof)
}

res_en <- run_one(train, "EN", "EN · integrated")
res_rf <- run_one(train, "RF", "RF · integrated")

cv_all  <- dplyr::bind_rows(res_en$cv,  res_rf$cv)
oof_all <- dplyr::bind_rows(res_en$oof, res_rf$oof)
write.csv(cv_all,  file.path(OUT, paste0("cv_performance", SUFFIX, ".csv")), row.names = FALSE)
write.csv(oof_all, file.path(OUT, paste0("oof_predictions", SUFFIX, ".csv")), row.names = FALSE)
cat("\nCV ROC-AUC (early-integrated, 87):\n")
print(as.data.frame(cv_all %>% dplyr::filter(.metric == "roc_auc") %>%
                      dplyr::select(model_set, mean, std_err)), row.names = FALSE)

plan(sequential)
message("\nDone. Results in ", OUT, "/")
