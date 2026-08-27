# ============================================================
# 06_symptombg_87_model.R  — symptoms+background-only model on the 87 cohort
# ------------------------------------------------------------
# Trains EN + RF on the 87 patients using ONLY the selected symptom/background
# variables (from 05_select_features). One of the three modelling approaches;
# 09_compare_performance will read this + 07 (+ 08) to compare.
#
# Selection was done on the 411 cohort (05), so the 87 stays a clean testbed.
# Uses the shared cv_folds87 folds via retarget_folds().
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

OUT    <- "results/06_symptombg_87"
SUFFIX <- ""                    # "_loose" for the sensitivity variable set
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

sanitize <- function(x) make.names(gsub("[^A-Za-z0-9_]+", "_", x), unique = TRUE)

## ---- build the reduced symptoms+background table (aligned to cv_folds87) -----
# ml87_proteomics.csv defines the canonical Patient order that cv_folds87 was
# built on; the symptom table must match it row-for-row.
canon  <- read.csv("processed_data/ml87_proteomics.csv", check.names = FALSE)$Patient
symp87 <- read.csv("processed_data/symptoms_bg_test87.csv", check.names = FALSE)
names(symp87) <- sanitize(names(symp87))
top_s  <- read.csv(file.path("results/task3_variable_selection",
                             paste0("top_variables_symptoms_background", SUFFIX, ".csv")))$Variable
keep_s <- intersect(top_s, names(symp87))
if (length(setdiff(top_s, names(symp87))))
  warning("selected symptom vars not found: ", paste(setdiff(top_s, names(symp87)), collapse = ", "))

reduced <- symp87 %>%
  dplyr::select(Patient, Lungcancer, dplyr::all_of(keep_s)) %>%
  dplyr::slice(match(canon, Patient))
stopifnot(nrow(reduced) == length(canon),
          identical(as.character(reduced$Patient), as.character(canon)))
write.csv(reduced, file.path("processed_data", paste0("ml87_top_symptbg", SUFFIX, ".csv")),
          row.names = FALSE)
cat(sprintf("symptoms+background 87 table: %d patients x %d predictors\n",
            nrow(reduced), length(keep_s)))

## ---- fit EN + RF on the shared folds -----------------------------------------
en_metrics <- metric_set(roc_auc, sens, yardstick::spec, bal_accuracy, j_index)

train <- reduced %>% dplyr::select(-Patient) %>%
  mutate(Lungcancer = factor(Lungcancer, levels = c("Yes", "No")))   # positive first

# CV folds: built ONCE here, shared by EN and RF (no .rds, no retarget).
# Same 87 patients/order/label as 07/08 + same seed => identical folds across the
# three 87-cohort approaches, so 09 can compare them fairly.
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

res_en <- run_one(train, "EN", "EN · symptoms+background")
res_rf <- run_one(train, "RF", "RF · symptoms+background")

cv_all  <- dplyr::bind_rows(res_en$cv,  res_rf$cv)
oof_all <- dplyr::bind_rows(res_en$oof, res_rf$oof)
write.csv(cv_all,  file.path(OUT, paste0("cv_performance", SUFFIX, ".csv")), row.names = FALSE)
write.csv(oof_all, file.path(OUT, paste0("oof_predictions", SUFFIX, ".csv")), row.names = FALSE)
cat("\nCV ROC-AUC (symptoms+background, 87):\n")
print(as.data.frame(cv_all %>% dplyr::filter(.metric == "roc_auc") %>%
                      dplyr::select(model_set, mean, std_err)), row.names = FALSE)

plan(sequential)
message("\nDone. Results in ", OUT, "/")
