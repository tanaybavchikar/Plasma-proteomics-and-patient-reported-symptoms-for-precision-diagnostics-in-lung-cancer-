# ============================================================
# 03_symptombg_411_model.R  — EN + RF on symptoms+background (411), for feature selection
# merged: en_symptbg.R + rf_symptbg.R
# Run with working directory = LC_Data.
# NOTE: faithful merge of the scripts named above; only data paths repointed
# to raw_data/ and metadata/. Preserved originals are in code/archive/.
# ============================================================

# =============================================================================
# en_symptbg.R  (Task 2: elastic net on SYMPTOMS + BACKGROUND)
# -----------------------------------------------------------------------------
# Same model/grid/metrics as ml_en_tidymodels.R (the symptoms-only EN), but on
# symptoms_bg_train.csv (110 symptoms + 24 background factors). Reuses the SHARED
# cv_folds.rds row partitions via retarget_folds() so this is directly
# comparable to the symptoms-only EN on identical CV splits.
#
# Differences from the symptoms-only script (all intentional):
#   - data = symptoms_bg_train.csv
#   - folds retargeted onto this table (cv_folds.rds was built on symptoms_train,
#     same 411 patients / same row order, so the partitions carry over exactly)
#   - step_normalize(Age) now actually does something (Age is a real predictor)
#   - glmnet SIGN FIX carried over: coefficients negated to read as log-odds of
#     cancer (glmnet parameterises toward the 2nd factor level = "No").
#
# Run with working dir = LC_Data, after build_ml_tables_symptbg.R (+ make_folds.R
# for cv_folds.rds).
# =============================================================================

## ---- Packages ---------------------------------------------------------------
required <- c("tidymodels", "glmnet", "readxl", "future")
to_install <- required[!required %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)
library(tidymodels); library(glmnet); library(readxl); library(future)

plan(multisession, workers = 18)   # server cores (PROJECT_STATUS.md) — parallel only affects speed

RES <- "results/en_symptbg"
if (!dir.exists(RES)) dir.create(RES, recursive = TRUE)
set.seed(900)
source("../code/helpers/varimp_across_folds.R")   # fold_importance / summarise_importance / plot_importance
# (folds built inline below — retarget_folds/folds87_utils no longer used)

## ---- Data -------------------------------------------------------------------
train_full <- read.csv("processed_data/symptoms_bg_train.csv", check.names = FALSE)
names(train_full) <- make.names(gsub("[^A-Za-z0-9_]+", "_", names(train_full)), unique = TRUE)
train <- train_full %>%
  dplyr::select(-Patient) %>%
  mutate(Lungcancer = factor(Lungcancer, levels = c("Yes", "No")))   # positive class first
cat("Training rows:", nrow(train), " | predictors:", ncol(train) - 1, "\n")
print(table(train$Lungcancer))

## ---- Recipe -----------------------------------------------------------------
en_recipe <- recipe(Lungcancer ~ ., data = train) %>%
  step_zv(all_predictors()) %>%
  step_normalize(any_of("Age"))   # Age now present -> real z-scoring; symptoms are 0/1 (glmnet standardises too)

## ---- Model: elastic-net logistic regression ---------------------------------
en_spec <- logistic_reg(penalty = tune(), mixture = tune()) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

en_wf <- workflow() %>% add_recipe(en_recipe) %>% add_model(en_spec)

## ---- Tuning grid (identical to the symptoms-only EN) ------------------------
en_grid <- expand_grid(mixture = seq(0, 1, by = 0.05),               # 21 alphas
                       penalty = 10^seq(-4, 0.5, length.out = 40))   # 40 lambdas

## ---- CV folds: built ONCE here, shared by EN and RF (no .rds, no retarget) ---
# Fixed seed => the RF section further down reuses this exact `folds` object.
set.seed(900)
folds <- vfold_cv(train, v = 10, repeats = 5, strata = Lungcancer)

## ---- Metrics ----------------------------------------------------------------
en_metrics <- metric_set(roc_auc, sens, yardstick::spec, bal_accuracy, j_index)

## ---- Tune -------------------------------------------------------------------
en_res <- tune_grid(
  en_wf, resamples = folds, grid = en_grid, metrics = en_metrics,
  control = control_grid(save_pred = TRUE, parallel_over = "everything"))

cv_metrics <- collect_metrics(en_res)
write.csv(cv_metrics, file.path(RES, "cv_metrics_all.csv"), row.names = FALSE)
cat("\nTop 5 by CV ROC-AUC:\n"); print(show_best(en_res, metric = "roc_auc", n = 5))

autoplot(en_res, metric = "roc_auc") +
  labs(title = "EN tuning (symptoms+background): CV ROC-AUC across penalty and mixture")
ggsave(file.path(RES, "en_tuning_plot.png"), width = 8, height = 5, dpi = 150)

## ---- Best model -------------------------------------------------------------
best_en  <- select_best(en_res, metric = "roc_auc")
cat("\nBest hyperparameters:\n"); print(best_en)
final_en <- finalize_workflow(en_wf, best_en) %>% fit(data = train)
saveRDS(final_en, file.path(RES, "final_en_model.rds"))

best_cv <- cv_metrics %>%
  semi_join(best_en, by = c("penalty", "mixture")) %>%
  dplyr::select(.metric, mean, std_err, n)
write.csv(best_cv, file.path(RES, "cv_performance_best.csv"), row.names = FALSE)
cat("\nCV performance of the best EN model:\n"); print(best_cv)

## ---- (A) Youden-optimal threshold performance -------------------------------
oof <- collect_predictions(en_res, parameters = best_en)
roc_pts <- roc_curve(oof, truth = Lungcancer, .pred_Yes) %>%
  mutate(youden = sensitivity + specificity - 1)
best_thr <- roc_pts$.threshold[which.max(roc_pts$youden)]
oof <- oof %>%
  mutate(pred_cls = factor(ifelse(.pred_Yes >= best_thr, "Yes", "No"), levels = c("Yes", "No")))
youden_perf <- tibble(
  threshold    = best_thr,
  sens         = yardstick::sens_vec(oof$Lungcancer, oof$pred_cls),
  spec         = yardstick::spec_vec(oof$Lungcancer, oof$pred_cls),
  bal_accuracy = yardstick::bal_accuracy_vec(oof$Lungcancer, oof$pred_cls))
write.csv(youden_perf, file.path(RES, "performance_youden_threshold.csv"), row.names = FALSE)
cat("\nPerformance at Youden-optimal threshold:\n"); print(youden_perf)

## ---- (B) Selected features: non-zero EN coefficients ------------------------
glmnet_fit <- extract_fit_engine(final_en)
co <- coef(glmnet_fit, s = best_en$penalty)
# SIGN FIX (see ml_en_tidymodels.R): glmnet targets the 2nd level ("No"), so raw
# coefficients are log-odds of CONTROL. Negate to express as cancer log-odds.
coef_df <- tibble(term = rownames(co), coefficient = -as.numeric(co)) %>%
  filter(term != "(Intercept)", coefficient != 0) %>%
  mutate(odds_ratio = exp(coefficient),
         direction  = ifelse(coefficient > 0, "higher cancer risk", "lower cancer risk")) %>%
  arrange(desc(abs(coefficient)))

dict113 <- read_excel("metadata/113 Merged file_cleanced.xlsx",
                      sheet = "Variable_explaination", col_names = c("var", "meaning")) %>%
  mutate(term = make.names(gsub("[^A-Za-z0-9_]+", "_", var)))
coef_df <- coef_df %>%
  left_join(dplyr::select(dict113, term, meaning), by = "term") %>%
  relocate(meaning, .after = term)
write.csv(coef_df, file.path(RES, "en_selected_features.csv"), row.names = FALSE)
cat("\nSelected features (non-zero coefficients):", nrow(coef_df), "\n")
print(head(coef_df, 15))

coef_df %>% slice_max(abs(coefficient), n = 20) %>%
  ggplot(aes(reorder(coalesce(meaning, term), coefficient), coefficient,
             fill = coefficient > 0)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("steelblue", "firebrick"), guide = "none") +
  labs(title = "EN (symptoms+background): top selected features (coef = log-odds of cancer)",
       x = NULL, y = "coefficient") + theme_minimal()
ggsave(file.path(RES, "en_selected_features.png"), width = 8, height = 6, dpi = 150)

## ---- (C) Label-shuffle sanity check -----------------------------------------
set.seed(1)
train_shuf <- train %>% mutate(Lungcancer = sample(Lungcancer))
shuf_res <- fit_resamples(
  finalize_workflow(en_wf, best_en),
  resamples = vfold_cv(train_shuf, v = 5, strata = Lungcancer),
  metrics = metric_set(roc_auc))
auc_shuf <- collect_metrics(shuf_res) %>% filter(.metric == "roc_auc") %>% pull(mean)
cat(sprintf("\nLabel-shuffle sanity check: CV ROC-AUC = %.3f (should be ~0.5)\n", auc_shuf))

## ---- (D) Variable importance across ALL folds -------------------------------
en_wf_final <- finalize_workflow(en_wf, best_en)
extract_en <- function(fit, split) {
  g  <- extract_fit_engine(fit)
  co <- coef(g, s = best_en$penalty)
  tibble(Variable = rownames(co), importance = abs(as.numeric(co))) %>%
    filter(Variable != "(Intercept)")
}
en_imp_sum <- fold_importance(en_wf_final, folds, extract_en) %>%
  summarise_importance() %>%
  left_join(dplyr::select(dict113, Variable = term, meaning), by = "Variable") %>%
  mutate(label = coalesce(meaning, Variable))
write.csv(en_imp_sum, file.path(RES, "varimp_across_folds.csv"), row.names = FALSE)

ggsave(file.path(RES, "varimp_mean.png"),
       plot_importance(en_imp_sum, "mean_imp",
                       "EN (symptoms+background): mean |coef| across folds", label_col = "label"),
       width = 8, height = 7, dpi = 150)
ggsave(file.path(RES, "varimp_selfreq.png"),
       plot_importance(en_imp_sum, "selection_freq",
                       "EN (symptoms+background): selection frequency across folds", label_col = "label"),
       width = 8, height = 7, dpi = 150)

plan(sequential)
message("\nDone. EN (symptoms+background) results saved in ", RES, "/")


# ================= rf_symptbg.R =================

# =============================================================================
# rf_symptbg.R  (Task 2: random forest on SYMPTOMS + BACKGROUND)
# -----------------------------------------------------------------------------
# Same model/grid/metrics as RF_tinymodels.R (the symptoms-only RF), but on
# symptoms_bg_train.csv (110 symptoms + 24 background factors), reusing the
# SHARED cv_folds.rds partitions via retarget_folds() so it's directly
# comparable to the symptoms-only RF on identical CV splits.
#
# Run with working dir = LC_Data, after build_ml_tables_symptbg.R (+ make_folds.R).
# =============================================================================

#----- PACKAGES -----#
required <- c("tidymodels", "ranger", "vip", "readxl", "future")
to_install <- required[!required %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)
library(tidymodels); library(ranger); library(vip); library(readxl); library(future)

plan(multisession, workers = 18)   # server cores (PROJECT_STATUS.md)

RES <- "results/rf_symptbg"
if (!dir.exists(RES)) dir.create(RES, recursive = TRUE)
set.seed(900)
source("../code/helpers/varimp_across_folds.R")
# (folds built inline below — retarget_folds/folds87_utils no longer used)

#----- DATA -----#
train_full <- read.csv("processed_data/symptoms_bg_train.csv", check.names = FALSE)
names(train_full) <- make.names(gsub("[^A-Za-z0-9_]+", "_", names(train_full)), unique = TRUE)
train <- train_full %>% dplyr::select(-Patient) %>%
  mutate(Lungcancer = factor(Lungcancer, levels = c("Yes", "No")))   # positive class first
cat("Training rows:", nrow(train), " | predictors:", ncol(train) - 1, "\n")
print(table(train$Lungcancer))

#----- RECIPE / MODEL / WORKFLOW -----#
rf_recipe <- recipe(Lungcancer ~ ., data = train) %>%
  step_zv(all_predictors())                       # trees need no scaling

rf_spec <- rand_forest(mtry = tune(), min_n = tune(), trees = tune()) %>%
  set_engine("ranger", importance = "permutation", num.threads = 1) %>%
  set_mode("classification")

rf_wf <- workflow() %>% add_recipe(rf_recipe) %>% add_model(rf_spec)

#----- CV folds: reuse the SAME folds built once in the EN section above -----#
# (identical object shared by both models; nothing rebuilt, no file loaded)

#----- GRID + METRICS (same as symptoms-only RF) -----#
rf_grid <- grid_regular(mtry(range = c(2L, 40L)),
                        min_n(range = c(2L, 40L)),
                        trees(range = c(500L, 2000L)),
                        levels = 5)                    # 5^3 = 125 combos
rf_metrics <- metric_set(roc_auc, sens, yardstick::spec, bal_accuracy, j_index)

#----- TUNE -----#
rf_res <- tune_grid(rf_wf, resamples = folds, grid = rf_grid,
                    metrics = rf_metrics,
                    control = control_grid(save_pred = TRUE, parallel_over = "everything"))

#----- RESULTS -----#
cv_metrics <- collect_metrics(rf_res)
write.csv(cv_metrics, file.path(RES, "cv_metrics_all.csv"), row.names = FALSE)
cat("\nTop 5 by CV ROC-AUC:\n"); print(show_best(rf_res, metric = "roc_auc", n = 5))

autoplot(rf_res, metric = "roc_auc") +
  labs(title = "RF tuning (symptoms+background): CV ROC-AUC across mtry and min_n")
ggsave(file.path(RES, "rf_tuning_plot.png"), width = 8, height = 5, dpi = 150)

best_rf  <- select_best(rf_res, metric = "roc_auc")
final_rf <- finalize_workflow(rf_wf, best_rf) %>% fit(data = train)
saveRDS(final_rf, file.path(RES, "final_rf_model.rds"))

best_cv <- cv_metrics %>% semi_join(best_rf, by = c("mtry", "min_n", "trees")) %>%
  dplyr::select(.metric, mean, std_err, n)
write.csv(best_cv, file.path(RES, "cv_performance_best.csv"), row.names = FALSE)
cat("\nCV performance of the best RF:\n"); print(best_cv)

#----- VARIABLE IMPORTANCE ACROSS ALL FOLDS -----#
extract_rf <- function(fit, split) {
  extract_fit_parsnip(fit) %>% vip::vi() %>%
    dplyr::transmute(Variable, importance = Importance)
}
dict113 <- read_excel("metadata/113 Merged file_cleanced.xlsx", sheet = "Variable_explaination",
                      col_names = c("var", "meaning")) %>%
  mutate(term = make.names(gsub("[^A-Za-z0-9_]+", "_", var)))

rf_imp_sum <- fold_importance(finalize_workflow(rf_wf, best_rf), folds, extract_rf) %>%
  summarise_importance() %>%
  left_join(dplyr::select(dict113, Variable = term, meaning), by = "Variable") %>%
  mutate(label = coalesce(meaning, Variable))
write.csv(rf_imp_sum, file.path(RES, "varimp_across_folds.csv"), row.names = FALSE)

ggsave(file.path(RES, "varimp_mean.png"),
       plot_importance(rf_imp_sum, "mean_imp",
                       "RF (symptoms+background): mean importance across folds", label_col = "label"),
       width = 8, height = 7, dpi = 150)
ggsave(file.path(RES, "varimp_selfreq.png"),
       plot_importance(rf_imp_sum, "selection_freq",
                       "RF (symptoms+background): selection frequency across folds", label_col = "label"),
       width = 8, height = 7, dpi = 150)

plan(sequential)
message("\nDone. RF (symptoms+background) results saved in ", RES, "/")
