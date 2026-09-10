# ============================================================
# 04_LCP2_proteomics_model.R  — EN + RF on LCP2 proteins, for feature selection
# merged: en_LCP2_proteomics.R + rf_LCP2_proteomics.R + analyse_topprot_sex_vs_cancer.R
# Run with working directory = LC_Data.
# =============================================================================
# en_LCP2_proteomics.R
# -----------------------------------------------------------------------------
# Elastic-net (EN / RLR), PROTEOMICS ONLY, on the LCP2 cohort (~171 patients
# after excluding S_1010, LOD-filtered proteins from clean_LCP2_olink.R).
# Mirrors en_87_proteomics.R exactly in structure, adapted for LCP2's column
# names and label:
#   - data   = ml_LCP2_proteomics.csv (surviving prot_* features, not ~2923)
#   - id col = Sample.ID (not Patient)
#   - label  = LC, levels c("Primary_LC", "No_Cancer") -- positive = Primary_LC
#   - folds  = cv_folds_LCP2.rds, built directly on this table (no retargeting
#     needed yet -- LCP2 has only one modelling table so far)
#   - feature dictionary = prot_name_lookup_LCP2.csv (Variable <-> Assay)
#
# EN's grid does NOT depend on the number of predictors, so it's reused as-is
# from the 87-cohort script (mixture/penalty grid is p-independent). Only RF's
# mtry needs resizing for the new protein count -- see rf_LCP2_proteomics.R.
#
# Run with working dir = LC_Data, after build_ml_tables_LCP2.R + make_folds_LCP2.R.
# =============================================================================

## ---- Packages ---------------------------------------------------------------
required <- c("tidymodels", "glmnet", "future")
to_install <- required[!required %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)
library(tidymodels); library(glmnet); library(future)

# ---- Parallel backend --------------------------------------------------------
# Each individual glmnet fit is cheap (whole penalty path in one call per
# mixture value), but tune_grid still has to repeat that across 21 mixtures x
# 100 resamples (10-fold x 10-repeat) = 2100 fits. Parallelising across
# resamples uses the server's idle cores instead of running that sequentially.
plan(multisession, workers = 18)   # server cores (PROJECT_STATUS.md) — NOT multicore

RES <- "results/en_LCP2_proteomics"
if (!dir.exists(RES)) dir.create(RES, recursive = TRUE)
set.seed(900)
source("../code/helpers/varimp_across_folds.R")   # fold_importance / summarise_importance / plot_importance

## ---- Data -------------------------------------------------------------------
train_full <- read.csv("processed_data/ml_LCP2_proteomics.csv", check.names = FALSE)
train <- train_full %>%
  dplyr::select(-Sample.ID, -Patient.ID) %>%
  mutate(LC = factor(LC, levels = c("Primary_LC", "No_Cancer")))   # cancer = positive, first

cat("Training rows:", nrow(train), " | predictors:", ncol(train) - 1, "\n")
print(table(train$LC))

## ---- Recipe -------------------------------------------------------------------
en_recipe <- recipe(LC ~ ., data = train) %>%
  step_zv(all_predictors()) %>%
  step_normalize(any_of("Age"))   # no-op here; kept for consistency with the other scripts

## ---- Model: elastic-net logistic regression ---------------------------------
en_spec <- logistic_reg(penalty = tune(), mixture = tune()) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

en_wf <- workflow() %>% add_recipe(en_recipe) %>% add_model(en_spec)

## ---- Tuning grid (p-independent -- reused as-is) -----------------------------
en_grid <- expand_grid(mixture = seq(0, 1, by = 0.05),               # 21 alphas
                       penalty = 10^seq(-4, 0.5, length.out = 40))   # 40 lambdas

## ---- CV folds: built ONCE here, shared by EN and RF (no .rds) -----------------
# Fixed seed => the RF section further down reuses this exact `folds` object.
set.seed(900)
folds <- vfold_cv(train, v = 10, repeats = 10, strata = LC)

## ---- Metrics ------------------------------------------------------------------
en_metrics <- metric_set(roc_auc, sens, yardstick::spec, bal_accuracy, j_index)

## ---- Tune over the grid by CV --------------------------------------------------
en_res <- tune_grid(
  en_wf,
  resamples = folds,
  grid      = en_grid,
  metrics   = en_metrics,
  control   = control_grid(save_pred = TRUE, parallel_over = "everything"))

## ---- Inspect tuning results ----------------------------------------------------
cv_metrics <- collect_metrics(en_res)
write.csv(cv_metrics, file.path(RES, "cv_metrics_all.csv"), row.names = FALSE)

cat("\nTop 5 hyperparameter combinations by CV ROC-AUC:\n")
print(show_best(en_res, metric = "roc_auc", n = 5))

autoplot(en_res, metric = "roc_auc") +
  labs(title = "EN tuning (LCP2 proteomics): CV ROC-AUC across penalty and mixture")
ggsave(file.path(RES, "en_tuning_plot.png"), width = 8, height = 5, dpi = 150)

## ---- Pick the best, finalize, and fit on the full LCP2 cohort -------------------
best_en  <- select_best(en_res, metric = "roc_auc")
cat("\nBest hyperparameters:\n"); print(best_en)

final_en <- finalize_workflow(en_wf, best_en) %>% fit(data = train)
saveRDS(final_en, file.path(RES, "final_en_model.rds"))

best_cv <- cv_metrics %>%
  semi_join(best_en, by = intersect(names(best_en), names(cv_metrics))) %>%
  dplyr::select(.metric, mean, std_err, n)
write.csv(best_cv, file.path(RES, "cv_performance_best.csv"), row.names = FALSE)
cat("\nCross-validated performance of the best EN model:\n"); print(best_cv)

## ---- (A) Performance at the Youden-optimal threshold ---------------------------
oof <- collect_predictions(en_res, parameters = best_en)
roc_pts <- roc_curve(oof, truth = LC, .pred_Primary_LC) %>%
  mutate(youden = sensitivity + specificity - 1)
best_thr <- roc_pts$.threshold[which.max(roc_pts$youden)]

oof <- oof %>%
  mutate(pred_cls = factor(ifelse(.pred_Primary_LC >= best_thr, "Primary_LC", "No_Cancer"),
                           levels = c("Primary_LC", "No_Cancer")))
youden_perf <- tibble(
  threshold    = best_thr,
  sens         = yardstick::sens_vec(oof$LC, oof$pred_cls),
  spec         = yardstick::spec_vec(oof$LC, oof$pred_cls),
  bal_accuracy = yardstick::bal_accuracy_vec(oof$LC, oof$pred_cls))
write.csv(youden_perf, file.path(RES, "performance_youden_threshold.csv"), row.names = FALSE)
cat("\nPerformance at Youden-optimal threshold:\n"); print(youden_perf)

## ---- (B) Out-of-fold predicted probability per sample (for later comparison) --
sample_ids <- train_full$Sample.ID
oof_sample <- collect_predictions(en_res, parameters = best_en) %>%
  dplyr::group_by(.row) %>%
  dplyr::summarise(prob_LC = mean(.pred_Primary_LC),
                   obs = dplyr::first(LC), .groups = "drop") %>%
  dplyr::arrange(.row) %>%
  dplyr::mutate(Sample.ID = sample_ids[.row], model = "EN") %>%
  dplyr::select(Sample.ID, model, obs, prob_LC)
write.csv(oof_sample, file.path(RES, "oof_predictions.csv"), row.names = FALSE)

## ---- (C) Selected proteins: non-zero EN coefficients ---------------------------
glmnet_fit <- extract_fit_engine(final_en)
co <- coef(glmnet_fit, s = best_en$penalty)
# SIGN CONVENTION: glmnet parameterises toward the SECOND factor level (here
# "No_Cancer"), so its raw coefficients are log-odds of CONTROL, not cancer.
# We flip the sign so a positive coefficient means higher cancer risk, matching
# the "positive class first" convention every metric above already uses.
# (Exact, not a hack: for a 2-class logit, log-odds(cancer) = -log-odds(control).)
coef_df <- tibble(term = rownames(co), coefficient = -as.numeric(co)) %>%
  filter(term != "(Intercept)", coefficient != 0) %>%
  mutate(odds_ratio = exp(coefficient),
         direction  = ifelse(coefficient > 0, "higher cancer risk", "lower cancer risk")) %>%
  arrange(desc(abs(coefficient)))

# Attach Olink Assay names (prot_name_lookup_LCP2.csv)
prot_lookup <- read.csv("processed_data/prot_name_lookup_LCP2.csv")
coef_df <- coef_df %>%
  left_join(prot_lookup, by = c("term" = "Variable")) %>%
  relocate(Assay, .after = term)
write.csv(coef_df, file.path(RES, "en_selected_features.csv"), row.names = FALSE)
cat("\nSelected proteins (non-zero coefficients):", nrow(coef_df), "\n")
print(head(coef_df, 15))

coef_df %>% slice_max(abs(coefficient), n = 20) %>%
  ggplot(aes(reorder(coalesce(Assay, term), coefficient), coefficient,
             fill = coefficient > 0)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("steelblue", "firebrick"), guide = "none") +
  labs(title = "EN (LCP2 proteomics): top selected proteins (coefficient = log-odds of cancer)",
       x = NULL, y = "coefficient") + theme_minimal()
ggsave(file.path(RES, "en_selected_features.png"), width = 8, height = 6, dpi = 150)

## ---- (D) Leakage sanity check: shuffle labels -> AUC should be ~0.5 ------------
set.seed(1)
train_shuf <- train %>% mutate(LC = sample(LC))
shuf_res <- fit_resamples(
  finalize_workflow(en_wf, best_en),
  resamples = vfold_cv(train_shuf, v = 5, strata = LC),
  metrics = metric_set(roc_auc))
auc_shuf <- collect_metrics(shuf_res) %>% filter(.metric == "roc_auc") %>% pull(mean)
cat(sprintf("\nLabel-shuffle sanity check: CV ROC-AUC = %.3f (should be ~0.5)\n", auc_shuf))

## ---- (E) Variable importance across ALL CV folds -------------------------------
en_wf_final <- finalize_workflow(en_wf, best_en)

extract_en <- function(fit, split) {
  g  <- extract_fit_engine(fit)
  co <- coef(g, s = best_en$penalty)
  tibble(Variable = rownames(co), importance = abs(as.numeric(co))) %>%
    filter(Variable != "(Intercept)")
}

en_imp_long <- fold_importance(en_wf_final, folds, extract_en)
en_imp_sum  <- summarise_importance(en_imp_long) %>%
  left_join(prot_lookup, by = "Variable") %>%
  mutate(label = coalesce(Assay, Variable))
write.csv(en_imp_sum, file.path(RES, "varimp_across_folds.csv"), row.names = FALSE)

ggsave(file.path(RES, "varimp_mean.png"),
       plot_importance(en_imp_sum, "mean_imp",
                       "EN (LCP2 proteomics): mean |coef| across folds", label_col = "label"),
       width = 8, height = 7, dpi = 150)
ggsave(file.path(RES, "varimp_median.png"),
       plot_importance(en_imp_sum, "median_imp",
                       "EN (LCP2 proteomics): median |coef| across folds", label_col = "label"),
       width = 8, height = 7, dpi = 150)
ggsave(file.path(RES, "varimp_selfreq.png"),
       plot_importance(en_imp_sum, "selection_freq",
                       "EN (LCP2 proteomics): selection frequency across folds", label_col = "label"),
       width = 8, height = 7, dpi = 150)
cat("\nVariable importance across folds saved (", nrow(en_imp_sum), " proteins).\n")

plan(sequential)   # release the parallel workers
message("\nDone. EN (LCP2 proteomics) results saved in ", RES, "/")


# ================= rf_LCP2_proteomics.R =================

# =============================================================================
# rf_LCP2_proteomics.R  —  Random forest, PROTEOMICS ONLY, LCP2 cohort
# -----------------------------------------------------------------------------
# Mirrors rf_87_proteomics.R exactly in structure, adapted for LCP2's column
# names/label, with ONE deliberate difference: mtry's range is NOT hardcoded
# to c(10, 200) (that was sized for the 87-cohort's ~2923 proteins). Instead
# it's computed from the actual surviving protein count here, so it stays
# sensible no matter what the LOD filter's final protein count turns out to
# be. See the "mtry range" section below for the reasoning.
#
#   - data   = ml_LCP2_proteomics.csv (surviving prot_* features)
#   - id col = Sample.ID (not Patient)
#   - label  = LC, levels c("Primary_LC", "No_Cancer") -- positive = Primary_LC
#   - folds  = cv_folds_LCP2.rds, built directly on this table
#   - feature dictionary = prot_name_lookup_LCP2.csv (Variable <-> Assay)
#
# Run with working dir = LC_Data, after build_ml_tables_LCP2.R + make_folds_LCP2.R.
# =============================================================================

#----- PACKAGES -----#
required <- c("tidymodels", "ranger", "vip", "future")
to_install <- required[!required %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)
library(tidymodels); library(ranger); library(vip); library(future)

# ---- Parallel backend --------------------------------------------------------
plan(multisession, workers = 18)   # server cores (PROJECT_STATUS.md) — NOT multicore

RES <- "results/rf_LCP2_proteomics"
if (!dir.exists(RES)) dir.create(RES, recursive = TRUE)
set.seed(900)
source("../code/helpers/varimp_across_folds.R")   # fold_importance / summarise_importance / plot_importance

#----- DATA -----#
train_full <- read.csv("processed_data/ml_LCP2_proteomics.csv", check.names = FALSE)
train <- train_full %>%
  dplyr::select(-Sample.ID, -Patient.ID) %>%
  mutate(LC = factor(LC, levels = c("Primary_LC", "No_Cancer")))   # positive class first
cat("Training rows:", nrow(train), " | predictors:", ncol(train) - 1, "\n")
print(table(train$LC))

#----- RECIPE / MODEL / WORKFLOW -----#
rf_recipe <- recipe(LC ~ ., data = train) %>%
  step_zv(all_predictors())                       # trees need no scaling

rf_spec <- rand_forest(mtry = tune(), min_n = tune(), trees = tune()) %>%
  set_engine("ranger", importance = "permutation", num.threads = 1) %>%
  set_mode("classification")

rf_wf <- workflow() %>% add_recipe(rf_recipe) %>% add_model(rf_spec)

#----- CV folds: reuse the SAME folds built once in the EN section above --------
# (identical object shared by both models; nothing rebuilt, no file loaded)

#----- mtry range: scaled to the ACTUAL surviving protein count ----------------
# The 87-cohort script hardcoded c(10, 200) for its ~2923 proteins. Reusing
# that blindly here would be wrong if LOD-filtering leaves a very different
# protein count. Instead: center the range on sqrt(p) (a common default
# reference point for classification RF -- ranger's own unspecified-mtry
# default is floor(sqrt(p))), then explore a band from sqrt(p)/3 up to
# 3*sqrt(p), same spirit as the 87-cohort's "widen beyond the default so
# tuning has room to explore" but expressed as a p-dependent formula instead
# of a hardcoded number.
p <- ncol(train) - 1
mtry_lower <- max(2L, floor(sqrt(p) / 3))
mtry_upper <- min(p, ceiling(sqrt(p) * 3))
cat(sprintf("\np = %d proteins -> mtry range [%d, %d] (sqrt(p) = %.1f)\n",
            p, mtry_lower, mtry_upper, sqrt(p)))

#----- GRID + METRICS -----#
rf_grid <- grid_regular(mtry(range = c(mtry_lower, mtry_upper)),
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
  labs(title = "RF tuning (LCP2 proteomics): CV ROC-AUC across mtry and min_n")
ggsave(file.path(RES, "rf_tuning_plot.png"), width = 8, height = 5, dpi = 150)

best_rf  <- select_best(rf_res, metric = "roc_auc")
final_rf <- finalize_workflow(rf_wf, best_rf) %>% fit(data = train)
saveRDS(final_rf, file.path(RES, "final_rf_model.rds"))

best_cv <- cv_metrics %>%
  semi_join(best_rf, by = intersect(names(best_rf), names(cv_metrics))) %>%
  dplyr::select(.metric, mean, std_err, n)
write.csv(best_cv, file.path(RES, "cv_performance_best.csv"), row.names = FALSE)
cat("\nCV performance of the best RF:\n"); print(best_cv)

#----- OUT-OF-FOLD PREDICTIONS PER SAMPLE (for comparison against EN) --------
sample_ids <- train_full$Sample.ID
oof_sample <- collect_predictions(rf_res, parameters = best_rf) %>%
  dplyr::group_by(.row) %>%
  dplyr::summarise(prob_LC = mean(.pred_Primary_LC),
                   obs = dplyr::first(LC), .groups = "drop") %>%
  dplyr::arrange(.row) %>%
  dplyr::mutate(Sample.ID = sample_ids[.row], model = "RF") %>%
  dplyr::select(Sample.ID, model, obs, prob_LC)
write.csv(oof_sample, file.path(RES, "oof_predictions.csv"), row.names = FALSE)

#----- VARIABLE IMPORTANCE ACROSS ALL FOLDS -----#
extract_rf <- function(fit, split) {
  extract_fit_parsnip(fit) %>% vip::vi() %>%
    dplyr::transmute(Variable, importance = Importance)
}

prot_lookup <- read.csv("processed_data/prot_name_lookup_LCP2.csv")

rf_imp_sum <- fold_importance(finalize_workflow(rf_wf, best_rf), folds, extract_rf) %>%
  summarise_importance() %>%
  left_join(prot_lookup, by = "Variable") %>%
  mutate(label = coalesce(Assay, Variable))
write.csv(rf_imp_sum, file.path(RES, "varimp_across_folds.csv"), row.names = FALSE)

ggsave(file.path(RES, "varimp_mean.png"),
       plot_importance(rf_imp_sum, "mean_imp",
                       "RF (LCP2 proteomics): mean importance across folds", label_col = "label"),
       width = 8, height = 7, dpi = 150)
ggsave(file.path(RES, "varimp_median.png"),
       plot_importance(rf_imp_sum, "median_imp",
                       "RF (LCP2 proteomics): median importance across folds", label_col = "label"),
       width = 8, height = 7, dpi = 150)
ggsave(file.path(RES, "varimp_selfreq.png"),
       plot_importance(rf_imp_sum, "selection_freq",
                       "RF (LCP2 proteomics): selection frequency across folds", label_col = "label"),
       width = 8, height = 7, dpi = 150)

plan(sequential)   # release the parallel workers
message("\nDone. RF (LCP2 proteomics) results saved in ", RES, "/")


# ================= analyse_topprot_sex_vs_cancer.R (optional check) =================

# =============================================================================
# analyse_topprot_sex_vs_cancer.R
# -----------------------------------------------------------------------------
# For the LCP2 EN-selected proteins, disentangle a SEX effect from a CANCER
# effect on NPX -- to check whether the most-frequently-selected proteins are
# picking up lung-cancer biology or just sex biology (the HPA "female-enriched"
# worry). Sex is balanced across LC in LCP2 (46/46 F, 39/40 M), so sex cannot
# CONFOUND the case/control comparison; this script asks the sharper question of
# whether each selected protein's variance is driven by cancer or by sex.
#
# For each top protein it fits, on the ~171 LCP2 samples:
#   npx ~ Sex + LC            -> effect of each, ADJUSTED for the other (+ p)
#   npx ~ Sex * LC            -> interaction (does the cancer effect differ by sex)
#   eta^2 (variance explained) for Sex vs LC, from the additive model
#   effect sizes reported in NPX units AND as fold-change (1 NPX = 2x)
# then BH-corrects across the tested proteins and classifies each as
# "cancer-driven" or "sex-driven" by which explains more variance.
#
# Reads EN results that already exist after en_LCP2_proteomics.R has run.
# Run with working dir = LC_Data, AFTER en_LCP2_proteomics.R.
# =============================================================================

## ---- Packages ---------------------------------------------------------------
required <- c("tidyverse", "broom", "ggrepel")
to_install <- required[!required %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)
library(tidyverse); library(broom); library(ggrepel)

RES_EN <- "results/en_LCP2_proteomics"
OUT    <- file.path(RES_EN, "sex_vs_cancer")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

TOP_N <- 20   # how many top proteins to examine (bump if you want more)

## ---- 1) Pick the top proteins by CV selection frequency ----------------------
# "Most frequently selected" = selection_freq across the 100 fold-fits, which is
# exactly the stability notion you described. Falls back to mean |coef|.
vi <- read.csv(file.path(RES_EN, "varimp_across_folds.csv"))
top <- vi %>%
  arrange(desc(selection_freq), desc(mean_imp)) %>%
  slice_head(n = TOP_N)
top_vars <- top$Variable
labels   <- setNames(dplyr::coalesce(top$label, top$Variable), top$Variable)
cat("Examining top", length(top_vars), "proteins by selection frequency.\n")

## ---- 2) Modelling table + Sex from metadata ----------------------------------
ml   <- read.csv("processed_data/ml_LCP2_proteomics.csv", check.names = FALSE)
meta <- read.csv("metadata/LCP2_metadata_Olink.csv") %>% dplyr::select(Sample.ID, Sex)

df <- ml %>%
  left_join(meta, by = "Sample.ID") %>%
  mutate(LC  = factor(LC,  levels = c("No_Cancer", "Primary_LC")),  # ref = control
         Sex = factor(Sex, levels = c("Male", "Female")))           # ref = Male
# (Sex ref = Male so a positive Sex coefficient means "higher in females".)

stopifnot(all(top_vars %in% names(df)))   # every selected protein is a column

## ---- 3) Per-protein sex-vs-cancer decomposition ------------------------------
analyse_one <- function(v) {
  d <- df %>% transmute(npx = .data[[v]], Sex, LC) %>% tidyr::drop_na()

  # marginal (unadjusted) group means -> effect sizes in NPX (= log2 fold change)
  sex_delta <- mean(d$npx[d$Sex == "Female"])     - mean(d$npx[d$Sex == "Male"])
  lc_delta  <- mean(d$npx[d$LC  == "Primary_LC"]) - mean(d$npx[d$LC  == "No_Cancer"])

  # additive model: each effect adjusted for the other
  fit <- lm(npx ~ Sex + LC, data = d)
  co  <- broom::tidy(fit)
  p_sex <- co$p.value[co$term == "SexFemale"]
  p_lc  <- co$p.value[co$term == "LCPrimary_LC"]

  # variance explained (eta^2); design is balanced so SS type doesn't matter
  a   <- anova(fit)
  ss  <- setNames(a$`Sum Sq`, rownames(a))
  tot <- sum(ss)

  # interaction: does the cancer effect differ by sex?
  p_int <- anova(lm(npx ~ Sex * LC, data = d))["Sex:LC", "Pr(>F)"]

  tibble(
    Variable        = v,
    Assay           = labels[[v]],
    cancer_delta_NPX = lc_delta,  cancer_FC = 2^lc_delta,  p_cancer = p_lc,  eta2_cancer = ss["LC"]  / tot,
    sex_delta_NPX    = sex_delta, sex_FC    = 2^sex_delta, p_sex    = p_sex, eta2_sex    = ss["Sex"] / tot,
    p_interaction    = p_int
  )
}

res <- purrr::map_dfr(top_vars, analyse_one) %>%
  mutate(
    q_cancer      = p.adjust(p_cancer, "BH"),
    q_sex         = p.adjust(p_sex,    "BH"),
    q_interaction = p.adjust(p_interaction, "BH"),
    driver        = ifelse(eta2_cancer >= eta2_sex, "cancer-driven", "sex-driven")
  ) %>%
  arrange(desc(eta2_cancer))

write.csv(res, file.path(OUT, "top_protein_sex_vs_cancer.csv"), row.names = FALSE)
cat("\nSex-vs-cancer decomposition of the top selected proteins:\n")
print(res %>% dplyr::select(Assay, driver, eta2_cancer, eta2_sex,
                            cancer_FC, sex_FC, q_cancer, q_sex, q_interaction), n = TOP_N)
cat(sprintf("\n%d/%d top proteins are cancer-driven (cancer explains more variance than sex).\n",
            sum(res$driver == "cancer-driven"), nrow(res)))

## ---- 4) Money plot: variance explained by cancer vs by sex -------------------
# Above the dashed 1:1 line = cancer explains more of that protein's variance
# than sex does (what you want your selected proteins to be).
p_scatter <- ggplot(res, aes(eta2_sex, eta2_cancer, label = Assay)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = driver), size = 2.4) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = Inf) +
  scale_colour_manual(values = c("cancer-driven" = "firebrick", "sex-driven" = "steelblue"),
                      name = NULL) +
  labs(title = "LCP2 top EN proteins: variance explained by cancer vs sex",
       subtitle = "Above the dashed line = cancer explains more than sex (real disease signal)",
       x = expression("Variance explained by Sex ("*eta^2*")"),
       y = expression("Variance explained by Cancer ("*eta^2*")")) +
  coord_equal() + theme_minimal()
ggsave(file.path(OUT, "sex_vs_cancer_variance.png"), p_scatter, width = 8, height = 7, dpi = 150)

## ---- 5) Boxplots: NPX by cancer status WITHIN each sex, top 9 proteins --------
long <- df %>%
  dplyr::select(Sex, LC, all_of(top_vars)) %>%
  pivot_longer(all_of(top_vars), names_to = "Variable", values_to = "npx") %>%
  mutate(Assay = labels[Variable])

p_box <- ggplot(long %>% filter(Variable %in% head(top_vars, 9)),
                aes(interaction(LC, Sex), npx, fill = LC)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ Assay, scales = "free_y") +
  scale_fill_manual(values = c("No_Cancer" = "grey70", "Primary_LC" = "firebrick"), name = NULL) +
  labs(title = "Top selected proteins: NPX by cancer status within each sex",
       x = NULL, y = "NPX") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT, "top_protein_boxplots.png"), p_box, width = 11, height = 8, dpi = 150)

message("\nDone. Sex-vs-cancer outputs in ", OUT, "/  ",
        "(top_protein_sex_vs_cancer.csv, sex_vs_cancer_variance.png, top_protein_boxplots.png)")
