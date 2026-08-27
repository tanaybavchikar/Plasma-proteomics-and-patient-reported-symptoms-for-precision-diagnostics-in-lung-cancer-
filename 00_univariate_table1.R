# ============================================================
# 00_univariate_table1.R  — Task 1: univariate description ('Table 1')
# from task1_univariate/table1_univariate.R
# Run with working directory = LC_Data.
# NOTE: faithful merge of the scripts named above; only data paths repointed
# to raw_data/ and metadata/. Preserved originals are in code/archive/.
# ============================================================

# =============================================================================
# table1_univariate.R
# -----------------------------------------------------------------------------
# Task 1 (Noora): the descriptive "Table 1" + univariate case-vs-control tests
# for the 411-patient TRAINING cohort. NO machine learning here -- this is the
# per-variable description every clinical/proteomics paper leads with, and the
# sanity baseline for which variables individually separate cancer from control.
#
# For every variable it reports frequency overall / in cases / in controls, and
# tests case-vs-control:
#   - binary variable (symptom or background yes/no): 2x2 table, ODDS RATIO
#     (event = cancer, so OR > 1 = more common in cancer) with a Woolf 95% CI
#     (Haldane-Anscombe +0.5 if any zero cell), tested with FISHER's exact test.
#   - continuous (Age): WILCOXON rank-sum (Mann-Whitney), reported as median [IQR]
#     per group; no OR.
# Then BH / FDR correction across ALL tests (one family) -> q-values, because
# ~130 univariate tests will throw false positives at raw p<0.05.
#
# Variables are tagged BACKGROUND (Age, Gender, the Q* history/demographic block)
# vs SYMPTOM (the module series Br_/Co_/Ph_/Pa_/Fa_/Vo_/App_/Sm_/Fe_/Oth_), so
# the outputs can be read/plotted per group.
#
# Data: the 411 training patients are taken from symptoms_clean.csv (which has
# the background block) by keeping only Patient IDs present in symptoms_train.csv
# (symptoms_train.csv itself is symptom-only, so it can't supply background).
#
# Outputs (results/table1_univariate/):
#   table1_univariate.csv        full table (freqs, OR/CI, test, p, q, meaning)
#   table1_background.png        background-factor frequencies, cases vs controls
#   table1_symptoms_sig.png      most significant symptoms, cases vs controls
#   table1_age_boxplot.png       Age by group (Wilcoxon)
#   table1_volcano.png           log2(OR) vs -log10(q), all binary variables
#
# Run with working dir = LC_Data.
# =============================================================================

library(tidyverse)
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
have_readxl <- requireNamespace("readxl", quietly = TRUE)

PROC <- "processed_data"
OUT  <- "results/table1_univariate"
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

TOP_SYMPTOMS_PLOT <- 30   # how many symptoms to show in the symptom bar plot (by p)

sanitize <- function(x) make.names(gsub("[^A-Za-z0-9_]+", "_", x))

## ---- 1) Data: full clean table restricted to the 411 TRAINING patients -------
clean     <- read.csv(file.path(PROC, "symptoms_clean.csv"), check.names = FALSE)
train_ids <- read.csv(file.path(PROC, "symptoms_train.csv"), check.names = FALSE)$Patient

dat <- clean %>% dplyr::filter(Patient %in% train_ids)
stopifnot(nrow(dat) == length(train_ids))          # every training patient found
names(dat) <- sanitize(names(dat))
cat(sprintf("Table 1 cohort: %d training patients\n", nrow(dat)))

## ---- 2) Outcome + variable roles ---------------------------------------------
# event = cancer (so OR > 1 means "more common in cancer"). control is the ref.
dat$outcome <- factor(ifelse(dat$Lung_cancer == 1, "cancer", "control"),
                      levels = c("control", "cancer"))
cat("Label balance: ",
    paste(names(table(dat$outcome)), table(dat$outcome), sep = "=", collapse = ", "), "\n")

exclude          <- c("Patient", "Lung_cancer", "Stage", "outcome")
symptom_prefixes <- c("Br_", "Co_", "Ph_", "Pa_", "Fa_", "Vo_", "App_", "Sm_", "Fe_", "Oth_")
vars <- setdiff(names(dat), exclude)

is_symptom <- function(v) any(startsWith(v, symptom_prefixes))
var_group  <- ifelse(vapply(vars, is_symptom, logical(1)), "symptom", "background")
names(var_group) <- vars

n_unique        <- vapply(vars, function(v) length(unique(na.omit(dat[[v]]))), integer(1))
continuous_vars <- vars[n_unique > 10]      # Age
binary_vars     <- vars[n_unique == 2]
const_vars      <- vars[n_unique < 2]
if (length(const_vars))
  cat("Dropping", length(const_vars), "constant variable(s):",
      paste(const_vars, collapse = ", "), "\n")

# Make sure every "binary" variable really is coded 0/1 (1 = present); warn+recode
for (v in binary_vars) {
  u <- sort(unique(na.omit(dat[[v]])))
  if (!identical(as.numeric(u), c(0, 1))) {
    dat[[v]] <- as.integer(dat[[v]] == max(u))
    warning(sprintf("Recoded %s from {%s} to 0/1 (1 = %s)", v, paste(u, collapse = ","), max(u)))
  }
}

## ---- 3) Human-readable meanings (optional; falls back to the variable name) ---
meaning_map <- tibble(Variable = vars, meaning = vars)
if (have_readxl) {
  dict <- tryCatch(
    readxl::read_excel("metadata/113 Merged file_cleanced.xlsx",
                       sheet = "Variable_explaination", col_names = c("var", "meaning")),
    error = function(e) NULL)
  if (!is.null(dict)) {
    dict <- dict %>% mutate(Variable = sanitize(var)) %>% distinct(Variable, .keep_all = TRUE)
    meaning_map <- meaning_map %>%
      left_join(dplyr::select(dict, Variable, meaning2 = meaning), by = "Variable") %>%
      mutate(meaning = dplyr::coalesce(meaning2, meaning)) %>%
      dplyr::select(Variable, meaning)
  }
}

## ---- 4) Per-variable tests ----------------------------------------------------
# binary: 2x2 -> OR (Haldane if needed) + Woolf CI + Fisher's exact test
binary_row <- function(v) {
  x <- dat[[v]]; y <- dat$outcome
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]

  a <- sum(x == 1 & y == "cancer");  cc <- sum(x == 0 & y == "cancer")   # cancer: pos, neg
  b <- sum(x == 1 & y == "control"); d  <- sum(x == 0 & y == "control")  # control: pos, neg
  freq_case <- a / (a + cc); freq_ctrl <- b / (b + d)

  aa <- a; bb <- b; ccc <- cc; dd <- d
  if (any(c(a, b, cc, d) == 0)) { aa <- a + .5; bb <- b + .5; ccc <- cc + .5; dd <- d + .5 }
  or  <- (aa * dd) / (bb * ccc)
  se  <- sqrt(1/aa + 1/bb + 1/ccc + 1/dd)
  ci_l <- exp(log(or) - 1.96 * se); ci_u <- exp(log(or) + 1.96 * se)
  p   <- fisher.test(matrix(c(a, cc, b, d), nrow = 2))$p.value

  tibble(Variable = v, group = var_group[[v]], type = "binary",
         n_case = a + cc, n_ctrl = b + d, pos_case = a, pos_ctrl = b,
         freq_case = freq_case, freq_ctrl = freq_ctrl,
         effect_type = "OR", effect = or, CI_low = ci_l, CI_high = ci_u,
         test = "Fisher exact", p = p,
         summary_case = sprintf("%d/%d (%.0f%%)", a, a + cc, 100 * freq_case),
         summary_ctrl = sprintf("%d/%d (%.0f%%)", b, b + d, 100 * freq_ctrl))
}

# continuous: median [IQR] per group + Wilcoxon rank-sum
continuous_row <- function(v) {
  x <- dat[[v]]; y <- dat$outcome
  wt <- wilcox.test(x ~ y)
  mc <- median(x[y == "cancer"],  na.rm = TRUE); ic <- IQR(x[y == "cancer"],  na.rm = TRUE)
  mk <- median(x[y == "control"], na.rm = TRUE); ik <- IQR(x[y == "control"], na.rm = TRUE)
  tibble(Variable = v, group = var_group[[v]], type = "continuous",
         n_case = sum(y == "cancer" & !is.na(x)), n_ctrl = sum(y == "control" & !is.na(x)),
         pos_case = NA_real_, pos_ctrl = NA_real_,
         freq_case = NA_real_, freq_ctrl = NA_real_,
         effect_type = "median_diff", effect = mc - mk, CI_low = NA_real_, CI_high = NA_real_,
         test = "Wilcoxon", p = wt$p.value,
         summary_case = sprintf("%.0f [%.0f]", mc, ic),
         summary_ctrl = sprintf("%.0f [%.0f]", mk, ik))
}

tab <- bind_rows(
  purrr::map_dfr(binary_vars,     binary_row),
  purrr::map_dfr(continuous_vars, continuous_row)
) %>%
  left_join(meaning_map, by = "Variable") %>%
  mutate(label = dplyr::coalesce(meaning, Variable),
         q = p.adjust(p, method = "BH")) %>%           # FDR across ALL tests (one family)
  relocate(label, meaning, .after = Variable) %>%
  arrange(group, p)

write.csv(tab, file.path(OUT, "table1_univariate.csv"), row.names = FALSE)
cat(sprintf("\nTable 1 written: %d variables (%d background, %d symptom). %d significant at q<0.05.\n",
            nrow(tab), sum(tab$group == "background"), sum(tab$group == "symptom"),
            sum(tab$q < 0.05, na.rm = TRUE)))
cat("\nMost significant variables overall:\n")
print(tab %>% arrange(p) %>%
        dplyr::select(label, group, summary_case, summary_ctrl, effect_type, effect, p, q) %>%
        head(15))

## ---- 5) Plots -----------------------------------------------------------------
# helper: dodged case-vs-control frequency bars for a set of binary rows
freq_barplot <- function(df, title) {
  df %>%
    dplyr::select(label, effect, freq_case, freq_ctrl) %>%
    tidyr::pivot_longer(c(freq_case, freq_ctrl), names_to = "grp", values_to = "freq") %>%
    mutate(grp = recode(grp, freq_case = "cancer", freq_ctrl = "control"),
           label = forcats::fct_reorder(label, effect)) %>%
    ggplot(aes(label, freq, fill = grp)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    coord_flip() +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c("cancer" = "firebrick", "control" = "grey60"), name = NULL) +
    labs(title = title, x = NULL, y = "Frequency (% of group)") +
    theme_minimal()
}

# (a) background factors (binary ones)
bg <- tab %>% dplyr::filter(group == "background", type == "binary")
ggsave(file.path(OUT, "table1_background.png"),
       freq_barplot(bg, "Background factors: frequency in cancer vs control"),
       width = 9, height = max(4, 0.32 * nrow(bg)), dpi = 150)

# (b) most significant symptoms
sy <- tab %>% dplyr::filter(group == "symptom", type == "binary") %>%
  arrange(p) %>% head(TOP_SYMPTOMS_PLOT)
ggsave(file.path(OUT, "table1_symptoms_sig.png"),
       freq_barplot(sy, sprintf("Top %d symptoms by significance: cancer vs control", nrow(sy))),
       width = 9, height = max(5, 0.32 * nrow(sy)), dpi = 150)

# (c) Age boxplot
if ("Age" %in% continuous_vars) {
  p_age <- wilcox.test(Age ~ outcome, data = dat)$p.value
  ggplot(dat, aes(outcome, Age, fill = outcome)) +
    geom_boxplot(width = 0.5, outlier.size = 0.6) +
    scale_fill_manual(values = c("control" = "grey60", "cancer" = "firebrick"), guide = "none") +
    labs(title = sprintf("Age by group (Wilcoxon p = %.3g)", p_age), x = NULL, y = "Age") +
    theme_minimal()
  ggsave(file.path(OUT, "table1_age_boxplot.png"), width = 5, height = 5, dpi = 150)
}

# (d) volcano: log2(OR) vs -log10(q) for all binary variables
volc <- tab %>% dplyr::filter(type == "binary") %>%
  mutate(log2OR = log2(effect), neglog10q = -log10(q),
         sig = ifelse(q < 0.05, ifelse(effect > 1, "up in cancer (q<0.05)", "down in cancer (q<0.05)"), "ns"))
ggplot(volc, aes(log2OR, neglog10q, colour = sig)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey80") +
  geom_point(size = 2, alpha = 0.8) +
  ggrepel::geom_text_repel(data = dplyr::filter(volc, q < 0.05),
                           aes(label = label), size = 2.7, max.overlaps = 20) +
  scale_colour_manual(values = c("up in cancer (q<0.05)" = "firebrick",
                                 "down in cancer (q<0.05)" = "steelblue", "ns" = "grey70"),
                      name = NULL) +
  labs(title = "Univariate case-vs-control associations (all binary variables)",
       x = "log2(odds ratio)  -- right = more common in cancer", y = "-log10(BH q-value)") +
  theme_minimal()
ggsave(file.path(OUT, "table1_volcano.png"), width = 9, height = 7, dpi = 150)

message("\nDone. Table 1 + plots in ", OUT, "/")


# =============================================================================
# MULTIVARIATE ANALYSIS  (commented out -- uncomment the block below to run it)
# -----------------------------------------------------------------------------
# Univariate above asks "is this variable associated with cancer on its own?".
# Multivariable logistic regression asks "is it associated INDEPENDENTLY of the
# others?" -- e.g. does 'living alone' still predict cancer once Age is in the
# model, or was it just a proxy for being older? Two complementary pieces:
#   (A) age + sex ADJUSTED OR for every binary variable (each fit in its own
#       little model alongside Age + Gender) -- the cleanest confounding check,
#       and a crude-vs-adjusted comparison so you can see which associations
#       shrink once age/sex are accounted for.
#   (B) ONE multivariable model on the univariate hits (p < 0.05) + Age/Gender,
#       giving mutually-adjusted ORs (a forest plot). Mind events-per-variable:
#       ~150 cancer events / (#predictors); keep #predictors under ~15 or the
#       ORs get unstable/overfit. That's why we screen to univariate hits first
#       rather than throwing all ~130 variables in at once.
# Note on sign: outcome levels are c("control","cancer"), so base glm() models
# P(cancer) and exp(coef) is ALREADY a cancer-oriented OR -- no sign flip needed
# here (this is base glm, not glmnet).
#
# library(broom)
#
# # ---- (A) Age + sex adjusted OR per binary variable ------------------------
# adj_or_row <- function(v) {
#   fit <- tryCatch(
#     glm(reformulate(c(v, "Age", "Gender"), response = "outcome"),
#         data = dat, family = binomial),
#     error = function(e) NULL)
#   if (is.null(fit)) return(NULL)
#   s  <- summary(fit)$coefficients
#   rn <- rownames(s)[2]                       # row 2 = the variable (after intercept)
#   est <- s[rn, "Estimate"]; se <- s[rn, "Std. Error"]; p <- s[rn, "Pr(>|z|)"]
#   tibble(Variable = v, adj_OR = exp(est),
#          adj_CI_low = exp(est - 1.96 * se), adj_CI_high = exp(est + 1.96 * se),
#          adj_p = p)
# }
# adj_vars <- setdiff(binary_vars, "Gender")   # Gender is a covariate, not a test target
# adj_tab  <- purrr::map_dfr(adj_vars, adj_or_row) %>%
#   mutate(adj_q = p.adjust(adj_p, "BH")) %>%
#   left_join(meaning_map, by = "Variable") %>%
#   arrange(adj_p)
# write.csv(adj_tab, file.path(OUT, "table1_multivariate_agesex_adjusted.csv"), row.names = FALSE)
#
# # crude (univariate) OR vs age/sex-adjusted OR -- does the association survive?
# compare <- tab %>% dplyr::select(Variable, label, crude_OR = effect, crude_p = p) %>%
#   inner_join(dplyr::select(adj_tab, Variable, adj_OR, adj_p), by = "Variable") %>%
#   arrange(crude_p)
# write.csv(compare, file.path(OUT, "table1_crude_vs_adjusted.csv"), row.names = FALSE)
# cat("\nCrude vs age/sex-adjusted OR (top 15 by crude p):\n"); print(head(compare, 15))
#
# # ---- (B) One multivariable model on the univariate hits + Age/Gender ------
# cand  <- tab %>% dplyr::filter(type == "binary", p < 0.05) %>% pull(Variable)
# preds <- unique(c("Age", "Gender", setdiff(cand, "Gender")))
# cat(sprintf("\nMultivariable model: %d predictors, %d cancer events (%.1f events/predictor)\n",
#             length(preds), sum(dat$outcome == "cancer"),
#             sum(dat$outcome == "cancer") / length(preds)))
# mv_fit <- glm(reformulate(preds, response = "outcome"), data = dat, family = binomial)
# mv_tab <- broom::tidy(mv_fit, conf.int = TRUE, exponentiate = TRUE) %>%
#   dplyr::filter(term != "(Intercept)") %>%
#   mutate(Variable = term) %>%
#   left_join(meaning_map, by = "Variable") %>%
#   mutate(label = dplyr::coalesce(meaning, Variable))
# write.csv(mv_tab, file.path(OUT, "table1_multivariable_model.csv"), row.names = FALSE)
#
# # forest plot of mutually-adjusted ORs (log scale; OR>1 = higher cancer odds)
# ggplot(mv_tab, aes(estimate, reorder(label, estimate))) +
#   geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
#   geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
#   scale_x_log10() +
#   labs(title = "Multivariable logistic regression: adjusted odds ratios",
#        subtitle = "OR > 1 = higher cancer odds, adjusted for every other term",
#        x = "Adjusted odds ratio (log scale)", y = NULL) +
#   theme_minimal()
# ggsave(file.path(OUT, "table1_multivariable_forest.png"), width = 8, height = 6, dpi = 150)
# =============================================================================
