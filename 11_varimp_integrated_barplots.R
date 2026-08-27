# ============================================================
# 11_varimp_integrated_barplots.R  — importance bars for the INTEGRATED models
# ------------------------------------------------------------
# ONE comparative figure: the variable importance from the integrated
# (proteins + symptoms+bg) models trained on ALL 87 patients, faceted by run
# (EN/RF x impfrac/selfreq) so you can compare them side by side. Each bar is
# coloured by modality (protein vs symptom/bg) to show which modality the
# combined model leans on.
#
# Importance is SCALED within each run (imp / max imp) so EN (|coef|) and RF
# (permutation importance) share a common 0..1 axis and the panels are
# comparable. TOP_N keeps each panel readable.
#
# Run with working directory = LC_Data, after models_87.R.
# ============================================================

library(tidyverse)
if (!requireNamespace("tidytext", quietly = TRUE)) install.packages("tidytext")
library(tidytext)   # reorder_within / scale_y_reordered (per-facet ordering)

DIR   <- "results/08_early_integrated_87"
TOP_N <- 15         # top features to show PER panel
mod_cols <- c("protein" = "#E76B8A", "symptom/bg" = "#2E8C8A")

classify <- function(v) ifelse(startsWith(v, "prot_"), "protein", "symptom/bg")
pretty   <- function(v) sub("^prot_", "", v)

runs <- tidyr::expand_grid(model = c("EN", "RF"), crit = c("_impfrac", "_selfreq")) %>%
  dplyr::mutate(tag = paste0("integrated_", model, crit),
                run = paste0(model, " · ", sub("_", "", crit)))

dat <- purrr::pmap_dfr(runs, function(model, crit, tag, run) {
  f <- file.path(DIR, paste0("varimp_final_", tag, ".csv"))
  if (!file.exists(f)) { message("SKIP (missing): ", f); return(NULL) }
  read.csv(f) %>%
    dplyr::filter(importance > 0) %>%
    dplyr::mutate(scaled = importance / max(importance, na.rm = TRUE),
                  modality = classify(Variable), label = pretty(Variable), run = run) %>%
    dplyr::arrange(dplyr::desc(scaled)) %>%
    dplyr::slice_head(n = TOP_N)
})
if (is.null(dat) || !nrow(dat))
  stop("No varimp_final_integrated_*.csv in ", DIR, ". Run models_87.R first.")

dat$run <- factor(dat$run, levels = c("EN · impfrac", "EN · selfreq",
                                      "RF · impfrac", "RF · selfreq"))

p <- ggplot(dat, aes(x = scaled,
                     y = reorder_within(label, scaled, run),
                     fill = modality)) +
  geom_col(width = 0.74) +
  scale_y_reordered() +
  scale_fill_manual(values = mod_cols, name = NULL) +
  facet_wrap(~ run, scales = "free_y", ncol = 2) +
  labs(title = "Integrated model — variable importance (all-87 fit)",
       subtitle = "scaled within each run (share of the top feature); proteins vs symptoms/bg",
       x = "scaled importance (top feature = 1)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))

out <- file.path(DIR, "varimp_final_integrated_compare.png")
ggsave(out, p, width = 11, height = 9, dpi = 150)
cat("wrote", out, "\n")
message("\nDone. Comparative integrated-importance figure in ", out)
