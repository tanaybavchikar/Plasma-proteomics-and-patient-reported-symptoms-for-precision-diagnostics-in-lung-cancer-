# ============================================================
# 02_create_ML_data.R  — build all ML tables + CV folds
# merged: build_ml_tables + build_ml_tables87 + build_ml_tables_symptbg + build_ml_tables_LCP2 + make_folds + make_folds87 + make_folds_LCP2
# Run with working directory = LC_Data.
# NOTE: faithful merge of the scripts named above; only data paths repointed
# to raw_data/ and metadata/. Preserved originals are in code/archive/.
# ============================================================

# =============================================================================
# build_ml_tables.R
# -----------------------------------------------------------------------------
# Builds the train/test tables for the symptoms model (LC vs control).
#
# Design:
#   TRAIN = the larger cohort (411) -- patients WITHOUT proteomics.
#   TEST  = the 87 patients who also have proteomics (Olink, minus LCP44).
#
# Outcome: the 505 file's OWN `Lung_cancer` column (1 = cancer, 0 = control).
#   It is complete (498/498) and agrees perfectly with the 442 EORTC label, so we
#   use it directly -- no metadata fallback, no dropped patients.
#
# *** LEAKAGE FIX: `Lung_cancer` and `Stage` are the diagnosis itself. They must
#     NOT be predictors -- we extract the label, then DROP both columns. (Leaving
#     them in gives a fake ROC-AUC of 1.0.) ***
#
# Run with working dir = LC_Data. Outputs to processed_data/.
# =============================================================================

library(tidyverse)
library(readxl)
PROC <- "processed_data"

# ---- Cleaned symptom cohort (498) + features --------------------------------
sym <- read.csv(file.path(PROC, "symptoms_clean.csv"), check.names = FALSE)

# ---- Label from the 505's own Lung_cancer column ----------------------------
sym <- sym %>%
  mutate(Lungcancer = factor(ifelse(Lung_cancer == 1, "Yes", "No"),
                             levels = c("No", "Yes")))

# ---- Remove leakage / non-feature diagnosis columns -------------------------
sym <- sym %>% dplyr::select(-Lung_cancer, -Stage)

# ---- Keep ONLY the symptom features listed in the 113 file's tab 2 ----------
# Drops Age, Gender, background Q's, and the 16 over-granular Pain columns that
# aren't in the variable dictionary -> a pure curated-symptoms model (110 feats).
dict113 <- read_excel("metadata/113 Merged file_cleanced.xlsx",
                      sheet = "Variable_explaination", col_names = c("var", "desc"))
keep_feats <- intersect(names(sym), dict113$var)
sym <- sym %>% dplyr::select(Patient, Lungcancer, all_of(keep_feats))

# Make names ASCII + syntactically valid here (once), so every downstream script
# and the shared CV folds use identical clean names (Swedish a/a/o break models).
names(sym) <- make.names(gsub("[^A-Za-z0-9_]+", "_", names(sym)), unique = TRUE)
message("Symptom features kept (in tab 2): ", length(keep_feats))

# ---- Subject IDs of the proteomics (Olink) patients -------------------------
olink_lcp <- read_delim("raw_data/VB-3207_NPX_2022-09-15.csv", delim = ";",
                        show_col_types = FALSE) %>%
  filter(str_starts(SampleID, "LCP"), SampleID != "LCP44") %>%
  distinct(SampleID) %>% pull(SampleID)

prot_subjects <- read_excel("metadata/113 Merged file_cleanced.xlsx", sheet = "Data") %>%
  filter(!is.na(Patient)) %>%
  filter(LCP_ID %in% olink_lcp) %>%
  pull(Patient)

# ---- Split: large cohort (train) vs the 87 proteomics patients (test) -------
dat <- sym %>% mutate(in_proteomics = Patient %in% prot_subjects)

train <- dat %>% filter(!in_proteomics) %>% dplyr::select(-in_proteomics)  # large cohort
test  <- dat %>% filter(in_proteomics)  %>% dplyr::select(-in_proteomics)  # 87 proteomics

message("TRAIN (large cohort): ", nrow(train),
        " | TEST (87 proteomics): ", nrow(test))
message("TRAIN balance: ", paste(names(table(train$Lungcancer)),
        table(train$Lungcancer), sep = "=", collapse = ", "))
message("TEST balance:  ", paste(names(table(test$Lungcancer)),
        table(test$Lungcancer), sep = "=", collapse = ", "))

write.csv(train, file.path(PROC, "symptoms_train.csv"),  row.names = FALSE)
write.csv(test,  file.path(PROC, "symptoms_test87.csv"), row.names = FALSE)

# ---- Sanity checks ----------------------------------------------------------
stopifnot(
  nrow(train) == 411,
  nrow(test)  == 87,
  length(keep_feats) == 110,                           # tab-2 symptom features only
  !any(c("Lung_cancer", "Stage", "Age", "Gender") %in% names(train)),  # no leak / non-features
  length(intersect(train$Patient, test$Patient)) == 0,
  all(c("No", "Yes") %in% as.character(train$Lungcancer)),
  all(c("No", "Yes") %in% as.character(test$Lungcancer))
)
message("All table-building checks passed (no leakage columns).")


# ================= build_ml_tables87.R =================

# =============================================================================
# build_ml_tables87.R
# -----------------------------------------------------------------------------
# Builds THREE modelling tables on the SAME 87 proteomics patients, so that
# symptoms-only, proteomics-only, and combined models are compared fairly:
#
#   ml87_symptoms.csv    Patient, Lungcancer, 110 symptom features
#   ml87_proteomics.csv  Patient, Lungcancer, ~2923 Olink proteins (prot_*)
#   ml87_combined.csv    Patient, Lungcancer, 110 symptoms + ~2923 proteins
#
# All three are written in the SAME row order (sorted by Patient) so a single
# row-index CV-folds object (make_folds87.R -> cv_folds87.rds) applies to all.
#
# Data notes:
#   * Olink NPX is already log2-normalised and QC-cleaned in
#     processed_data/LCP1_olink_clean_wide.txt (proteins in rows, LCP in cols)
#     -> we transpose to patients-in-rows.
#   * Two ID systems: Olink uses LCP IDs, the symptom table uses numeric Patient
#     IDs. We bridge them with the 113 file's "Data" sheet (LCP_ID <-> Patient),
#     exactly as build_ml_tables.R does.
#   * Protein names are sanitised and prefixed "prot_" so (a) odd names like
#     HLA-A / ERVV-1 don't break model code, and (b) in the combined table you
#     can always tell a protein from a symptom. A lookup (sanitised <-> Assay)
#     is written to processed_data/prot_name_lookup.csv.
#   * Below-LOD handling: Olink reports values below the limit of detection; a
#     protein detected in almost no one is mostly noise. PROT_MAX_MISSING drops
#     proteins whose fraction-below-LOD (MissingFreq) exceeds the cutoff. Default
#     = 1 (keep all; let elastic net shrink them). The console prints how many
#     you'd drop at 0.5 / 0.25 so you can tighten it later.
#
# Run with working dir = LC_Data, after build_ml_tables.R (needs symptoms_test87).
# =============================================================================

library(tidyverse)
library(readxl)
PROC <- "processed_data"

## ---- tunable: below-LOD (MissingFreq) filter --------------------------------
PROT_MAX_MISSING <- 1.0   # keep all proteins. Set e.g. 0.5 to drop noisy assays.

## ---- 1) Symptoms-only table for the 87 --------------------------------------
# symptoms_test87.csv already has Patient, Lungcancer, and the 110 sanitised
# tab-2 symptom features on exactly these 87 patients (from build_ml_tables.R).
sym87 <- read.csv(file.path(PROC, "symptoms_test87.csv"), check.names = FALSE) %>%
  arrange(Patient)
stopifnot(nrow(sym87) == 87, all(c("Patient", "Lungcancer") %in% names(sym87)))
label87 <- sym87 %>% dplyr::select(Patient, Lungcancer)   # the shared label

## ---- 2) LCP <-> Patient crosswalk (113 "Data" sheet) ------------------------
# NOTE: readxl reads numeric cells as DOUBLE, never integer, whereas read.csv()
# gives INTEGER for whole numbers. Left as-is, prot87$Patient would be double and
# sym87$Patient integer, so the identical() check at the bottom fails on TYPE even
# though every value and the ordering match. Coerce here, at the source.
xwalk <- read_excel("metadata/113 Merged file_cleanced.xlsx", sheet = "Data") %>%
  filter(!is.na(Patient), !is.na(LCP_ID)) %>%
  mutate(Patient = as.integer(Patient)) %>%
  distinct(LCP_ID, Patient)

## ---- 3) Olink wide -> patients-in-rows --------------------------------------
olink_wide <- read_delim(file.path(PROC, "LCP1_olink_clean_wide.txt"),
                         delim = "\t", show_col_types = FALSE)
# first column is the protein name ("Assay"); the rest are LCP sample columns
assay_col <- names(olink_wide)[1]
proteins  <- olink_wide[[assay_col]]

# optional below-LOD filter using MissingFreq from the long file ---------------
miss <- read_delim(file.path(PROC, "LCP1_olink_clean_long.txt"),
                   delim = "\t", show_col_types = FALSE) %>%
  distinct(Assay, MissingFreq)
mfreq <- tibble(Assay = proteins) %>%
  left_join(miss, by = "Assay") %>%
  mutate(MissingFreq = ifelse(is.na(MissingFreq), 0, MissingFreq))
cat(sprintf("Proteins total: %d | would drop at MissingFreq>0.5: %d | >0.25: %d\n",
            nrow(mfreq), sum(mfreq$MissingFreq > 0.5), sum(mfreq$MissingFreq > 0.25)))
keep_prot <- mfreq$MissingFreq <= PROT_MAX_MISSING
cat(sprintf("Keeping %d proteins (PROT_MAX_MISSING = %.2f).\n",
            sum(keep_prot), PROT_MAX_MISSING))

olink_wide <- olink_wide[keep_prot, ]
proteins   <- proteins[keep_prot]

# sanitise + prefix protein names; keep a lookup back to the original Assay
prot_clean <- paste0("prot_", make.names(gsub("[^A-Za-z0-9_]+", "_", proteins),
                                         unique = TRUE))
prot_lookup <- tibble(Variable = prot_clean, Assay = proteins)
write.csv(prot_lookup, file.path(PROC, "prot_name_lookup.csv"), row.names = FALSE)

# transpose: rows = LCP samples, cols = proteins
mat <- t(as.matrix(olink_wide[, -1]))            # LCP in rows now
colnames(mat) <- prot_clean
prot_lcp <- tibble(LCP_ID = rownames(mat)) %>% bind_cols(as_tibble(mat))

## ---- 4) Attach Patient IDs + label ------------------------------------------
prot87 <- prot_lcp %>%
  inner_join(xwalk, by = "LCP_ID") %>%
  inner_join(label87, by = "Patient") %>%
  relocate(Patient, Lungcancer) %>%
  dplyr::select(-LCP_ID) %>%
  arrange(Patient)

# every Olink sample should map to a labelled patient, and cover all 87
if (nrow(prot_lcp) != nrow(prot87))
  warning(sprintf("%d Olink samples did not map to a labelled Patient",
                  nrow(prot_lcp) - nrow(prot87)))

## ---- 5) Combined table (symptoms + proteins) --------------------------------
combined87 <- sym87 %>%
  inner_join(dplyr::select(prot87, -Lungcancer), by = "Patient") %>%
  arrange(Patient)

## ---- 6) Write ---------------------------------------------------------------
write.csv(sym87,       file.path(PROC, "ml87_symptoms.csv"),   row.names = FALSE)
write.csv(prot87,      file.path(PROC, "ml87_proteomics.csv"), row.names = FALSE)
write.csv(combined87,  file.path(PROC, "ml87_combined.csv"),   row.names = FALSE)

## ---- 7) Report + sanity checks ----------------------------------------------
n_sym  <- ncol(sym87) - 2
n_prot <- ncol(prot87) - 2
n_comb <- ncol(combined87) - 2
message(sprintf("symptoms: 87 x %d | proteomics: %d x %d | combined: %d x %d",
                n_sym, nrow(prot87), n_prot, nrow(combined87), n_comb))
message("label balance (87): ",
        paste(names(table(label87$Lungcancer)), table(label87$Lungcancer),
              sep = "=", collapse = ", "))

stopifnot(
  nrow(sym87) == 87, nrow(prot87) == 87, nrow(combined87) == 87,
  n_comb == n_sym + n_prot,                                   # nothing lost in the merge
  identical(sym87$Patient, prot87$Patient),                  # SAME order for shared folds
  identical(sym87$Patient, combined87$Patient),
  all(c("No", "Yes") %in% as.character(label87$Lungcancer))
)
message("All 87-cohort table checks passed.")


# ================= build_ml_tables_symptbg.R =================

# =============================================================================
# build_ml_tables_symptbg.R  (Task 2: symptoms + background)
# -----------------------------------------------------------------------------
# The symptom model tables (symptoms_train.csv / symptoms_test87.csv) carry only
# the ~110 symptom features. Task 2 folds the BACKGROUND factors in as extra
# predictors: Age, Gender, and the Q* history/demographic block (living alone,
# education, prior flu/antibiotics, comorbidities Q7a-l, weight change, smoking).
# Those live in symptoms_clean.csv, so we join them on Patient.
#
# CRITICAL: the output keeps the SAME patients in the SAME row order as the
# symptom tables. That's what lets en_symptbg.R / rf_symptbg.R reuse the shared
# cv_folds.rds partitions (via retarget_folds) -- so "symptoms" vs
# "symptoms+background" is compared on identical CV splits, and any AUC change is
# due to the added features, not different folds.
#
# Writes: symptoms_bg_train.csv (411) and symptoms_bg_test87.csv (87).
# Run with working dir = LC_Data, after build_ml_tables.R.
# =============================================================================

library(tidyverse)

PROC  <- "processed_data"
clean <- read.csv(file.path(PROC, "symptoms_clean.csv"), check.names = FALSE)

# Background = Age, Gender, and every Q* column (symptoms are Br_/Co_/Ph_/Pa_/
# Fa_/Vo_/App_/Sm_/Fe_/Oth_, none start with Q, so this is a clean split).
bg_cols <- c("Age", "Gender", grep("^Q", names(clean), value = TRUE))
stopifnot(length(bg_cols) >= 3, !anyNA(clean[bg_cols]))
cat("Background factors added (", length(bg_cols), "):\n", paste(bg_cols, collapse = ", "), "\n\n")

add_background <- function(infile, outfile) {
  d   <- read.csv(file.path(PROC, infile), check.names = FALSE)
  out <- d %>% left_join(clean %>% dplyr::select(Patient, all_of(bg_cols)), by = "Patient")
  # guards: same rows, same order, no missing background
  stopifnot(nrow(out) == nrow(d),
            all(out$Patient == d$Patient),
            !anyNA(out[bg_cols]))
  write.csv(out, file.path(PROC, outfile), row.names = FALSE)
  cat(sprintf("%-22s -> %-25s %d rows x %d predictors (%d symptom + %d background)\n",
              infile, outfile, nrow(out), ncol(out) - 2, ncol(d) - 2, length(bg_cols)))
}

add_background("symptoms_train.csv",  "symptoms_bg_train.csv")
add_background("symptoms_test87.csv", "symptoms_bg_test87.csv")

message("\nDone. symptoms+background tables written -> feed into en_symptbg.R / rf_symptbg.R.")


# ================= build_ml_tables_LCP2.R =================

# =============================================================================
# build_ml_tables_LCP2.R
# -----------------------------------------------------------------------------
# Builds the proteomics-only modelling table for the LCP2 cohort:
#   ml_LCP2_proteomics.csv   Sample.ID, Patient.ID, LC, <surviving proteins> (prot_*)
#
# Input is the CLEANED + FILTERED long file from clean_LCP2_olink.R (S_1010
# already excluded, LOD-based protein filter already applied) -- this script
# only pivots it wide and attaches the label. It does not redo any QC.
#
# Mirrors build_ml_tables87.R's conventions (prot_ prefix, sanitised names,
# a Variable<->Assay lookup, sorted rows, sanity checks) but for a single
# proteomics-only cohort -- there's no symptoms/combined table for LCP2 yet
# (that's the later integration step, blocked on checking Patient ID overlap
# with the symptoms cohort -- see PROJECT_STATUS.md).
#
# Run with working dir = LC_Data, after clean_LCP2_olink.R.
# =============================================================================

library(tidyverse)

PROC <- "processed_data"

## ---- 1) Cleaned + filtered long data (from clean_LCP2_olink.R) ---------------
long <- read.delim(file.path(PROC, "LCP2_olink_clean_filtered_long.txt"),
                   sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("Input: %d rows | %d samples | %d proteins\n",
            nrow(long), dplyr::n_distinct(long$Sample.ID), dplyr::n_distinct(long$OlinkID)))

## ---- 2) Sanitise protein names + prefix "prot_"; keep a lookup ---------------
assay_map <- long %>% distinct(OlinkID, Assay)
assay_map$Variable <- paste0("prot_", make.names(gsub("[^A-Za-z0-9_]+", "_", assay_map$Assay),
                                                 unique = TRUE))
write.csv(assay_map, file.path(PROC, "prot_name_lookup_LCP2.csv"), row.names = FALSE)

## ---- 3) Pivot wide: rows = Sample.ID, cols = proteins ------------------------
# Pivot on OlinkID (guaranteed unique per protein) rather than Assay directly,
# then rename to the sanitised prot_ names via assay_map -- avoids any risk of
# duplicate column names if two assays ever shared a display name.
wide <- long %>%
  dplyr::select(Sample.ID, OlinkID, NPX) %>%
  pivot_wider(id_cols = Sample.ID, names_from = OlinkID, values_from = NPX)

id_cols_order <- names(wide)[-1]                              # OlinkIDs, in column order
new_names <- assay_map$Variable[match(id_cols_order, assay_map$OlinkID)]
names(wide)[-1] <- new_names

## ---- 4) Attach Patient.ID + LC label from the metadata -----------------------
meta <- read.csv("metadata/LCP2_metadata_Olink.csv", check.names = FALSE) %>%
  dplyr::select(Sample.ID, Patient.ID, LC) %>%
  mutate(LC = factor(LC, levels = c("Primary_LC", "No_Cancer")))   # cancer = positive, first

prot_LCP2 <- meta %>%
  inner_join(wide, by = "Sample.ID") %>%
  relocate(Sample.ID, Patient.ID, LC) %>%
  arrange(Sample.ID)

if (nrow(prot_LCP2) != nrow(wide))
  warning(sprintf("%d proteomics samples did not match a metadata row",
                  nrow(wide) - nrow(prot_LCP2)))

## ---- 5) Write ------------------------------------------------------------------
write.csv(prot_LCP2, file.path(PROC, "ml_LCP2_proteomics.csv"), row.names = FALSE)

## ---- 6) Report + sanity checks --------------------------------------------------
n_prot <- ncol(prot_LCP2) - 3   # minus Sample.ID, Patient.ID, LC
message(sprintf("ml_LCP2_proteomics.csv: %d samples x %d proteins", nrow(prot_LCP2), n_prot))
message("Label balance: ", paste(names(table(prot_LCP2$LC)), table(prot_LCP2$LC),
                                 sep = "=", collapse = ", "))

stopifnot(
  n_prot == dplyr::n_distinct(long$OlinkID),
  all(c("Primary_LC", "No_Cancer") %in% as.character(prot_LCP2$LC)),
  !anyNA(prot_LCP2$LC),
  nrow(prot_LCP2) == dplyr::n_distinct(long$Sample.ID)
)
message("All LCP2 table checks passed.")


# (CV folds are NOT built or saved here anymore. Each model script — 03/04/06/07 —
#  builds its folds inline with vfold_cv + a fixed seed, so nothing is written to
#  disk and no fold .rds carries the dataset around.)


# (87-cohort folds also built inline — in 06/07/08, v=5 repeats=10, same seed,
#  so all three 87 approaches share identical splits without any saved file.)


# (LCP2 folds also built inline — in 04, v=10 repeats=10, same seed.)

message("\nDone — all ML tables written to processed_data/. ",
        "CV folds are built inline in the model scripts (03/04/06/07), not saved to disk.")
