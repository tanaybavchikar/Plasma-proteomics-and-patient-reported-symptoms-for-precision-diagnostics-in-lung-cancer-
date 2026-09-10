# ============================================================
# 01_clean_data.R  — clean symptom, LCP1 (MS+Olink), LCP2 (Olink) data
# merged: clean_symptom_data.R + clean_data_LCP1.R + clean_LCP2_olink.R
# Run with working directory = LC_Data.
# ============================================================

# =============================================================================
# clean_symptom_data.R
# -----------------------------------------------------------------------------
# Produces SAVED, cleaned symptom-cohort files (the symptom analogue of
# clean_data_LCP1.R). Until now the exclusions only lived inside the analysis
# scripts; this writes them to disk so every downstream task reads the same
# cleaned data.
#
# Cleaning steps (from patients_to_exclude.docx):
#   1. Patient 1511 is KEPT as 1511 (the CORRECT id). Noora confirmed the MS data
#      mislabeled this patient as 1551; the symptom files already use 1511, so
#      nothing is relabeled here.
#   2. ALWAYS exclude: previous cancer (7) + incomplete symptoms (1207, not in file).
#   3. For STAGE analyses only: additionally exclude the 7 missing-staging patients.
#
# Outputs (to processed_data/):
#   symptoms_clean.csv        498 patients  (use for LC-vs-control work)
#   symptoms_clean_stage.csv  491 patients  (use only for stage-based analyses)
#
# Run with the working directory = LC_Data (where the Excel file lives).
# =============================================================================

library(tidyverse)
library(readxl)

if (!dir.exists("processed_data")) dir.create("processed_data")

# ---- Exclusion lists (subject/Patient IDs) ----------------------------------
excl_prev_cancer <- c(1080, 1538, 1600, 1606, 1618, 1703, 4195)  # always
excl_incomplete  <- c(1207)                                       # always (not in 505)
excl_no_stage    <- c(1603, 1792, 1854, 1855, 4207, 4302, 4378)   # stage analyses only

# ---- Read raw symptom file --------------------------------------------------
sym_raw <- read_excel("raw_data/505 pat bakgrund och nuvarande symtom.xlsx")

# ---- 1. Patient ID: keep 1511 (the correct id) ------------------------------
# Noora confirmed the MS data mislabeled this patient as 1551; the symptom files
# already use the correct 1511, so nothing is relabeled here.

# ---- 2. Main cleaned cohort: drop always-on exclusions ----------------------
symptoms_clean <- sym_raw %>%
  filter(!Patient %in% c(excl_prev_cancer, excl_incomplete))

# ---- 3. Stage-analysis cohort: additionally drop missing-staging ------------
symptoms_clean_stage <- symptoms_clean %>%
  filter(!Patient %in% excl_no_stage)

# ---- Save -------------------------------------------------------------------
write.csv(symptoms_clean,       "processed_data/symptoms_clean.csv",       row.names = FALSE)
write.csv(symptoms_clean_stage, "processed_data/symptoms_clean_stage.csv", row.names = FALSE)

message("Symptom cleaning done: ",
        nrow(symptoms_clean), " patients (main), ",
        nrow(symptoms_clean_stage), " patients (stage).")

# ---- Sanity checks (halt if counts drift) -----------------------------------
stopifnot(
  nrow(symptoms_clean)       == 498,
  nrow(symptoms_clean_stage) == 491,
  !any(c(excl_prev_cancer, excl_incomplete) %in% symptoms_clean$Patient),
  1511 %in% symptoms_clean$Patient,       # correct ID kept (NOT relabeled)
  !(1551 %in% symptoms_clean$Patient)     # mislabel must not appear on symptom side
)
message("All symptom sanity checks passed.")


# ================= from clean_data_LCP1.R =================

# =============================================================================
# clean_data_LCP1.R
# -----------------------------------------------------------------------------
# Adapted cleaning script for the LCP1 cohort, using ONLY the files present in
# the LC_Data folder (no external metadata / no "patients to exclude" file).
#
# MS columns are labeled by PATIENT ID (SampleID_1), mapped from the metadata
# file lc-metadata_MSsamples.txt (MS_ID == gene-table column name + "X__POOL_").
#
# Patient to exclude: 1618 (= LCP44), removed for previous cancer
# (patients_to_exclude.docx); 114 - 1 = 113 patients.
# Mislabel fix: Noora confirmed patient 1511 was mislabeled as 1551 in the MS
# metadata, so we correct 1551 -> 1511 here. (Olink uses LCP IDs, so it is
# unaffected and still excludes LCP44.)
#
# Run this with the LC_Data.Rproj open, so the working directory is LC_Data/.
# Outputs are written to LC_Data/processed_data/.
# =============================================================================

# ---- 0. Packages ------------------------------------------------------------
# Install anything missing, then load. (This covers the "update R packages" step.)
required <- c("tidyverse")   # tidyverse = dplyr, tidyr, stringr, readr, etc.
to_install <- required[!required %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)

library(tidyverse)

# Patients to exclude (previous cancer; the SAME person in two ID systems):
EXCLUDE_MS    <- "1618"    # patient ID (SampleID_1) used for the MS columns
EXCLUDE_OLINK <- "LCP44"   # LCP ID used in the Olink file

# Create output folder
if (!dir.exists("processed_data")) dir.create("processed_data")


# =============================================================================
# 1. MASS SPECTROMETRY (gene-centric) data
# =============================================================================

# Read the raw gene table (tab-separated). check.names = FALSE keeps the
# original column names (which contain the LCP IDs we need).
genes_raw <- read.delim("raw_data/LCP1_genes_table.txt", check.names = FALSE)

# The first 4 columns are annotation: Gene Name, Gene ID, Protein ID(s),
# Description. Keep these aside.
gene_ann <- genes_raw[, 1:4]

# Abundance columns are the TMT channels, e.g. "1562_LCP27_setA_tmt16plex_127N"
# (named PatientID_SampleID_TMTset_TMTtype_TMTtag, per the MS column dictionary).
#
# IMPORTANT: every channel ALSO has a sibling "... - Quanted PSM count" column
# that matches "tmtNplex" too. Those carry PSM counts, NOT abundance, so we must
# exclude them (along with the per-set PSM/Peptide/q-value/MS1-area summaries).
# Forgetting this silently averages counts into the abundances.
abundance_cols <- grep("tmt\\d+plex", colnames(genes_raw), value = TRUE)
abundance_cols <- abundance_cols[
  !str_detect(abundance_cols,
              regex("PSM|Peptide|Quanted|count|q.?value|area", ignore_case = TRUE))
]

# Drop the internal-standard channels (column name starts with "IS"); their
# values are all 0 because every sample is quantified relative to them.
sample_cols <- abundance_cols[!str_detect(abundance_cols, "^IS")]

# Build a tidy frame: Gene.Name + sample abundance columns only.
genes <- genes_raw %>%
  dplyr::select(Gene.Name = `Gene Name`, all_of(sample_cols))

# ---- Map each sample column to its PATIENT ID via the metadata --------------
# In the metadata, MS_ID equals the gene-table column name with an "X__POOL_"
# prefix, and SampleID_1 is the patient ID. So we build a lookup and relabel.
meta <- read.delim("metadata/lc-metadata_MSsamples.txt", check.names = FALSE)
meta$colkey <- sub("^X__POOL_", "", meta$MS_ID)
col_to_pat  <- setNames(as.character(meta$SampleID_1), meta$colkey)

pat_of_col <- col_to_pat[sample_cols]          # patient ID for each sample column
stopifnot(!anyNA(pat_of_col))                  # every column must map to a patient

# Correct Noora's mislabel: this patient is 1511, the MS metadata calls it 1551.
pat_of_col[pat_of_col == "1551"] <- "1511"

# ---- Merge duplicate (replicate) samples ------------------------------------
# 6 patients were measured twice across TMT sets. Average each patient's
# replicate channels (na.rm = TRUE), grouping columns by PATIENT ID.
unique_ids <- unique(pat_of_col)
merged_mat <- sapply(unique_ids, function(id) {
  cols <- sample_cols[pat_of_col == id]
  rowMeans(genes[, cols, drop = FALSE], na.rm = TRUE)
})
# rowMeans returns NaN when all values are NA -> convert back to NA.
merged_mat[is.nan(merged_mat)] <- NA

genes_merged <- bind_cols(
  dplyr::select(genes, Gene.Name),
  as.data.frame(merged_mat, check.names = FALSE)   # keep numeric IDs (no "X" prefix)
)

# Drop genes that are entirely NA across all samples.
genes_merged <- genes_merged %>%
  dplyr::filter(!if_all(-Gene.Name, is.na))
# Keep annotation rows in sync with the genes we retained.
gene_ann <- dplyr::filter(gene_ann, `Gene Name` %in% genes_merged$Gene.Name)

# ---- Exclude the flagged patient (1618 = previous cancer) -------------------
genes_clean <- genes_merged %>% dplyr::select(-any_of(EXCLUDE_MS))

# Final layout: genes in rows, samples in columns (already the case).
message("MS: ", nrow(genes_clean), " genes x ",
        ncol(genes_clean) - 1, " patients (expected 113).")

write.table(genes_clean, "processed_data/LCP1_MS_genecentric_clean.txt",
            sep = "\t", row.names = FALSE)
write.csv(gene_ann, "processed_data/LCP1_gene_annotations.csv", row.names = FALSE)


# =============================================================================
# 2. OLINK Explore (NPX) data
# =============================================================================

# The CSV is semicolon-separated and already in long NPX format
# (one row per sample x assay). We read it directly.
# (Alternative: OlinkAnalyze::read_NPX("raw_data/VB-3207_NPX_2022-09-15.csv"),
#  which gives the same long table plus Olink-specific helpers.)
olink <- read_delim("raw_data/VB-3207_NPX_2022-09-15.csv", delim = ";",
                    show_col_types = FALSE)

# ---- Drop technical control samples -----------------------------------------
# The plate controls (Control_1_Run122, Control_2_RUn123, ...) are not patients.
olink <- olink %>% dplyr::filter(!str_detect(SampleID, "^Control"))

# ---- Randomly pick ONE assay for proteins measured by several assays --------
# A few proteins (UniProt) are measured by more than one OlinkID (assay).
# We keep one assay per protein, chosen at random with a fixed seed for
# reproducibility. set.seed(1) reproduces the supervisor's selection:
#   OID31014, OID30225, OID30563, OID20074, OID20911, OID20153
olink_proteins   <- distinct(olink, OlinkID, UniProt)
duplicate_uniprot <- olink_proteins$UniProt[duplicated(olink_proteins$UniProt)]

duplicates <- list()
for (u in unique(duplicate_uniprot)) {
  duplicates[[u]] <- unique(olink_proteins$OlinkID[olink_proteins$UniProt == u])
}

set.seed(1)
keep <- c()
for (oids in duplicates) keep <- c(keep, sample(oids, size = 1))
remove <- unlist(duplicates)[!unlist(duplicates) %in% keep]

olink <- olink %>% dplyr::filter(!OlinkID %in% remove)

# ---- Exclude the flagged patient (LCP44) ------------------------------------
olink_clean <- olink %>% dplyr::filter(SampleID != EXCLUDE_OLINK)

message("Olink: ", n_distinct(olink_clean$UniProt), " proteins x ",
        n_distinct(olink_clean$SampleID),
        " patients (expected 2923 proteins x 87 patients).")

# Save the cleaned long table.
write.table(olink_clean, "processed_data/LCP1_olink_clean_long.txt",
            sep = "\t", row.names = FALSE)

# ---- Wide matrix: proteins in rows, samples in columns ----------------------
# This is the layout the agenda asks for. We use the Assay (gene symbol) as the
# row name; switch to UniProt if you prefer accessions.
olink_wide <- olink_clean %>%
  dplyr::select(Assay, SampleID, NPX) %>%
  pivot_wider(names_from = SampleID, values_from = NPX)

write.table(olink_wide, "processed_data/LCP1_olink_clean_wide.txt",
            sep = "\t", row.names = FALSE)


# =============================================================================
# 3. Sanity checks (fail loudly if the cohort sizes drift)
# =============================================================================
stopifnot(
  ncol(genes_clean) - 1 == 113,                       # MS patients (114 - excluded 1618)
  n_distinct(olink_clean$UniProt) == 2923,            # Olink proteins
  n_distinct(olink_clean$SampleID) == 87,             # Olink patients
  !(EXCLUDE_MS %in% colnames(genes_clean)),           # 1618 gone (MS)
  "1511" %in% colnames(genes_clean),                  # mislabel corrected to 1511
  !("1551" %in% colnames(genes_clean)),               # mislabel 1551 not present
  !(EXCLUDE_OLINK %in% olink_clean$SampleID)          # LCP44 gone (Olink)
)
message("All sanity checks passed.")

# ================= from clean_LCP2_olink.R =================

# =============================================================================
# clean_LCP2_olink.R
# -----------------------------------------------------------------------------
# Cleans the new LCP2 Olink Explore HT export (long format, 172 samples x 5416
# proteins, 2 plates) before it goes into build_ml_tables_LCP2.R.
#
# Three things happen here, in order:
#   1. Read the raw long file correctly (header + tab separator declared
#      explicitly; base R's row-name auto-detection then absorbs the extra
#      unnamed leading column the export has, since the header row has one
#      fewer field than the data rows).
#   2. Programmatically find and exclude the failed sample: per Sample.ID,
#      compute the fraction of that sample's SampleQC flags that are "FAIL",
#      and drop whichever sample sits far above everyone else (confirmed on
#      this file: S_1010, ~80.5% FAIL, every other sample 0%).
#   3. Compute, per protein, the proportion of (remaining) samples below LOD,
#      plot it, and filter out proteins over the PROT_MAX_BELOW_LOD cutoff.
#      APOE is a separate case: its LOD is NA in every sample (AssayQC WARN,
#      a protein-level issue unrelated to the S_1010 sample-level issue), so
#      its "proportion below LOD" is undefined (NaN), not zero. It's dropped
#      here as a precaution (DROP_UNDEFINED_LOD_PROTEINS) rather than kept and
#      quietly wrong.
#   4. (Optional) Restrict to proteins also measured in the LCP1 (87-cohort)
#      Olink panel, so the two cohorts share a common feature space (needed to
#      validate/replicate/combine an LCP2 model on LCP1). LCP1 and LCP2 are
#      DIFFERENT Olink products (Explore 3072-era vs Explore HT): their OlinkID
#      assay identifiers do NOT overlap at all, so the intersection is taken on
#      UniProt (the stable protein identifier), never OlinkID.
#
# Run with working dir = LC_Data.
# =============================================================================

library(tidyverse)

RAW  <- "raw_data/LCP2_olink.txt"
PROC <- "processed_data"
if (!dir.exists(PROC)) dir.create(PROC, recursive = TRUE)

## ---- tunables ----------------------------------------------------------------
SAMPLE_FAIL_PROP_CUTOFF   <- 0.5    # flag a sample if >50% of its assays FAILed
PROT_MAX_BELOW_LOD        <- 0.75   # keep proteins with <=75% of samples below LOD
DROP_UNDEFINED_LOD_PROTEINS <- TRUE # drop proteins where LOD is NA for every sample (e.g. APOE)
RESTRICT_TO_LCP1_OVERLAP  <- TRUE   # keep only proteins also measured in the LCP1 Olink panel
LCP1_REF <- "processed_data/LCP1_olink_clean_long.txt"  # source of the LCP1 protein list (has UniProt)
LCP1_MATCH_KEY <- "UniProt"         # join key across the two Olink products (NOT OlinkID -- see header)

## ---- 1) Read the raw long file ------------------------------------------------
LCP2_Olink <- read.table(RAW, header = TRUE, sep = "\t",
                         quote = "\"", stringsAsFactors = FALSE)
cat(sprintf("Read %d rows, %d columns from %s\n", nrow(LCP2_Olink), ncol(LCP2_Olink), RAW))
stopifnot(is.numeric(LCP2_Olink$NPX), is.numeric(LCP2_Olink$LOD))   # catches a broken read early

n_samples_raw  <- dplyr::n_distinct(LCP2_Olink$Sample.ID)
n_proteins_raw <- dplyr::n_distinct(LCP2_Olink$OlinkID)
cat(sprintf("Samples: %d | Proteins: %d\n", n_samples_raw, n_proteins_raw))

## ---- 2) Find and exclude the failed sample ------------------------------------
sample_qc_summary <- LCP2_Olink %>%
  group_by(Sample.ID) %>%
  summarise(
    n_assays  = n(),
    n_fail    = sum(SampleQC == "FAIL"),
    prop_fail = n_fail / n_assays,
    .groups = "drop"
  ) %>%
  arrange(desc(prop_fail))

cat("\nTop 5 samples by SampleQC FAIL proportion:\n")
print(head(sample_qc_summary, 5))

failed_samples <- sample_qc_summary %>%
  filter(prop_fail > SAMPLE_FAIL_PROP_CUTOFF) %>%
  pull(Sample.ID)
cat(sprintf("\nExcluding %d sample(s) as failed: %s\n",
            length(failed_samples), paste(failed_samples, collapse = ", ")))

cleaned_LCP2_Olink <- LCP2_Olink %>% filter(!(Sample.ID %in% failed_samples))
stopifnot(nrow(LCP2_Olink) - nrow(cleaned_LCP2_Olink) == length(failed_samples) * n_proteins_raw)

## ---- 3) Per-protein proportion below LOD --------------------------------------
detect_LOD <- function(data) {
  data %>%
    group_by(OlinkID, UniProt, Assay) %>%
    summarise(
      n_samples      = n(),
      n_na_lod       = sum(is.na(LOD)),
      n_valid        = n_samples - n_na_lod,
      n_below_lod    = sum(NPX < LOD, na.rm = TRUE),
      prop_below_lod = n_below_lod / n_valid,   # NaN if n_valid == 0 (e.g. APOE)
      .groups = "drop"
    )
}

proteins_below_LOD <- detect_LOD(cleaned_LCP2_Olink)

cat("\nProteins with an undefined (NaN) below-LOD proportion (LOD never computable):\n")
print(proteins_below_LOD %>% filter(n_na_lod > 0) %>% dplyr::select(OlinkID, Assay, n_na_lod))

## ---- Plot: proportion below LOD, all proteins --------------------------------
p_lod <- ggplot(proteins_below_LOD, aes(prop_below_lod)) +
  geom_histogram(binwidth = 0.03, boundary = 0) +
  geom_vline(xintercept = PROT_MAX_BELOW_LOD, color = "red", linetype = "dashed") +
  labs(title = "LCP2: proportion of samples below LOD, per protein",
       x = "Proportion of samples below LOD", y = "Number of proteins") +
  theme_minimal()
ggsave(file.path(PROC, "LCP2_prop_below_LOD_hist.png"), p_lod, width = 8, height = 5, dpi = 150)
print(p_lod)

## ---- Load the LCP1 protein list (for the cross-cohort overlap filter) ---------
# LCP1 and LCP2 are different Olink products, so we intersect on UniProt, not
# OlinkID (OlinkIDs do not overlap between products at all). Using the CLEANED
# LCP1 long table means the shared set is "usable in both cohorts", not just
# "measured in both".
lcp1_proteins <- character(0)
if (RESTRICT_TO_LCP1_OVERLAP) {
  lcp1_ref <- read.table(LCP1_REF, header = TRUE, sep = "\t",
                         quote = "\"", stringsAsFactors = FALSE)
  lcp1_proteins <- unique(lcp1_ref[[LCP1_MATCH_KEY]])
  n_shared <- sum(unique(proteins_below_LOD[[LCP1_MATCH_KEY]]) %in% lcp1_proteins)
  cat(sprintf("\nLCP1 overlap: %d LCP2 proteins share a %s with LCP1 (of %d LCP2 proteins).\n",
              n_shared, LCP1_MATCH_KEY, dplyr::n_distinct(proteins_below_LOD$OlinkID)))
}

## ---- Apply the filter, with an explicit reason per exclusion ------------------
# Order matters: the first matching condition wins, so a protein failing several
# checks is labelled by the most fundamental one. Putting the LOD checks before
# the overlap check makes "not in LCP1 panel" mean specifically "would have been
# kept on quality, but isn't shared" -- i.e. the marginal proteins lost to the
# cross-cohort restriction.
proteins_below_LOD <- proteins_below_LOD %>%
  mutate(
    reason_excluded = case_when(
      n_na_lod > 0 & DROP_UNDEFINED_LOD_PROTEINS               ~ "no computable LOD",
      prop_below_lod > PROT_MAX_BELOW_LOD                      ~ "too many below LOD",
      RESTRICT_TO_LCP1_OVERLAP & !(.data[[LCP1_MATCH_KEY]] %in% lcp1_proteins) ~ "not in LCP1 panel",
      TRUE                                                     ~ "keep"
    )
  )

cat("\nProtein exclusion breakdown:\n")
print(table(proteins_below_LOD$reason_excluded))

proteins_to_keep <- proteins_below_LOD %>%
  filter(reason_excluded == "keep") %>%
  pull(OlinkID)

final_LCP2_Olink <- cleaned_LCP2_Olink %>% filter(OlinkID %in% proteins_to_keep)

cat(sprintf("\nFinal: %d samples x %d proteins (of %d raw proteins)\n",
            dplyr::n_distinct(final_LCP2_Olink$Sample.ID),
            dplyr::n_distinct(final_LCP2_Olink$OlinkID),
            n_proteins_raw))

## ---- Write outputs -------------------------------------------------------------
write.csv(proteins_below_LOD, file.path(PROC, "LCP2_proteins_below_LOD.csv"), row.names = FALSE)
write.table(final_LCP2_Olink, file.path(PROC, "LCP2_olink_clean_filtered_long.txt"),
           sep = "\t", row.names = FALSE, quote = FALSE)

message("\nDone. Cleaned + filtered long table written to ",
        file.path(PROC, "LCP2_olink_clean_filtered_long.txt"),
        " -- feed this into build_ml_tables_LCP2.R next.")
