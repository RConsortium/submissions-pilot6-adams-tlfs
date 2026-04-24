#************************************************************************
# Purpose:     Table 14.6.02 - Frequency of Normal and Abnormal (Beyond
#              Normal Range) Laboratory Values During Treatment
# Input:       adlbc.json, adlbh.json, adsl.json
# Output:      data/tables/t-14-6-02.rds (ARD), data/tables/t-14-6-02.pdf
# Population:  Safety (SAFFL = "Y")
# Filter:      On-treatment records (AENTMTFL = "Y")
#************************************************************************

# Setup -----------------------------------------------------------------
## Load libraries -------------------------------------------------------
library(dplyr)
library(tidyr)
library(datasetjson)
library(gtsummary)
library(flextable)

## Load datasets --------------------------------------------------------
adlbc <- read_dataset_json(file.path(path$adam, "adlbc.json"), decimals_as_floats = TRUE)
adlbh <- read_dataset_json(file.path(path$adam, "adlbh.json"), decimals_as_floats = TRUE)
adsl  <- read_dataset_json(file.path(path$adam, "adsl.json"),  decimals_as_floats = TRUE)

# Note: ADLBC uses ANL01FL = "Y" and ADLBH uses AENTMTFL = "Y" as analysis flags.

# Input checks -----------------------------------------------------------
# Validate required columns are present before derivation/rendering.
required_adsl <- c("SAFFL", "TRT01A")
required_adlbc <- c("SAFFL", "ANL01FL", "USUBJID", "TRTA", "PARAM", "PARAMN", "PARCAT1", "LBNRIND")
required_adlbh <- c("SAFFL", "AENTMTFL", "USUBJID", "TRTA", "PARAM", "PARAMN", "PARCAT1", "LBNRIND")

missing_adsl <- setdiff(required_adsl, names(adsl))
missing_adlbc <- setdiff(required_adlbc, names(adlbc))
missing_adlbh <- setdiff(required_adlbh, names(adlbh))

if (length(missing_adsl) > 0) {
  stop(sprintf("ADSL is missing required columns: %s", paste(missing_adsl, collapse = ", ")))
}
if (length(missing_adlbc) > 0) {
  stop(sprintf("ADLBC is missing required columns: %s", paste(missing_adlbc, collapse = ", ")))
}
if (length(missing_adlbh) > 0) {
  stop(sprintf("ADLBH is missing required columns: %s", paste(missing_adlbh, collapse = ", ")))
}

# Helper function overview ----------------------------------------------
# derive_lab_counts(df, anl_flag): returns list(counts, pvals) for one source dataset.
# pivot_wide(counts, pvals): reshapes counts into treatment-by-range columns and appends p-values.
# tidy_section(wide_df, section_label): normalizes columns for display and tags section labels.
# build_gt(data): renders the formatted gt table object used for PDF output.

# Derivations -----------------------------------------------------------

## Treatment N counts (Safety population) ------------------------------
trt_n <- adsl %>%
  filter(SAFFL == "Y") %>%
  count(TRT01A, name = "N") %>%
  arrange(TRT01A)

## Helper: derive per-arm Low/Normal/High counts + % + p-value ---------
# anl_flag: character name of analysis flag column
derive_lab_counts <- function(df, anl_flag = "ANL01FL") {
  # Subset to analysis population: safety + analysis flag
  filtered <- df[df$SAFFL == "Y" & !is.na(df[[anl_flag]]) & df[[anl_flag]] == "Y", ]

  # Worst post-baseline record per subject per parameter:
  # Low (1) > High (2) > Normal (3)
  analysis <- filtered %>%
    mutate(
      worst_rank = case_when(
        LBNRIND == "LOW"  ~ 1L,
        LBNRIND == "HIGH" ~ 2L,
        TRUE              ~ 3L
      )
    ) %>%
    group_by(USUBJID, TRTA, PARAM, PARAMN, PARCAT1) %>%
    slice_min(worst_rank, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      LBNRIND_grp = case_when(
        LBNRIND == "LOW"  ~ "Low",
        LBNRIND == "HIGH" ~ "High",
        TRUE              ~ "Normal"
      )
    )

  # Parameter-specific N per arm (subjects with at least one analysis record)
  # Used as denominator for percentages per the table specification.
  param_n <- analysis %>%
    count(PARAM, TRTA, name = "N_param")

  # Counts per PARAM x TRT x range group
  counts <- analysis %>%
    count(PARAM, PARAMN, PARCAT1, TRTA, LBNRIND_grp) %>%
    # Ensure all combinations exist (fill 0 for missing cells)
    complete(
      nesting(PARAM, PARAMN, PARCAT1),
      TRTA = unique(analysis$TRTA),
      LBNRIND_grp = c("Low", "Normal", "High"),
      fill = list(n = 0L)
    ) %>%
    # Join parameter-specific N for percentage denominator
    left_join(param_n, by = c("PARAM", "TRTA")) %>%
    mutate(
      pct  = ifelse(!is.na(N_param) & N_param > 0, n / N_param * 100, 0),
      cell = sprintf("%d(%d%%)", n, round(pct))
    )

  # Chi-square p-value per PARAM (Low+High vs Normal, across arms)
  pvals <- analysis %>%
    mutate(abnormal = ifelse(LBNRIND_grp %in% c("Low", "High"), "Abnormal", "Normal")) %>%
    group_by(PARAM) %>%
    summarise(
      p_val = tryCatch(
        {
          tbl <- table(TRTA, abnormal)
          if (nrow(tbl) < 2 || ncol(tbl) < 2) NA_real_
          else chisq.test(tbl, simulate.p.value = TRUE)$p.value
        },
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(p_fmt = ifelse(is.na(p_val), "", sprintf("%.3f", p_val)))

  list(counts = counts, pvals = pvals)
}

## Run derivations for CHEM and HEM ------------------------------------
chem <- derive_lab_counts(adlbc, anl_flag = "ANL01FL")
hem  <- derive_lab_counts(adlbh, anl_flag = "AENTMTFL")

## Pivot wide: one row per PARAM, columns = TRT x LBNRIND_grp ----------
pivot_wide <- function(counts, pvals) {
  counts %>%
    select(PARAM, PARAMN, PARCAT1, TRTA, LBNRIND_grp, cell) %>%
    pivot_wider(
      names_from  = c(TRTA, LBNRIND_grp),
      values_from = cell,
      names_sep   = "||"
    ) %>%
    left_join(pvals, by = "PARAM") %>%
    arrange(PARAMN)
}

chem_wide <- pivot_wide(chem$counts, chem$pvals)
hem_wide  <- pivot_wide(hem$counts,  hem$pvals)

## Build display data ---------------------------------------------------
# Expected treatment arm order
trt_order <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
lnr_order <- c("Low", "Normal", "High")

col_order <- c(
  "PARAM",
  paste(rep(trt_order, each = 3), lnr_order, sep = "||"),
  "p_fmt"
)

tidy_section <- function(wide_df, section_label) {
  # Select only columns that exist
  existing_cols <- intersect(col_order, names(wide_df))
  wide_df <- wide_df %>%
    select(PARAM, PARCAT1, PARAMN, all_of(setdiff(existing_cols, "PARAM")))

  # Fill any missing TRT x LBNRIND combos with blank
  for (col in setdiff(col_order[-1], names(wide_df))) {
    wide_df[[col]] <- ""
  }

  wide_df %>%
    select(all_of(c("PARAM", "PARCAT1", "PARAMN", col_order[-1]))) %>%
    mutate(section = section_label) %>%
    arrange(desc(PARAMN))
}

display_data <- bind_rows(
  tidy_section(chem_wide, "CHEMISTRY"),
  tidy_section(hem_wide,  "HEMATOLOGY")
)

# Output locations ------------------------------------------------------
table_output <- if (is.null(path$table_output)) {
  file.path(getwd(), "data/tables/pdf")
} else {
  path$table_output
}

table_ard <- if (is.null(path$table_ard)) {
  file.path(getwd(), "data/tables/ard")
} else {
  path$table_ard
}

dir.create(table_output, recursive = TRUE, showWarnings = FALSE)
dir.create(table_ard, recursive = TRUE, showWarnings = FALSE)

# ARD output (intermediate) --------------------------------------------
saveRDS(display_data, file.path(table_ard, "t-14-6-02.rds"))


# Table rendering with gtsummary ---------------------------------------

build_gtsummary <- function(data) {
  # Remove grouping/section columns for display
  display <- data %>% select(-PARCAT1, -PARAMN, -section)
  # gtsummary expects a data.frame with variables as columns
  tbl <- gtsummary::tbl_summary(
    display,
    by = NULL,
    statistic = list(all_categorical() ~ "{n}"),
    label = list(PARAM ~ "Parameter", p_fmt ~ "p-val [1]")
  )
  tbl
}

gtsum_tbl <- build_gtsummary(display_data)

# Save as RDS (ARD output)
saveRDS(display_data, file.path(table_ard, "t-14-6-02.rds"))

# Save as PDF using flextable (via as_flex_table)
ft <- as_flex_table(gtsum_tbl)
pdf_file <- file.path(table_output, "t-14-6-02.pdf")
tmp_img <- tempfile(fileext = ".png")
flextable::save_as_image(ft, path = tmp_img)
if (requireNamespace("magick", quietly = TRUE)) {
  img <- magick::image_read(tmp_img)
  magick::image_write(img, path = pdf_file, format = "pdf")
} else {
  warning("magick package required to save PDF. Saved as PNG instead.")
  file.copy(tmp_img, sub(".pdf$", ".png", pdf_file))
}

message(
  "Table 14.6.02 saved: ",
  file.path(table_ard, "t-14-6-02.rds"),
  " and ",
  pdf_file
)
