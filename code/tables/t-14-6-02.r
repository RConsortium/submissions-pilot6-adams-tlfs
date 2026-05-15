#************************************************************************
# Purpose:     Table 14.6.02 - Frequency of Normal and Abnormal (Beyond
#              Normal Range) Laboratory Values During Treatment
# Input:       adlbc.json, adlbh.json, adsl.json
# Output:      data/tables/t-14-6-02.rds (ARD), data/tables/t-14-6-02.pdf
# Population:  Safety population
# Filter:      On-treatment analysis records
#************************************************************************

# Setup -----------------------------------------------------------------
## Load libraries -------------------------------------------------------
library(dplyr)
library(tidyr)
library(datasetjson)
library(cards)
library(gtsummary)
library(gt)
library(docorator)

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
      pct = ifelse(!is.na(N_param) & N_param > 0, n / N_param * 100, 0),
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
  overall_freq <- counts %>%
    group_by(PARAM, PARAMN, PARCAT1) %>%
    summarise(overall_n = sum(n), .groups = "drop")

  counts %>%
    select(PARAM, PARAMN, PARCAT1, TRTA, LBNRIND_grp, cell) %>%
    pivot_wider(
      names_from = c(TRTA, LBNRIND_grp),
      values_from = cell,
      names_sep = "||"
    ) %>%
    left_join(overall_freq, by = c("PARAM", "PARAMN", "PARCAT1")) %>%
    left_join(pvals, by = "PARAM") %>%
    arrange(desc(overall_n), PARAMN)
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
    select(all_of(c("PARAM", "PARCAT1", "PARAMN", "overall_n", col_order[-1]))) %>%
    mutate(section = section_label) %>%
    arrange(desc(overall_n), desc(PARAMN)) %>%
    select(-overall_n)
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
# Build a lightweight gtsummary representation from display_data solely for
# ARD serialization; PDF rendering continues to use the custom gt layout below.
t_14_6_02_summary_tbl <- display_data %>%
  mutate(section = factor(section, levels = c("CHEMISTRY", "HEMATOLOGY"))) %>%
  tbl_summary(
    by = section,
    include = PARAM,
    statistic = list(all_categorical() ~ "{n}"),
    missing = "no"
  ) %>%
  modify_header(label = "")

t_14_6_02_ard <- gather_ard(t_14_6_02_summary_tbl)
if (nrow(t_14_6_02_ard) == 0) {
  stop(
    "ARD output is empty for t-14-6-02. Check ADLBC/ADLBH analysis records and display_data derivation filters."
  )
}

saveRDS(t_14_6_02_ard, file.path(table_ard, "t-14-6-02.rds"))

# Table rendering with gt -----------------------------------------------

## Column label helpers -------------------------------------------------
# Build N labels for header
n_label <- function(trt_name) {
  n <- trt_n$N[trt_n$TRT01A == trt_name]
  if (length(n) == 0) trt_name
  else sprintf("%s (N=%d)", trt_name, n)
}

trt_labels <- setNames(
  lapply(trt_order, n_label),
  trt_order
)

# Abbreviate treatment names for display
trt_display <- c(
  "Placebo"              = "Placebo",
  "Xanomeline Low Dose"  = "Xan. Low",
  "Xanomeline High Dose" = "Xan. High"
)

## Build section rows (blank separator + section header) ----------------
section_rows <- display_data %>%
  group_by(section) %>%
  summarise(first_paramn = min(PARAMN), .groups = "drop") %>%
  arrange(first_paramn)

## gt table build -------------------------------------------------------
build_gt <- function(data) {
  trt_col_map <- setNames(col_order[-length(col_order)], col_order[-length(col_order)])

  gt_tbl <- data %>%
    select(-PARCAT1, -PARAMN, -section) %>%
    gt(rowname_col = "PARAM") %>%
    tab_header(
      title = md("**Table 14-6.02**"),
      subtitle = md(
        "**Frequency of Normal and Abnormal (Beyond Normal Range)**  \n**Laboratory Values During Treatment**"
      )
    ) %>%
    tab_stubhead(label = "") %>%
    cols_label(
      !!paste("Placebo",              "Low",    sep = "||") := "Low",
      !!paste("Placebo",              "Normal", sep = "||") := "Normal",
      !!paste("Placebo",              "High",   sep = "||") := "High",
      !!paste("Xanomeline Low Dose",  "Low",    sep = "||") := "Low",
      !!paste("Xanomeline Low Dose",  "Normal", sep = "||") := "Normal",
      !!paste("Xanomeline Low Dose",  "High",   sep = "||") := "High",
      !!paste("Xanomeline High Dose", "Low",    sep = "||") := "Low",
      !!paste("Xanomeline High Dose", "Normal", sep = "||") := "Normal",
      !!paste("Xanomeline High Dose", "High",   sep = "||") := "High",
      p_fmt = "p-val\n[1]"
    ) %>%
    tab_spanner(
      label = sprintf("Placebo (N=%d)", trt_n$N[trt_n$TRT01A == "Placebo"]),
      columns = starts_with("Placebo||")
    ) %>%
    tab_spanner(
      label = sprintf(
        "Xan. Low (N=%d)",
        trt_n$N[trt_n$TRT01A == "Xanomeline Low Dose"]
      ),
      columns = starts_with("Xanomeline Low Dose||")
    ) %>%
    tab_spanner(
      label = sprintf(
        "Xan. High (N=%d)",
        trt_n$N[trt_n$TRT01A == "Xanomeline High Dose"]
      ),
      columns = starts_with("Xanomeline High Dose||")
    ) %>%
    tab_row_group(
      label = md("**HEMATOLOGY**  \n----------"),
      rows = data$section == "HEMATOLOGY"
    ) %>%
    tab_row_group(
      label = md("**CHEMISTRY**  \n----------"),
      rows = data$section == "CHEMISTRY"
    ) %>%
    opt_table_font(font = "Courier New") %>%
    cols_align(align = "right", columns = -PARAM) %>%
    cols_width(
      PARAM ~ px(200),
      starts_with("Placebo||") ~ px(80),
      starts_with("Xanomeline") ~ px(80)
    ) %>%
    tab_options(
      table.font.size = px(7),
      data_row.padding = px(0.5),
      row_group.font.weight = "bold",
      table.border.top.style = "hidden",
      table.border.bottom.style = "hidden"
    )

  gt_tbl
}

gt_tbl <- build_gt(display_data)

## Save PDF output ------------------------------------------------------
gt_tbl %>%
  as_docorator(
    display_name = "t-14-6-02",
    display_loc = table_output,
    save_object = FALSE,
    header = fancyhead(
      fancyrow(left = "Protocol: CDISCPILOT01", center = NA, right = doc_pagenum()),
      fancyrow(left = "Population: Safety", center = NA, right = NA)
    ),
    footer = fancyfoot(
      fancyrow(left = "[1] Chi-square p-value (treatment arm comparison of abnormal vs normal)."),
      fancyrow(left = doc_path(), center = NA, right = doc_datetime())
    ),
    geometry = geom_set(
      paperwidth = "11in",
      paperheight = "8.5in",
      left = "0.5in",
      right = "0.5in",
      top = "0.5in",
      bottom = "0.5in",
      headheight = "24pt",
      footskip = "24pt"
    )
  ) %>%
  render_pdf(table_output)

message(
  "Table 14.6.02 saved: ",
  file.path(table_ard, "t-14-6-02.rds"),
  " and ",
  file.path(table_output, "t-14-6-02.pdf")
)
