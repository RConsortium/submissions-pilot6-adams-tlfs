#************************************************************************
# Purpose:     Table 14.6.03 - Frequency of Normal and Abnormal (Clinically Significant Change from
#              Previous Visit) Laboratory Values During Treatment
# Input:       adlbc.json, adlbh.json, adsl.json
# Output:      data/tables/ard/t_14_6_03.rds (ARD), data/tables/pdf/t_14_6_03.pdf
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

## Load utility functions -------------------------------------------------
source(file.path("code", "utils", "doc_relative_path.r"))
source(file.path("code", "utils", "table_functions.r"))

## Load datasets --------------------------------------------------------
adlbc <- read_dataset_json(file.path(path$adam, "adlbc.json")) %>%
  rename_with(tolower)
adlbh <- read_dataset_json(file.path(path$adam, "adlbh.json")) %>%
  rename_with(tolower)
adsl <- read_dataset_json(file.path(path$adam, "adsl.json")) %>%
  rename_with(tolower)

# Note: anl01fl does not flag every visit in ADLBC, so scheduled-visit records are derived from visitnum/avisitn.
# Each record is compared with the preceding scheduled record for the same subject and parameter.

# Input checks -----------------------------------------------------------
# Validate required columns are present before derivation/rendering.
required_adsl <- c("saffl", "trt01a")
required_adlb <- c(
  "saffl", "ablfl", "ady", "adt", "avisitn", "visitnum", "usubjid", "trta", "param", "paramcd", "paramn",
  "parcat1", "aval", "a1lo", "a1hi"
)
required_adlbc <- required_adlb
required_adlbh <- required_adlb

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
# derive_lab_counts(df): returns list(counts, pvals) for one source dataset.
# pivot_wide(counts, pvals): reshapes counts into treatment-by-range columns and appends p-values.
# tidy_section(wide_df, section_label): normalizes columns for display and tags section labels.
# build_gt(data): renders the formatted gt table object used for PDF output.

# Derivations -----------------------------------------------------------

## Treatment N counts (Safety population) ------------------------------
trt_n <- adsl %>%
  filter(saffl == "Y") %>%
  count(trt01a, name = "n") %>%
  arrange(trt01a)

## Helper: derive per-arm Low/Normal/High counts + % + p-value ---------
# df: ADLBC or ADLBH. Per the SAP, a value is abnormal (Low/High) when it is beyond the threshold range
# (< 0.5 * LLN or > 1.5 * ULN) or differs from the preceding scheduled visit by more than 50% of the normal range
# (|change| > 0.5 * (ULN - LLN)); an increase is High and a decrease is Low.
derive_lab_counts <- function(df) {
  # Safety population, scheduled visits only (unscheduled and duplicate end-of-treatment records are excluded)
  filtered <- df[df$saffl == "Y" & df$visitnum == floor(df$visitnum) & !(df$avisitn %in% 99), ]

  # Flag each post-baseline on-treatment record against the preceding scheduled record (baseline/screening included)
  flagged <- filtered %>%
    group_by(usubjid, trta, param, paramcd, paramn, parcat1) %>%
    arrange(adt, visitnum, .by_group = TRUE) %>%
    mutate(
      chg_prev = aval - lag(aval),
      sig_chg = abs(chg_prev) > 0.5 * (a1hi - a1lo),
      is_low = (sig_chg & chg_prev < 0) | aval < 0.5 * a1lo,
      is_high = (sig_chg & chg_prev > 0) | aval > 1.5 * a1hi,
      is_low = dplyr::coalesce(is_low, FALSE),
      is_high = dplyr::coalesce(is_high, FALSE)
    ) %>%
    ungroup() %>%
    filter(!(ablfl %in% "Y"), (ady > 0) %in% TRUE)

  # Classify each subject by the first abnormal record (Low or High); subjects with none are Normal
  analysis <- flagged %>%
    group_by(usubjid, trta, param, paramcd, paramn, parcat1) %>%
    arrange(adt, visitnum, .by_group = TRUE) %>%
    slice(if (any(is_low | is_high)) which(is_low | is_high)[1] else 1L) %>%
    ungroup() %>%
    mutate(
      lbnrind_grp = case_when(
        is_low ~ "Low",
        is_high ~ "High",
        TRUE ~ "Normal"
      )
    )

  # Parameter-specific N per arm (subjects with at least one analysis record)
  # Used as denominator for percentages per the table specification.
  param_n <- analysis %>%
    count(param, trta, name = "n_param")

  # Counts per param x trt x range group
  counts <- analysis %>%
    count(param, paramcd, paramn, parcat1, trta, lbnrind_grp) %>%
    # Ensure all combinations exist (fill 0 for missing cells)
    complete(
      nesting(param, paramcd, paramn, parcat1),
      trta = unique(analysis$trta),
      lbnrind_grp = c("Low", "Normal", "High"),
      fill = list(n = 0L)
    ) %>%
    # Join parameter-specific N for percentage denominator
    left_join(param_n, by = c("param", "trta")) %>%
    mutate(
      proportion = ifelse(!is.na(n_param) & n_param > 0, n / n_param, 0),
      # Bare "0" for empty cells, padded with figure spaces so the 0 sits where the n digits end, just left of "("
      cell = ifelse(n == 0, paste0("0", strrep(" ", 5)), sprintf("%d(%s%%)", n, format_percent(proportion)))
    )

  # Fisher exact p-value per param (Low/Normal/High, across arms)
  pvals <- analysis %>%
    group_by(param) %>%
    summarise(
      p_val = tryCatch(
        {
          tbl <- table(trta, lbnrind_grp)
          if (nrow(tbl) < 2 || ncol(tbl) < 2) {
            NA_real_
          } else {
            fisher.test(tbl, workspace = 2e7)$p.value
          }
        },
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(p_fmt = ifelse(is.na(p_val), "", sprintf("%.3f", p_val)))

  list(counts = counts, pvals = pvals)
}

## Run derivations for CHEM and HEM ------------------------------------
chem <- derive_lab_counts(adlbc)
hem <- derive_lab_counts(adlbh)

## Pivot wide: one row per param, columns = trt x lbnrind_grp ----------
pivot_wide <- function(counts, pvals) {
  overall_freq <- counts %>%
    group_by(param, paramcd, paramn, parcat1) %>%
    summarise(overall_n = sum(n), .groups = "drop")

  counts %>%
    select(param, paramcd, paramn, parcat1, trta, lbnrind_grp, cell) %>%
    pivot_wider(
      names_from = c(trta, lbnrind_grp),
      values_from = cell,
      names_sep = "||"
    ) %>%
    left_join(overall_freq, by = c("param", "paramcd", "paramn", "parcat1")) %>%
    left_join(pvals, by = "param") %>%
    arrange(paramcd)
}

chem_wide <- pivot_wide(chem$counts, chem$pvals)
hem_wide <- pivot_wide(hem$counts, hem$pvals)

## Build display data ---------------------------------------------------
# Expected treatment arm order
trt_order <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
lnr_order <- c("Low", "Normal", "High")

col_order <- c(
  "param",
  paste(rep(trt_order, each = 3), lnr_order, sep = "||"),
  "p_fmt"
)

tidy_section <- function(wide_df, section_label) {
  # Select only columns that exist
  existing_cols <- intersect(c(col_order, "overall_n"), names(wide_df))
  wide_df <- wide_df %>%
    select(param, paramcd, parcat1, paramn, all_of(setdiff(existing_cols, "param")))

  # Fill any missing trt x lbnrind combos with blank
  for (col in setdiff(col_order[-1], names(wide_df))) {
    wide_df[[col]] <- ""
  }

  wide_df %>%
    select(all_of(c("param", "paramcd", "parcat1", "paramn", "overall_n", col_order[-1]))) %>%
    mutate(section = section_label) %>%
    arrange(paramcd) %>%
    select(-overall_n)
}

display_data <- bind_rows(
  tidy_section(chem_wide, "CHEMISTRY"),
  tidy_section(hem_wide, "HEMATOLOGY")
)

# ARD output (intermediate) --------------------------------------------
# Build a lightweight gtsummary representation from display_data solely for
# ARD serialization; PDF rendering continues to use the custom gt layout below.
t_14_6_03_summary_tbl <- display_data %>%
  mutate(section = factor(section, levels = c("CHEMISTRY", "HEMATOLOGY"))) %>%
  tbl_summary(
    by = section,
    include = param,
    statistic = list(all_categorical() ~ "{n}"),
    missing = "no"
  ) %>%
  modify_header(label = "")

t_14_6_03_ard <- gather_ard(t_14_6_03_summary_tbl)
ard_tbl <- if (is.list(t_14_6_03_ard) && "tbl_summary" %in% names(t_14_6_03_ard)) {
  t_14_6_03_ard$tbl_summary
} else {
  t_14_6_03_ard
}

if (!is.data.frame(ard_tbl) || nrow(ard_tbl) == 0) {
  stop(
    "ARD output is empty for t_14_6_03. Check ADLBC/ADLBH analysis records and display_data derivation filters."
  )
}

saveRDS(ard_tbl, file.path(path$table_ard, "t_14_6_03.rds"))

# Table rendering with gt -----------------------------------------------

## Column label helpers -------------------------------------------------
# Build N labels for header
n_label <- function(trt_name) {
  n <- trt_n$n[trt_n$trt01a == trt_name]
  if (length(n) == 0) {
    trt_name
  } else {
    sprintf("%s (N=%d)", trt_name, n)
  }
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
  summarise(first_paramn = min(paramn), .groups = "drop") %>%
  arrange(first_paramn)

## gt table build -------------------------------------------------------
build_gt <- function(data) {
  trt_col_map <- setNames(col_order[-length(col_order)], col_order[-length(col_order)])

  gt_tbl <- data %>%
    select(-paramcd, -parcat1, -paramn, -section) %>%
    gt() %>%
    tab_header(
      title = md("**Table 14-6.03**"),
      subtitle = md(
        "**Frequency of Normal and Abnormal (Clinically Significant Change from Previous Visit)**  \n**Laboratory Values During Treatment**"
      )
    ) %>%
    cols_label(param = "") %>%
    cols_label(
      !!paste("Placebo", "Low", sep = "||") := "Low",
      !!paste("Placebo", "Normal", sep = "||") := "Normal",
      !!paste("Placebo", "High", sep = "||") := "High",
      !!paste("Xanomeline Low Dose", "Low", sep = "||") := "Low",
      !!paste("Xanomeline Low Dose", "Normal", sep = "||") := "Normal",
      !!paste("Xanomeline Low Dose", "High", sep = "||") := "High",
      !!paste("Xanomeline High Dose", "Low", sep = "||") := "Low",
      !!paste("Xanomeline High Dose", "Normal", sep = "||") := "Normal",
      !!paste("Xanomeline High Dose", "High", sep = "||") := "High",
      p_fmt = "p-val\n[1]"
    ) %>%
    tab_spanner(
      label = sprintf("Placebo (N=%d)", trt_n$n[trt_n$trt01a == "Placebo"]),
      columns = starts_with("Placebo||")
    ) %>%
    tab_spanner(
      label = sprintf(
        "Xan. Low (N=%d)",
        trt_n$n[trt_n$trt01a == "Xanomeline Low Dose"]
      ),
      columns = starts_with("Xanomeline Low Dose||")
    ) %>%
    tab_spanner(
      label = sprintf(
        "Xan. High (N=%d)",
        trt_n$n[trt_n$trt01a == "Xanomeline High Dose"]
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
    cols_align(align = "right", columns = -param) %>%
    cols_align(align = "left", columns = param) %>%
    cols_width(
      param ~ px(200),
      starts_with("Placebo||") ~ px(80),
      starts_with("Xanomeline") ~ px(80)
    ) %>%
    tab_options(
      table.font.size = px(7),
      data_row.padding = px(0.5),
      row_group.font.weight = "bold",
      table.border.top.style = "hidden",
      table.border.bottom.style = "hidden",
      table.border.left.style = "hidden",
      table.border.right.style = "hidden",
      stub.border.style = "none",
      table_body.vlines.style = "none",
      column_labels.vlines.style = "none",
      row_group.border.left.style = "none",
      row_group.border.right.style = "none"
    )

  gt_tbl
}

gt_tbl <- build_gt(display_data)

## Save PDF output ------------------------------------------------------
gt_tbl %>%
  as_docorator(
    display_name = "t_14_6_03",
    display_loc = path$table_output,
    save_object = FALSE,
    header = fancyhead(
      fancyrow(left = "Protocol: CDISCPILOT01", center = NA, right = doc_pagenum()),
      fancyrow(left = "Population: Safety", center = NA, right = NA)
    ),
    footer = fancyfoot(
      fancyrow(left = "[1] Fisher exact p-value (treatment arm comparison across low, normal and high categories)."),
      fancyrow(left = "Note: Percentages are based on the number of subjects with non-missing assessments (i.e., the total of the"),
      fancyrow(left = "subjects in the low, normal, and high categories) within each treatment group."),
      fancyrow(left = paste0("Source: ", doc_relative_path()), center = NA, right = doc_datetime())
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
  render_pdf(path$table_output)

message(
  "Table 14.6.03 saved: ",
  file.path(path$table_ard, "t_14_6_03.rds"),
  " and ",
  file.path(path$table_output, "t_14_6_03.pdf")
)
