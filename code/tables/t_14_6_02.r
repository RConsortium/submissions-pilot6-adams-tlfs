#************************************************************************
# Purpose:     Table 14.6.02 - Frequency of Normal and Abnormal (Beyond
#              Normal Range) Laboratory Values During Treatment
# Input:       adlbc.json, adlbh.json, adsl.json
# Output:      data/tables/ard/t_14_6_02.rds (ARD), data/tables/pdf/t_14_6_02.pdf
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
adlbc <- read_dataset_json(file.path(path$adam, "adlbc.json"), decimals_as_floats = TRUE) %>%
  rename_with(tolower)
adlbh <- read_dataset_json(file.path(path$adam, "adlbh.json"), decimals_as_floats = TRUE) %>%
  rename_with(tolower)
adsl <- read_dataset_json(file.path(path$adam, "adsl.json"), decimals_as_floats = TRUE) %>%
  rename_with(tolower)

# Note: ADLBC uses anl01fl = "Y" and ADLBH uses aentmtfl = "Y" as analysis flags.

# Input checks -----------------------------------------------------------
# Validate required columns are present before derivation/rendering.
required_adsl <- c("saffl", "trt01a")
required_adlbc <- c("saffl", "anl01fl", "usubjid", "trta", "param", "paramn", "parcat1", "lbnrind")
required_adlbh <- c("saffl", "aentmtfl", "usubjid", "trta", "param", "paramn", "parcat1", "lbnrind")

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
  filter(saffl == "Y") %>%
  count(trt01a, name = "n") %>%
  arrange(trt01a)

## Helper: derive per-arm Low/Normal/High counts + % + p-value ---------
# anl_flag: character name of analysis flag column
derive_lab_counts <- function(df, anl_flag = "anl01fl") {
  # Subset to analysis population: safety + analysis flag
  filtered <- df[df$saffl == "Y" & !is.na(df[[anl_flag]]) & df[[anl_flag]] == "Y", ]

  # Worst post-baseline record per subject per parameter:
  # Low (1) > High (2) > Normal (3)
  analysis <- filtered %>%
    mutate(
      worst_rank = case_when(
        lbnrind == "LOW" ~ 1L,
        lbnrind == "HIGH" ~ 2L,
        TRUE ~ 3L
      )
    ) %>%
    group_by(usubjid, trta, param, paramn, parcat1) %>%
    slice_min(worst_rank, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      lbnrind_grp = case_when(
        lbnrind == "LOW" ~ "Low",
        lbnrind == "HIGH" ~ "High",
        TRUE ~ "Normal"
      )
    )

  # Parameter-specific N per arm (subjects with at least one analysis record)
  # Used as denominator for percentages per the table specification.
  param_n <- analysis %>%
    count(param, trta, name = "n_param")

  # Counts per param x trt x range group
  counts <- analysis %>%
    count(param, paramn, parcat1, trta, lbnrind_grp) %>%
    # Ensure all combinations exist (fill 0 for missing cells)
    complete(
      nesting(param, paramn, parcat1),
      trta = unique(analysis$trta),
      lbnrind_grp = c("Low", "Normal", "High"),
      fill = list(n = 0L)
    ) %>%
    # Join parameter-specific N for percentage denominator
    left_join(param_n, by = c("param", "trta")) %>%
    mutate(
      proportion = ifelse(!is.na(n_param) & n_param > 0, n / n_param, 0),
      cell = sprintf("%d(%s%%)", n, format_percent(proportion))
    )

  # Chi-square p-value per param (Low+High vs Normal, across arms)
  pvals <- analysis %>%
    mutate(abnormal = ifelse(lbnrind_grp %in% c("Low", "High"), "Abnormal", "Normal")) %>%
    group_by(param) %>%
    summarise(
      p_val = tryCatch(
        {
          tbl <- table(trta, abnormal)
          if (nrow(tbl) < 2 || ncol(tbl) < 2) {
            NA_real_
          } else {
            chisq.test(tbl, simulate.p.value = TRUE)$p.value
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
chem <- derive_lab_counts(adlbc, anl_flag = "anl01fl")
hem <- derive_lab_counts(adlbh, anl_flag = "aentmtfl")

## Pivot wide: one row per param, columns = trt x lbnrind_grp ----------
pivot_wide <- function(counts, pvals) {
  overall_freq <- counts %>%
    group_by(param, paramn, parcat1) %>%
    summarise(overall_n = sum(n), .groups = "drop")

  counts %>%
    select(param, paramn, parcat1, trta, lbnrind_grp, cell) %>%
    pivot_wider(
      names_from = c(trta, lbnrind_grp),
      values_from = cell,
      names_sep = "||"
    ) %>%
    left_join(overall_freq, by = c("param", "paramn", "parcat1")) %>%
    left_join(pvals, by = "param") %>%
    arrange(desc(overall_n), paramn)
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
    select(param, parcat1, paramn, all_of(setdiff(existing_cols, "param")))

  # Fill any missing trt x lbnrind combos with blank
  for (col in setdiff(col_order[-1], names(wide_df))) {
    wide_df[[col]] <- ""
  }

  wide_df %>%
    select(all_of(c("param", "parcat1", "paramn", "overall_n", col_order[-1]))) %>%
    mutate(section = section_label) %>%
    arrange(desc(overall_n), desc(paramn)) %>%
    select(-overall_n)
}

display_data <- bind_rows(
  tidy_section(chem_wide, "CHEMISTRY"),
  tidy_section(hem_wide, "HEMATOLOGY")
)

# ARD output (intermediate) --------------------------------------------
# Build a lightweight gtsummary representation from display_data solely for
# ARD serialization; PDF rendering continues to use the custom gt layout below.
t_14_6_02_summary_tbl <- display_data %>%
  mutate(section = factor(section, levels = c("CHEMISTRY", "HEMATOLOGY"))) %>%
  tbl_summary(
    by = section,
    include = param,
    statistic = list(all_categorical() ~ "{n}"),
    missing = "no"
  ) %>%
  modify_header(label = "")

t_14_6_02_ard <- gather_ard(t_14_6_02_summary_tbl)
ard_tbl <- if (is.list(t_14_6_02_ard) && "tbl_summary" %in% names(t_14_6_02_ard)) {
  t_14_6_02_ard$tbl_summary
} else {
  t_14_6_02_ard
}

if (!is.data.frame(ard_tbl) || nrow(ard_tbl) == 0) {
  stop(
    "ARD output is empty for t_14_6_02. Check ADLBC/ADLBH analysis records and display_data derivation filters."
  )
}

saveRDS(ard_tbl, file.path(path$table_ard, "t_14_6_02.rds"))

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
    select(-parcat1, -paramn, -section) %>%
    gt(rowname_col = "param") %>%
    tab_header(
      title = md("**Table 14-6.02**"),
      subtitle = md(
        "**Frequency of Normal and Abnormal (Beyond Normal Range)**  \n**Laboratory Values During Treatment**"
      )
    ) %>%
    tab_stubhead(label = "") %>%
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
      table.border.bottom.style = "hidden"
    )

  gt_tbl
}

gt_tbl <- build_gt(display_data)

## Save PDF output ------------------------------------------------------
gt_tbl %>%
  as_docorator(
    display_name = "t_14_6_02",
    display_loc = path$table_output,
    save_object = FALSE,
    header = fancyhead(
      fancyrow(left = "Protocol: CDISCPILOT01", center = NA, right = doc_pagenum()),
      fancyrow(left = "Population: Safety", center = NA, right = NA)
    ),
    footer = fancyfoot(
      fancyrow(left = "[1] Chi-square p-value (treatment arm comparison of abnormal vs normal)."),
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
  "Table 14.6.02 saved: ",
  file.path(path$table_ard, "t_14_6_02.rds"),
  " and ",
  file.path(path$table_output, "t_14_6_02.pdf")
)
