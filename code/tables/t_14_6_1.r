#************************************************************************
# Purpose:     Generate Table 14.6.1 - Shifts of Hy's Law Values During Treatment
# Input:       ADLBC (+ ADSL for treatment N in headers)
# Output:      t_14_6_1.pdf and t_14_6_1.rds
#************************************************************************

# Note to Reviewer
# To rerun the code below, please refer to the ADRG appendix.
# After required packages are installed, the path variable needs to be defined
# in the .Rprofile file.

# Libraries ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(datasetjson)
library(gt)
library(docorator)

# Helpers ----------------------------------------------------------------
ULN_THRESHOLD <- 1.5
TREATMENT_VISIT_MIN <- 3
TREATMENT_VISIT_MAX <- 12
TRT_LEVELS <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
COL_SEPARATOR <- "__"

convert_blanks_to_na_local <- function(df) {
  df %>%
    mutate(across(where(is.character), ~ na_if(.x, "")))
}

derive_hys_shift <- function(adlbc) {
  adlbc_clean <- adlbc %>%
    convert_blanks_to_na_local() %>%
    mutate(AVISITN = as.numeric(AVISITN))

  adlbc_base <- adlbc_clean %>%
    filter(PARAMCD %in% c("ALT", "AST", "BILI"), AVISIT == "Baseline")

  adlbc_trt <- adlbc_clean %>%
    filter(
      PARAMCD %in% c("ALT", "AST", "BILI"),
      AVISITN >= TREATMENT_VISIT_MIN,
      AVISITN <= TREATMENT_VISIT_MAX
    )

  base_flags <- adlbc_base %>%
    group_by(USUBJID, TRTA, SAFFL) %>%
    summarise(
      trans_base_high = as.integer(any(PARAMCD %in% c("ALT", "AST") & R2A1HI > ULN_THRESHOLD, na.rm = TRUE)),
      trans_base_has_result = any(PARAMCD %in% c("ALT", "AST") & !is.na(R2A1HI)),
      bili_base_high = as.integer(any(PARAMCD == "BILI" & R2A1HI > ULN_THRESHOLD, na.rm = TRUE)),
      bili_base_has_result = any(PARAMCD == "BILI" & !is.na(R2A1HI)),
      .groups = "drop"
    ) %>%
    mutate(
      hylaw_base_high = as.integer(trans_base_high == 1 & bili_base_high == 1),
      hylaw_base_has_result = trans_base_has_result & bili_base_has_result
    )

  visit_flags <- adlbc_trt %>%
    group_by(USUBJID, TRTA, AVISITN) %>%
    summarise(
      trans_high = any(PARAMCD %in% c("ALT", "AST") & R2A1HI > ULN_THRESHOLD, na.rm = TRUE),
      trans_has_result = any(PARAMCD %in% c("ALT", "AST") & !is.na(R2A1HI)),
      bili_high = any(PARAMCD == "BILI" & R2A1HI > ULN_THRESHOLD, na.rm = TRUE),
      bili_has_result = any(PARAMCD == "BILI" & !is.na(R2A1HI)),
      .groups = "drop"
    ) %>%
    mutate(hylaw_high = trans_high & bili_high)

  trt_flags <- visit_flags %>%
    group_by(USUBJID, TRTA) %>%
    summarise(
      trans_any_high = any(trans_high, na.rm = TRUE),
      trans_any_normal = any(trans_has_result & !trans_high, na.rm = TRUE),
      hylaw_any_high = any(hylaw_high, na.rm = TRUE),
      hylaw_any_normal = any((trans_has_result & bili_has_result) & !hylaw_high, na.rm = TRUE),
      .groups = "drop"
    )

  status_from_base <- function(base_high, any_high, any_normal) {
    case_when(
      is.na(base_high) ~ NA_character_,
      base_high == 0 & any_high ~ "High",
      base_high == 0 & !any_high ~ "Normal",
      base_high == 1 & any_normal ~ "Normal",
      base_high == 1 & !any_normal ~ "High",
      TRUE ~ NA_character_
    )
  }

  subject_flags <- base_flags %>%
    left_join(trt_flags, by = c("USUBJID", "TRTA")) %>%
    filter(SAFFL == "Y") %>%
    mutate(
      trans_treatment_status = status_from_base(trans_base_high, trans_any_high, trans_any_normal),
      hylaw_treatment_status = status_from_base(hylaw_base_high, hylaw_any_high, hylaw_any_normal)
    )

  bind_rows(
    subject_flags %>%
      filter(trans_base_has_result) %>%
      transmute(
        PARAM = "Transaminase 1.5 x ULN",
        TRTA,
        baseline_status = if_else(trans_base_high == 1, "High", "Normal"),
        treatment_status = trans_treatment_status
      ),
    subject_flags %>%
      filter(hylaw_base_has_result) %>%
      transmute(
        PARAM = "Total Bili 1.5 x ULN and Transaminase 1.5 x ULN",
        TRTA,
        baseline_status = if_else(hylaw_base_high == 1, "High", "Normal"),
        treatment_status = hylaw_treatment_status
      )
  )
}

calc_cmh_pvalue <- function(param_df) {
  cmh_input <- param_df %>%
    filter(
      !is.na(TRTA),
      !is.na(baseline_status),
      !is.na(treatment_status)
    ) %>%
    mutate(
      TRTA = factor(TRTA, levels = TRT_LEVELS),
      baseline_status = factor(baseline_status, levels = c("Normal", "High")),
      treatment_status = factor(treatment_status, levels = c("Normal", "High"))
    )

  if (nrow(cmh_input) == 0) {
    return(NA_real_)
  }

  cmh_array <- xtabs(~ TRTA + treatment_status + baseline_status, data = cmh_input)
  strata_n <- apply(cmh_array, 3, sum)

  if (any(strata_n <= 1)) {
    return(NA_real_)
  }

  tryCatch(
    as.numeric(stats::mantelhaen.test(cmh_array)$p.value),
    error = function(e) NA_real_
  )
}

format_cell <- function(numerator, denominator) {
  if (is.na(denominator) || denominator == 0) {
    return("0")
  }
  pct <- round(100 * numerator / denominator)
  sprintf("%d (%d%%)", numerator, pct)
}

build_table_rows <- function(param_df, param_label) {
  base_levels <- c("Normal", "High")
  shift_levels <- c("Normal", "High")

  n_df <- param_df %>%
    count(TRTA, baseline_status, name = "n") %>%
    complete(TRTA = TRT_LEVELS, baseline_status = base_levels, fill = list(n = 0))

  shift_df <- param_df %>%
    count(TRTA, baseline_status, treatment_status, name = "n_shift") %>%
    complete(
      TRTA = TRT_LEVELS,
      baseline_status = base_levels,
      treatment_status = shift_levels,
      fill = list(n_shift = 0)
    ) %>%
    left_join(n_df, by = c("TRTA", "baseline_status")) %>%
    mutate(cell = mapply(format_cell, n_shift, n))

  to_wide <- function(df, value_col) {
    df %>%
      mutate(col = paste(TRTA, baseline_status, sep = COL_SEPARATOR)) %>%
      transmute(col, value = as.character(.data[[value_col]])) %>%
      pivot_wider(names_from = col, values_from = value)
  }

  n_row <- to_wide(n_df, "n") %>%
    mutate(PARAM = param_label, Shift = "n", p_value = sprintf("%.3f", calc_cmh_pvalue(param_df)))

  normal_row <- shift_df %>%
    filter(treatment_status == "Normal") %>%
    to_wide("cell") %>%
    mutate(PARAM = "", Shift = "Normal", p_value = "")

  high_row <- shift_df %>%
    filter(treatment_status == "High") %>%
    to_wide("cell") %>%
    mutate(PARAM = "", Shift = "High", p_value = "")

  bind_rows(n_row, normal_row, high_row) %>%
    select(
      PARAM,
      Shift,
      Placebo__Normal,
      Placebo__High,
      `Xanomeline Low Dose__Normal`,
      `Xanomeline Low Dose__High`,
      `Xanomeline High Dose__Normal`,
      `Xanomeline High Dose__High`,
      p_value
    )
}

empty_table_row <- function() {
  tibble::tibble(
    PARAM = "",
    Shift = "",
    Placebo__Normal = "",
    Placebo__High = "",
    `Xanomeline Low Dose__Normal` = "",
    `Xanomeline Low Dose__High` = "",
    `Xanomeline High Dose__Normal` = "",
    `Xanomeline High Dose__High` = "",
    p_value = ""
  )
}

# Load Data ---------------------------------------------------------------
adsl <- read_dataset_json(file.path(path$adam, "adsl.json")) %>%
  convert_blanks_to_na_local()

adlbc <- read_dataset_json(file.path(path$adam, "adlbc.json"), decimals_as_floats = TRUE) %>%
  convert_blanks_to_na_local()

# Build Table -------------------------------------------------------------
shift_subject <- derive_hys_shift(adlbc)

t_14_6_1 <- bind_rows(
  build_table_rows(
    shift_subject %>% filter(PARAM == "Transaminase 1.5 x ULN"),
    "Transaminase 1.5 x ULN"
  ),
  empty_table_row(),
  build_table_rows(
    shift_subject %>% filter(PARAM == "Total Bili 1.5 x ULN and Transaminase 1.5 x ULN"),
    "Total Bili 1.5 x ULN and Transaminase 1.5 x ULN"
  )
)

trt_n <- adsl %>%
  filter(SAFFL == "Y") %>%
  count(TRT01A, name = "N")

get_trt_n <- function(trt_label) {
  n_val <- trt_n %>%
    filter(TRT01A == trt_label) %>%
    pull(N)

  if (length(n_val) == 0) {
    return(0L)
  }

  as.integer(n_val[[1]])
}

table_output <- if (!is.null(path$table_output)) path$table_output else file.path(getwd(), "data/tables/pdf")
table_ard <- if (!is.null(path$table_ard)) path$table_ard else file.path(getwd(), "data/tables/ard")

dir.create(table_output, recursive = TRUE, showWarnings = FALSE)
dir.create(table_ard, recursive = TRUE, showWarnings = FALSE)

saveRDS(t_14_6_1, file.path(table_ard, "t_14_6_1.rds"))

gt_table <- t_14_6_1 %>%
  gt() %>%
  cols_label(
    PARAM = "",
    Shift = "",
    Placebo__Normal = "Normal",
    Placebo__High = "High",
    `Xanomeline Low Dose__Normal` = "Normal",
    `Xanomeline Low Dose__High` = "High",
    `Xanomeline High Dose__Normal` = "Normal",
    `Xanomeline High Dose__High` = "High",
    p_value = "P-value"
  ) %>%
  tab_spanner(
    label = md(paste0("Placebo<br>(N=", get_trt_n("Placebo"), ")")),
    columns = c(Placebo__Normal, Placebo__High)
  ) %>%
  tab_spanner(
    label = md(paste0("Xanomeline Low Dose<br>(N=", get_trt_n("Xanomeline Low Dose"), ")")),
    columns = c(`Xanomeline Low Dose__Normal`, `Xanomeline Low Dose__High`)
  ) %>%
  tab_spanner(
    label = md(paste0("Xanomeline High Dose<br>(N=", get_trt_n("Xanomeline High Dose"), ")")),
    columns = c(`Xanomeline High Dose__Normal`, `Xanomeline High Dose__High`)
  ) %>%
  tab_header(
    title = "Table 14-6.01",
    subtitle = "Shifts of Hy's Law Values During Treatment"
  ) %>%
  cols_align(
    align = "center",
    columns = c(
      Placebo__Normal,
      Placebo__High,
      `Xanomeline Low Dose__Normal`,
      `Xanomeline Low Dose__High`,
      `Xanomeline High Dose__Normal`,
      `Xanomeline High Dose__High`,
      p_value
    )
  ) %>%
  cols_align(align = "left", columns = c(PARAM, Shift))

gt_table %>%
  as_docorator(
    display_name = "t_14_6_1",
    display_loc = table_output,
    header = fancyhead(
      fancyrow(left = "Protocol: CDISCPILOT01", center = NA, right = doc_pagenum()),
      fancyrow(left = "Population: Safety Population", center = NA, right = NA)
    ),
    footer = fancyfoot(
      fancyrow(left = "CMH p-value from general association test stratified by baseline status."),
      fancyrow(left = doc_path(), center = NA, right = doc_datetime())
    )
  ) %>%
  render_pdf(table_output)

print(trt_n)
print(t_14_6_1)
