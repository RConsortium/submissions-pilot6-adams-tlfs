#*************************************************************************
# Purpose:     Generate Table 14.7.04 - Summary of Concomitant Medications
#              (Number of Subjects)
# Input:       ADSL, CM
# Output:      t14_7_4.pdf
#************************************************************************

# Note to Reviewer
# To rerun the code below, please refer ADRG appendix.
# After required package are installed, the path variable needs to be defined
# in the .Rprofile file

# Load required packages
library(dplyr)
library(datasetjson)
library(cards)
library(gtsummary)
library(gt)
library(docorator)
library(stringr)

# Import utility functions
source(file.path("code", "utils", "doc_relative_path.r"))

# ----------------------------------------------------------------------------
# Local script functions
# ----------------------------------------------------------------------------

# Local helper to replace blank strings with NA
convert_blanks_to_na_local <- function(df) {
  df %>%
    mutate(across(where(is.character), ~ na_if(.x, "")))
}

# Insert a blank row before each new value of `group_var` (must be ordered before calling).
insert_blank_rows_by <- function(df, group_var, insert_before_first = FALSE) {
  group_var <- rlang::ensym(group_var)

  df_marked <- df %>%
    mutate(
      .prev = lag(!!group_var),
      .insert = !is.na(!!group_var) & (is.na(.prev) | !!group_var != .prev) & (insert_before_first | row_number() > 1)
    ) %>%
    select(-.prev)

  blank <- df_marked[1, ] %>% mutate(across(everything(), ~NA))
  idx <- which(df_marked$.insert)
  out <- df_marked
  for (i in rev(idx)) {
    out <- bind_rows(out[seq_len(i - 1), ], blank, out[i:nrow(out), ])
  }

  out %>% select(-.insert)
}
# ----------------------------------------------------------------------------
# Load datasets
# ----------------------------------------------------------------------------

adsl <- read_dataset_json(file.path(path$adam, "adsl.json")) %>%
  convert_blanks_to_na_local() %>%
  rename_with(tolower)

cm <- read_dataset_json(file.path(path$sdtm, "cm.json")) %>%
  convert_blanks_to_na_local() %>%
  rename_with(tolower)

# ----------------------------------------------------------------------------
# Prepare data
# ----------------------------------------------------------------------------

# Define treatment groups
trt_levels <- c(
  "Placebo",
  "Xanomeline Low Dose",
  "Xanomeline High Dose"
)

# Define population labels
population <- "saffl"

# Add treatment group levels
adsl <- adsl %>%
  filter(saffl == "Y") %>%
  select(usubjid, saffl, trt01p)

cm <- cm %>%
  left_join(adsl, by = "usubjid")


# ----------------------------------------------------------------------------
# Get table summary data
# ----------------------------------------------------------------------------

tbl <- cm %>%
  tbl_hierarchical(
    variables = c(cmclas, cmdecod),
    by = trt01p,
    id = usubjid,
    denominator = adsl,
    digits = everything() ~ list(p = 0),
    overall_row = TRUE,
    label = list(..ard_hierarchical_overall.. = "Patients receiving at least one concomitant medication")
  ) %>%
  modify_post_fmt_fun(
    fmt_fun = ~ ifelse(. == "0 (0%)", "0", .),
    columns = all_stat_cols()
  )

tbl <- tbl %>%
  sort_hierarchical(
    sort = list(cmclas = "alphanumeric", cmdecod = "descending"),
    # note gtsummary version >=2.6.0 is required for new by_level feature
    # a factor by variable could not be used here, so left character
    by_level = "Placebo"
  ) %>%
  modify_header(label = "\nTherapeutic class, n (%)") %>%
  modify_header(!!!list(
    all_stat_cols() ~ "**{level}**  \n(N={n})"
  )) %>%
  remove_footnote_header(columns = everything()) %>%
  # reorder the columns, second column should be Low Dose, not High Dose
  modify_table_body(
    ~ .x %>% relocate(stat_3, .before = stat_2)
  ) %>%
  modify_table_body(~ insert_blank_rows_by(.x, group_var = group1_level, insert_before_first = TRUE))


# ----------------------------------------------------------------------------
# Output data
# ----------------------------------------------------------------------------

# Convert gtsummary object to ARD format
t_14_7_04_ard <- gather_ard(tbl)

# Save ARD as RDS
saveRDS(t_14_7_04_ard, file.path(path$table_ard, "t_14_7_04_ard.rds"))

# Save as PDF
gt_table <-
  tbl %>%
  as_gt() %>%
  cols_width(
    label ~ pct(60),
    all_stat_cols() ~ pct(40 / 3)
  ) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = list(
      cells_column_spanners(spanners = everything()),
      cells_column_labels(columns = all_stat_cols())
    )
  ) %>%
  tab_header(
    title = "Table 14-7.04",
    subtitle = "Summary of Concomitant Medications (Number of Subjects)"
  )

# Add header and footer using docorator
gt_table %>%
  as_docorator(
    display_name = "t_14_7_04",
    tbl_scale = FALSE,
    display_loc = path$table_output,
    header = fancyhead(
      fancyrow(left = "Protocol: CDISCPILOT01", center = NA, right = doc_pagenum()),
      fancyrow(left = "Population: Safety Subjects", center = NA, right = NA)
    ),
    footer = fancyfoot(
      fancyrow(left = paste0("Source: ", doc_relative_path()), center = NA, right = doc_datetime())
    ),
    save_object = FALSE
  ) %>%
  render_pdf(path$table_output)
