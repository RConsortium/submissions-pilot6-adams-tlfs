#************************************************************************
# Purpose:     Generate Table 14.4.1 - Summary of Planned Exposure to Study Drug
# Input:       ADSL
# Output:      t_14_4_01.rds, t_14_4_01.pdf
#************************************************************************

# Run from the project root with paths defined by .Rprofile.

# Load required packages
library(datasetjson)
library(dplyr)
library(gtsummary)
library(gt)
library(docorator)

# Import utility functions
source(file.path(path$utils, "doc_relative_path.r"))

# ----------------------------------------------------------------------------
# Load datasets
# ----------------------------------------------------------------------------

adsl <- read_dataset_json(
  file.path(path$adam, "adsl.json")
) %>%
  mutate(across(where(is.character), ~ na_if(.x, ""))) %>%
  rename_with(tolower)

# ----------------------------------------------------------------------------
# Validate inputs
# ----------------------------------------------------------------------------

required_vars <- c(
  "usubjid", "saffl", "comp24fl",
  "trt01a", "trt01an", "avgdd", "cumdose"
)

stopifnot(all(required_vars %in% names(adsl)))

stopifnot(
  !anyNA(adsl[required_vars]),
  all(trimws(adsl$usubjid) != ""),
  nrow(adsl) == n_distinct(adsl$usubjid),
  all(adsl$saffl %in% c("Y", "N")),
  all(adsl$comp24fl %in% c("Y", "N")),
  is.numeric(adsl$avgdd),
  is.numeric(adsl$cumdose)
)

# ----------------------------------------------------------------------------
# Prepare data
# ----------------------------------------------------------------------------

# Define treatment order
trt_levels <- c(
  "Placebo",
  "Xanomeline Low Dose",
  "Xanomeline High Dose"
)

stopifnot(
  all(adsl$trt01a[adsl$saffl == "Y"] %in% trt_levels)
)

# Prepare the safety population
adsl_safety <- adsl %>%
  filter(saffl == "Y") %>%
  mutate(
    trt01a = factor(trt01a, levels = trt_levels)
  )

# Check population counts
population_counts <- adsl_safety %>%
  group_by(trt01an, trt01a) %>%
  summarise(
    safety_n = n(),
    completers_week24_n = sum(comp24fl == "Y"),
    .groups = "drop"
  ) %>%
  arrange(trt01an)

population_counts

# Select Week 24 completers
adsl_completers <- adsl_safety %>%
  filter(comp24fl == "Y")

# ----------------------------------------------------------------------------
# Get table summary data
# ----------------------------------------------------------------------------

# Calculate exposure summaries
safety_summary <- adsl_safety %>%
  tbl_summary(
    by = trt01a,
    include = c(avgdd, cumdose),
    type = c(avgdd, cumdose) ~ "continuous2",
    statistic = all_continuous() ~ c(
      "{N_nonmiss}",
      "{mean}",
      "{sd}",
      "{median}",
      "{min}",
      "{max}"
    ),
    digits = all_continuous() ~ c(0, 1, 2, 1, 1, 1),
    label = list(
      avgdd = "Average daily dose (mg)",
      cumdose = "Cumulative dose at end of study [2]"
    ),
    missing = "no"
  ) %>%
  modify_header(label = "") %>%
  modify_header(all_stat_cols() ~ "**{level}**  \n(N={n})") %>%
  remove_footnote_header(columns = all_stat_cols())

safety_summary

# Calculate exposure summaries for completers
completer_summary <- adsl_completers %>%
  tbl_summary(
    by = trt01a,
    include = c(avgdd, cumdose),
    type = c(avgdd, cumdose) ~ "continuous2",
    statistic = all_continuous() ~ c(
      "{N_nonmiss}",
      "{mean}",
      "{sd}",
      "{median}",
      "{min}",
      "{max}"
    ),
    digits = all_continuous() ~ c(0, 1, 2, 1, 1, 1),
    label = list(
      avgdd = "Average daily dose (mg)",
      cumdose = "Cumulative dose at end of study [2]"
    ),
    missing = "no"
  ) %>%
  modify_header(label = "") %>%
  modify_header(all_stat_cols() ~ "**{level}**  \n(N={n})") %>%
  remove_footnote_header(columns = all_stat_cols())

completer_summary

# Combine population panels
t_14_4_01 <- tbl_merge(
  tbls = list(completer_summary, safety_summary),
  tab_spanner = c(
    "**Completers at Week 24**",
    "**Safety Population [1]**"
  )
) %>%
  modify_table_body(
    ~ .x %>%
      mutate(
        across(
          all_stat_cols(),
          ~ case_when(
            row_type == "level" & label == "N Non-missing" ~
              paste0(.x, "\u2007\u2007\u2007"),
            row_type == "level" & label != "SD" ~
              paste0(.x, "\u2007"),
            TRUE ~ .x
          )
        )
      )
  )

t_14_4_01

# ----------------------------------------------------------------------------
# Output data
# ----------------------------------------------------------------------------

# Save analysis results
t_14_4_01_ard <- gather_ard(t_14_4_01)

saveRDS(
  t_14_4_01_ard,
  file.path(path$table_ard, "t_14_4_01.rds")
)

# Format the table
gt_table <- t_14_4_01 %>%
  as_gt(auto_align = FALSE) %>%
  cols_align(align = "right", columns = all_stat_cols()) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_labels(columns = everything())
  ) %>%
  tab_header(
    title = "Table 14-4.01",
    subtitle = "Summary of Planned Exposure to Study Drug, as of End of Study"
  )

gt_table

# Render the review PDF
gt_table %>%
  as_docorator(
    display_name = "t_14_4_01",
    display_loc = path$table_output,
    header = fancyhead(
      fancyrow(left = "Protocol: CDISCPILOT01", center = NA, right = doc_pagenum()),
      fancyrow(left = "Population: Safety", center = NA, right = NA)
    ),
    footer = fancyfoot(
      fancyrow(left = "[1] Includes completers and early terminators."),
      fancyrow(left = "[2] End of Study refers to week 26/Early Termination."),
      fancyrow(
        left = paste0("Source: ", doc_relative_path()),
        center = NA,
        right = doc_datetime()
      )
    ),
    save_object = FALSE
  ) %>%
  render_pdf(path$table_output)
