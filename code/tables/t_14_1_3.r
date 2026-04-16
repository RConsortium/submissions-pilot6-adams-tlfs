#************************************************************************
# Purpose:     Generate Table 14.1.3 - Summary of Number of Subjects By Site
# Input:       ADSL
# Output:      t14_1_3.rtf
#************************************************************************

# Note to Reviewer
# To rerun the code below, please refer ADRG appendix.
# After required package are installed, the path variable needs to be defined
# in the .Rprofile file

# Load required packages
library(dplyr)
library(tidyr)
library(datasetjson)
library(gtsummary)
library(gt)
library(docorator)

# ----------------------------------------------------------------------------
# Load datasets
# ----------------------------------------------------------------------------

# Local helper to replace blank strings with NA
convert_blanks_to_na_local <- function(df) {
  df %>%
    mutate(across(where(is.character), ~ na_if(.x, "")))
}

adsl <- read_dataset_json(
  file.path(path$adam, "adsl.json"),
) %>%
  convert_blanks_to_na_local() %>%
  rename_with(tolower) %>%
  select(usubjid, siteid, sitegr1, trt01p, ittfl, efffl, comp24fl)

adsl <- read_dataset_json(
  file.path(path$adam, "adsltest.json")
)

# ----------------------------------------------------------------------------
# Prepare data
# ----------------------------------------------------------------------------

# Define treatment groups
trt_levels <- c(
  "Placebo",
  "Xanomeline Low Dose",
  "Xanomeline High Dose",
  "Total"
)

# Define population labels
population_def <- tribble(
  ~population_id, ~population_label,
  "ittfl", "Intent-To-Treat (ITT)",
  "efffl", "Efficacy",
  "comp24fl", "Completer Week 24",
)

# Add total treatment group
adsl_total <- bind_rows(
  adsl,
  adsl %>% mutate(trt01p = "Total")
) %>%
  mutate(trt01p = factor(as.character(trt01p), levels = trt_levels))

adsl_site_combined <- adsl_total %>%
  mutate(siteid = paste(sep = "#", sitegr1, siteid))

# ----------------------------------------------------------------------------
# Get table summary data
# ----------------------------------------------------------------------------

# Map population flag variable names to short labels used in column headers
pop_label_map <- c(ittfl = "ITT", efffl = "Eff", comp24fl = "Com")

# Pivot populations to long (one row per subject × population).
adsl_transposed_by_pop <- adsl_site_combined %>%
  pivot_longer(
    cols      = c(ittfl, efffl, comp24fl),
    names_to  = "population",
    values_to = "value"
  ) %>%
  filter(value == "Y") %>%
  # Map population variable names to short labels for column headers
  mutate(population =
      factor(recode(population, !!!pop_label_map),
        levels = unname(pop_label_map)
      )
  ) %>%
  arrange(sitegr1, siteid, population)

# tbl_summary: columns = trt_pop (combined factor), rows = siteid levels.
# statistic = "{n}" gives raw counts per cell (no denominator needed).
t_14_1_3 <- adsl_transposed_by_pop %>%
  tbl_strata(
    strata = trt01p,
    .tbl_fun =
      ~ .x %>%
      tbl_summary(
        by        = population,
        include   = siteid,
        statistic = all_categorical() ~ "{n}",
        digits    = all_categorical() ~ 0
      ) %>%
      modify_header(all_stat_cols() ~ "**{level}**") %>%
      remove_footnote_header(columns = all_stat_cols()),
    .header = "**{strata}**  \n(N = {n})"
  ) %>%
  # Split siteid back into sitegr1 and siteid for labeling
  modify_table_body(
    ~ .x %>%
      mutate(pooled_id = sapply(strsplit(label, "#"), `[`, 1)) %>%
      mutate(label = sapply(strsplit(label, "#"), `[`, 2)) %>%
      select(pooled_id, everything())
  ) %>%
  modify_header(
    label    = "**Site Id**",
    pooled_id  = "**Pooled Id**",
  )

# ----------------------------------------------------------------------------
# Output data
# ----------------------------------------------------------------------------

# Convert gtsummary object to ARD format
t_14_1_3_ard <- gather_ard(t_14_1_3)

# Save ARD as RDS
saveRDS(t_14_1_3_ard, file.path(path$table_ard, "t_14_1_3.rds"))

# Save as PDF
gt_table <-
  t_14_1_3 %>%
  as_gt() %>%
  cols_align(align = "right") %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_labels(columns = everything())
  ) %>%
  tab_header(
    title = "Table 14-1.03",
    subtitle = "Summary of Number of Subjects by Site"
  )

# Add header and footer using docorator
gt_table %>%
  as_docorator(
    display_name = "t_14_1_3",
    display_loc  = path$table_output,
    header = fancyhead(
      fancyrow(left = "Protocol: CDISCPILOT01", center = NA, right = doc_pagenum()),
      fancyrow(left = "Population: All Subjects", center = NA, right = NA)
    ),
    footer = fancyfoot(
      fancyrow(left =
          "NOTE: ITT: Number of subjects in the ITT population, Eff: Number of subjects in the Efficacy population;"
      ),
      fancyrow(left = "Com: Number of subjects completing Week 24."),
      fancyrow(left = doc_path(), center = NA, right = doc_datetime())
    )
  ) %>%
  render_pdf(path$table_output)