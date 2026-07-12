#************************************************************************
# Purpose:     Generate Table 14.6.1 - Summary Statistics for Continuous Laboratory Values
# Input:       ADSL, ADLBH, ADLBC
# Output:      t14_6_1.pdf
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
library(forcats)
library(purrr)
# Import utility functions
source(file.path("code", "utils", "doc_relative_path.r"))
source(file.path("code", "utils", "table_functions.r"))

# ----------------------------------------------------------------------------
# Load datasets
# ----------------------------------------------------------------------------

# Local helper to replace blank strings with NA
convert_blanks_to_na_local <- function(df) {
  df %>%
    mutate(across(where(is.character), ~ na_if(.x, "")))
}

# Convert all variable names to lowercase for code style consistency
adsl <- read_dataset_json(file.path(path$adam, "adsl.json")) %>%
  convert_blanks_to_na_local() %>%
  rename_with(tolower) %>%
  filter(saffl == "Y")

# Keep only records used for analysis;
adlb <- bind_rows(
  read_dataset_json(file.path(path$adam, "adlbc.json")) %>%
    convert_blanks_to_na_local() %>%
    rename_with(tolower) %>%
    filter(saffl == "Y" & ((avisitn > 0 & avisitn < 99) | (avisitn == 99 & aentmtfl == "Y") | ablfl == "Y")),
  read_dataset_json(file.path(path$adam, "adlbh.json")) %>%
    convert_blanks_to_na_local() %>%
    rename_with(tolower) %>%
    filter(!(paramcd %in% c("ANISO", "MACROCY", "POLYCHR", "POIKILO", "MICROCY"))) %>%
    filter(saffl == "Y" & ((avisitn > 0 & avisitn < 99) | (avisitn == 99 & aentmtfl == "Y") | ablfl == "Y"))
) %>%
  select(usubjid, trtp, trtpn, parcat1, paramcd, param, paramn, avisit, avisitn, aval, chg)

# ----------------------------------------------------------------------------
# Prepare data
# ----------------------------------------------------------------------------

adlb_updated <- adlb %>%
  # Update Parameter Category labels
  mutate(
    category = case_when(
      parcat1 == "CHEM" ~ "Chemistry",
      parcat1 == "HEM" ~ "Hematology",
      TRUE ~ parcat1
    ),
  ) %>%
  # Update End of Treatment visit label
  mutate(
    avisit = if_else(avisit == "End of Treatment", "End [1]", avisit)
  ) %>%
  # Add factors for strata variables
  mutate(
    trtp = fct_reorder(factor(trtp), trtpn, .na_rm = TRUE),
    param = fct_reorder(factor(param), paramn, .na_rm = TRUE),
    avisit = fct_reorder(factor(avisit), avisitn, .na_rm = TRUE)
  )

# ----------------------------------------------------------------------------
# Get table summary data
# ----------------------------------------------------------------------------

# Add spaces to align parentheses in SD values with other values in the same column
format_meansd <- function(str, width) {
  if (grepl(".*\\([\\d\\.]+\\)", str, perl = TRUE)) {
    # Extract the string inside the parentheses
    inside_parens <- sub(".*\\(([\\d\\.]+)\\).*", "\\1", str, perl = TRUE)
    # Calculate the number of spaces needed to be added
    num_spaces <- width - nchar(inside_parens)
    # Add spaces before the string inside the parentheses
    if (num_spaces < 0) {
      num_spaces <- 0
    }
    str <- sub("\\(([\\d\\.]+)\\)", paste0("(", strrep("\u2007", num_spaces), "\\1)"), str, perl = TRUE)
  } else {
    str
  }
}

# Function for calculating individual block trt/param/visit.
get_individual_summary <- function(data) {

  # Analysis Value
  aval_summary <-
    data %>%
    mutate(summary_var = aval) %>%
    tbl_wide_summary(
      type = summary_var ~ "continuous",
      include = summary_var,
      statistic = c("{N_nonmiss}", "{mean} ({sd})"),
      digits = summary_var ~ c(0, 1, 2),
    ) %>%
    modify_header(!!!list(
      stat_1 = "N",
      stat_2 = "Mean (SD)"
    )) %>%
    modify_table_body(
      ~ .x %>%
        mutate(
          stat_2 = format_meansd(stat_2, width = 6)
        )
    ) %>%
    remove_footnote_header(columns = everything())

  # Change from baseline
  chg_summary <-
    data %>%
    mutate(summary_var = chg) %>%
    tbl_wide_summary(
      type = summary_var ~ "continuous",
      include = summary_var,
      statistic = c("{mean} ({sd})"),
      digits = summary_var ~ c(1, 2),
    ) %>%
    modify_table_body(
      ~ .x %>%
        mutate(
          stat_1 = format_meansd(stat_1, width = 6)
        )
    ) %>%
    modify_header(stat_1 = "Mean (SD) [2]") %>%
    remove_footnote_header(columns = everything())

  tbl_merge(
    tbls = list(aval_summary, chg_summary),
    merge_vars = c("tbl_id1", "tbl_id1_lbl", "label", "param", "row_type", "variable", "var_type", "var_label")
  )
}

t_14_6_1 <- adlb_updated %>%
  # Level 1 horizontal stratification: Treatment
  tbl_strata2(
    strata = trtp,
    .combine_args = list(
      merge_vars = c("tbl_id1", "tbl_id1_lbl", "label", "variable", "var_type", "var_label", "row_type")
    ),
    .tbl_fun =
      ~ .x %>%
        # Level 2 vertical stratification: Parameter Category, Parameter, and Visit
        tbl_strata2(
          strata = c(category, param, avisit),
          .tbl_fun = ~ get_individual_summary(.x),
          .combine_with = "tbl_stack",
          .sep = "#"
        ) %>%
        modify_table_body(
          ~ .x %>%
            select(-groupname_col)
        )
  ) %>%
  modify_header(label = "") %>%
  modify_table_body(
    ~ .x %>%
      mutate(
        category = str_split(tbl_id1_lbl, "#", simplify = TRUE)[, 1],
        param = str_split(tbl_id1_lbl, "#", simplify = TRUE)[, 2],
        label = str_split(tbl_id1_lbl, "#", simplify = TRUE)[, 3],
      ) %>%
      # Set baseline values to blank
      mutate(
        stat_1_2_1 = if_else(label == "Baseline", "", stat_1_2_1),
        stat_1_2_2 = if_else(label == "Baseline", "", stat_1_2_2),
        stat_1_2_3 = if_else(label == "Baseline", "", stat_1_2_3)
      )
  ) %>%
  modify_indent(columns = label, rows = row_type == "label", indent = 4L)


# ----------------------------------------------------------------------------
# Output data
# ----------------------------------------------------------------------------

# Convert gtsummary object to ARD format
t_14_6_1_ard <- gather_ard(t_14_6_1)

# Save ARD as RDS
saveRDS(t_14_6_1_ard, file.path(path$table_ard, "t_14_6_1.rds"))

# Save as PDF

# Add spanning rows:
# For each parameter category
# For each parameter

t_14_6_1_grouped <- t_14_6_1 %>%
  modify_table_body(
    ~ .x %>%
      group_by(param, .add = TRUE) %>%
      # Parameter spanning row
      reframe(
        add_row(
          pick(everything()),
          label = first(param),
          row_type = "header",
          category = first(category),
          .before = 1
        )
      ) %>%
      # Parameter blank row
      reframe(
        add_row(
          pick(everything()),
          row_type = "header",
          category = first(category),
          .before = 1
        )
      ) %>%
      # Parameter category spanning row
      group_by(category, .add = TRUE) %>%
      reframe(
        add_row(
          pick(everything()),
          label = first(category),
          .before = 1
        )
      )
  )

# We have to split table into several tables to avoid mid-page split of categories
total_rows <- nrow(t_14_6_1_grouped$table_body)
# Find the location of the first Hematology category row as it needs to start on a new page
hematology_row <- which(t_14_6_1_grouped$table_body$label == "Hematology")[1]

t_14_6_1_by_page <-
  t_14_6_1_grouped %>%
  tbl_split_by_rows(row_numbers = c(
    seq(14, hematology_row, by = 12),
    seq(hematology_row + 12, total_rows, by = 12)
  ))

gt_tables_list <- map(t_14_6_1_by_page, ~ {
  .x %>%
    as_gt(auto_align = FALSE) %>%
    tab_style(
      style = cell_text(align = "center", v_align = "bottom", whitespace = "pre"),
      locations = cells_column_labels(columns = everything())
    ) %>%
    cols_align(
      align = "right",
      columns = all_stat_cols()
    ) %>%
    cols_align(
      align = "auto",
      columns = label
    ) %>%
    cols_width(
      label ~ "20%",
    ) %>%
    tab_header(
      title = "Table 14-6.01",
      subtitle = "Summary Statistics for Continuous Laboratory Values"
    )
})

# Add header and footer using docorator
gt_group(.list = gt_tables_list) %>%
  as_docorator(
    display_name = "t_14_6_1",
    display_loc = path$table_output,
    tbl_scale = FALSE,
    header = fancyhead(
      fancyrow(left = "Protocol: CDISCPILOT01", center = NA, right = doc_pagenum()),
      fancyrow(left = "Population: Safety", center = NA, right = NA)
    ),
    footer = fancyfoot(
      fancyrow(
        left = "[1] Last observed value while on treatment (prior to or at Week 24)"
      ),
      fancyrow(
        left = "[2] Change from Baseline"
      ),
      fancyrow(left = paste0("Source: ", doc_relative_path()), center = NA, right = doc_datetime())
    ),
    save_object = FALSE
  ) %>%
  render_pdf(path$table_output)
