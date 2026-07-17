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
    rename_with(tolower),
  read_dataset_json(file.path(path$adam, "adlbh.json")) %>%
    convert_blanks_to_na_local() %>%
    rename_with(tolower) %>%
    filter(!(paramcd %in% c("ANISO", "MACROCY", "POLYCHR", "POIKILO", "MICROCY")))
) %>%
  filter(!is.na(aval) & saffl == "Y" &
           ((avisitn > 0 & avisitn < 99) | (avisitn == 99 & aentmtfl == "Y") | ablfl == "Y")) %>%
  select(usubjid, trta, trtan, parcat1, paramcd, param, paramn, avisit, avisitn, aval, chg)

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
    avisit = if_else(avisit == "End of Treatment", "End [2]", avisit)
  ) %>%
  # Add factors for strata variables
  mutate(
    trta = fct_reorder(factor(trta), trtan, .na_rm = TRUE),
    param = fct_reorder(factor(param), paramn, .na_rm = TRUE),
    avisit = fct_reorder(factor(avisit), avisitn, .na_rm = TRUE)
  )

# ----------------------------------------------------------------------------
# Get table summary data
# ----------------------------------------------------------------------------

# Add spaces to align parentheses in SD values with other values in the same column
label_style_number_width <- function(width = 6, digits = 2, ...) {
  function(x) style_number(x, digits = digits, ...) %>% str_pad(width = width, side = "left", pad = " ")
}

# Function for calculating individual block trt/param/visit.
get_individual_summary <- function(data, strata) {

  parameter <- str_split(strata, "#", simplify = TRUE)[2]
  print(paste0("Processing parameter: ", parameter))

  # For Creatine Kinase we have to use larger width to fit SD values
  if (parameter == "Creatine Kinase (U/L)") {
    sd_width <- 6
  } else {
    sd_width <- 5
  }

  # Analysis Value N
  stat_n <-
    data %>%
    tbl_summary(
      by = trta,
      include = avisit,
      statistic = ~"{n}",
      label = list(avisit = parameter)
    ) %>%
    modify_header(all_stat_cols() ~ "N") %>%
    modify_spanning_header(all_stat_cols() ~ "{level}") %>%
    remove_footnote_header()

  # Analysis Value Mean (SD)
  stat_aval <-
    data %>%
    tbl_continuous(
      by = trta,
      variable = aval,
      include = avisit,
      statistic = ~"{mean} ({sd})",
      label = list(avisit = parameter),
      digits = ~list(mean = label_style_number(digits = 1), sd = label_style_number_width(sd_width))
    ) %>%
    modify_header(all_stat_cols() ~ "Mean (SD)") %>%
    modify_spanning_header(all_stat_cols() ~ "{level}") %>%
    remove_footnote_header()

  # Change from baseline Mean (SD)
  stat_chg <-
    data %>%
    tbl_continuous(
      by = trta,
      variable = chg,
      include = avisit,
      statistic = ~"{mean} ({sd})",
      label = list(avisit = parameter),
      digits = ~list(mean = label_style_number(digits = 1), sd = label_style_number_width(sd_width))
    ) %>%
    modify_table_body(
      ~ .x %>%
        mutate(
          across(all_stat_cols(), ~ifelse(label == "Baseline", NA, .))
        )
    ) %>%
    modify_header(all_stat_cols() ~ "Mean (SD) [1]") %>%
    modify_spanning_header(all_stat_cols() ~ "{level}") %>%
    remove_footnote_header()

  tbl_final <-
    list(stat_n, stat_aval, stat_chg) %>%
    tbl_merge(tab_spanner = FALSE) %>%
    modify_header(label = "") %>%
    modify_column_alignment(columns = all_stat_cols(), align = "right") %>%
    # reorder the column to group the by treatment
    modify_table_body(
      ~ .x %>%
        mutate(across(all_stat_cols(), as.character)) %>%
        {
          stat_col_order <- select(., all_stat_cols()) %>% names() %>% sort()
          relocate(., all_of(stat_col_order), .after = "label")
        }
    )
}

t_14_6_1 <- adlb_updated %>%
  tbl_strata2(
    strata = c(category, param),
    .tbl_fun = get_individual_summary,
    .combine_with = "tbl_stack",
    .sep = "#"
  ) %>%
  modify_table_body(
    ~ .x %>%
      select(-groupname_col)
  ) %>%
  modify_header(label = "") %>%
  modify_table_body(
    ~ .x %>%
      mutate(
        category = str_split(tbl_id1_lbl, "#", simplify = TRUE)[, 1],
        param = str_split(tbl_id1_lbl, "#", simplify = TRUE)[, 2],
      ) %>%
      arrange(category, param)
  ) %>%
  modify_indent(columns = label, rows = row_type == "level", indent = 2L) %>%
  # Add spanning rows for each parameter category
  modify_table_body(
    ~ .x %>%
      # Parameter category spanning row
      group_by(category, .add = TRUE) %>%
      reframe(
        add_row(
          pick(everything()),
          row_type = "header",
          label = NA,
          .before = 1
        )
      ) %>%
      group_by(category, .add = TRUE) %>%
      reframe(
        add_row(
          pick(everything()),
          row_type = "header",
          label = first(category),
          .before = 1
        )
      )
  )


# ----------------------------------------------------------------------------
# Output data
# ----------------------------------------------------------------------------

# Convert gtsummary object to ARD format
t_14_6_1_ard <- gather_ard(t_14_6_1)

# Save ARD as RDS
saveRDS(t_14_6_1_ard, file.path(path$table_ard, "t_14_6_1.rds"))

# Save as PDF

# We have to split table into several tables to avoid mid-page split of categories
total_rows <- nrow(t_14_6_1$table_body)
# Find the location of the first Hematology category row as it needs to start on a new page
hematology_row <- which(t_14_6_1$table_body$label == "Hematology")[1]

t_14_6_1_by_page <-
  t_14_6_1 %>%
  tbl_split_by_rows(row_numbers = c(
    seq(14, hematology_row, by = 12),
    seq(hematology_row + 13, total_rows, by = 12)
  ))

gt_tables_list <- map(t_14_6_1_by_page, ~ {
  # For creatine kinase, we need to use different label width to fit SD values
  is_creatine_kinase <- any(str_detect(.x$table_body$label, "Creatine Kinase \\(U/L\\)"))

  if (is.na(is_creatine_kinase) || !is_creatine_kinase) {
    label_width <- "20%"
  } else {
    label_width <- "15%"
  }
  col_width_formula <- as.formula(paste0("label ~ '", label_width, "'"))

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
    cols_width(col_width_formula) %>%
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
        left = "[1] Change from Baseline"
      ),
      fancyrow(
        left = "[2] Last observed value while on treatment (prior to or at Week 24)"
      ),
      fancyrow(left = paste0("Source: ", doc_relative_path()), center = NA, right = doc_datetime())
    ),
    save_object = FALSE
  ) %>%
  render_pdf(path$table_output)
