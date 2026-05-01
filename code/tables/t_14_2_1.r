#************************************************************************
# Purpose:     Generate Table 14.2.1 - Summary of Demographic and Baseline Characteristics
# Input:       ADSL
# Output:      t14_2_1.pdf
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

# Import utility functions
source(file.path("code", "utils", "doc_relative_path.r"))

# ----------------------------------------------------------------------------
# Load datasets
# ----------------------------------------------------------------------------

# Local helper to replace blank strings with NA
convert_blanks_to_na_local <- function(df) {
  df %>%
    mutate(across(where(is.character), ~ na_if(.x, "")))
}

adsl <- read_dataset_json(file.path(path$adam, "adsl.json")) %>%
  convert_blanks_to_na_local() %>%
  rename_with(tolower) %>%
  filter(ittfl == "Y")

# ----------------------------------------------------------------------------
# Prepare data
# ----------------------------------------------------------------------------


# BMI Category does not have a numeric vairable in ADSL, so order is defined manually
adsl_factorized <- adsl %>%
  mutate(
    trt01p = fct_reorder(factor(trt01p), trt01pn),
    race = fct_reorder(factor(race), racen),
    agegr1 = fct_reorder(factor(agegr1), agegr1n),
    bmiblgr1 = factor(bmiblgr1, levels = c("<25", "25-<30", ">=30"))
  )

# ----------------------------------------------------------------------------
# Get table summary data
# ----------------------------------------------------------------------------
summary_cols <- c(
  "age", "agegr1", "sex", "race", "mmsetot",
  "durdis", "educlvl", "weightbl", "heightbl", "bmibl", "bmiblgr1"
)

# Create summary table for demographic and baseline characteristics
t_14_2_1 <- adsl_factorized %>%
  tbl_summary(
    by = trt01p,
    type = all_continuous() ~ "continuous2",
    include = all_of(summary_cols),
    statistic = list(
      all_continuous() ~ c("{N_nonmiss}", "{mean}", "{sd}", "{median}", "{min}", "{max}"),
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no",
    digits = list(
      all_continuous() ~ c(0, 1, 2, 1, 1, 1),
      all_categorical() ~ list(
        p = label_style_number(digits = 0, width = 3, justify = "right")
      )
    ),
    label = list(
      age = "Age (y)",
      agegr1 = "Age (y)",
      sex = "Sex",
      race = "Race (Origin)",
      mmsetot = "MMSE",
      durdis = "Duration of Disease",
      educlvl = "Years of Education",
      weightbl = "Baseline Weight (kg)",
      heightbl = "Baseline Height (cm)",
      bmibl = "Baseline BMI",
      bmiblgr1 = "Baseline BMI"
    )
  ) %>%
  add_overall(last = TRUE) %>%
  add_p(
    pvalue_fun = label_style_number(digits = 4),
    test = list(
      all_continuous() ~ "oneway.test",
      all_categorical() ~ "chisq.test"
    ),
    test.args = list(
      all_continuous() ~ list(var.equal = TRUE)
    )
  ) %>%
  add_stat_label(location = "row", label = list(
    all_continuous() ~ c("n", "Mean", "SD", "Median", "Min", "Max"),
    all_categorical() ~ ""
  )) %>%
  modify_header(label = "") %>%
  modify_header(list(
    all_stat_cols() ~ "**{level}**  \n(N={n})",
    p.value ~ "**p-value  \n[1]**"
  )) %>%
  remove_footnote_header(columns = everything()) %>%
  modify_table_body(
    ~ .x %>%
      mutate(variable = factor(variable, levels = summary_cols)) %>%
      group_by(variable, .add = TRUE) %>%
      reframe(add_row(pick(everything()), .before = 1))
  )

tb <- t_14_2_1$table_body

# ----------------------------------------------------------------------------
# Output data
# ----------------------------------------------------------------------------

# Convert gtsummary object to ARD format
t_14_2_1_ard <- gather_ard(t_14_2_1)

# Save ARD as RDS
saveRDS(t_14_2_1_ard, file.path(path$table_ard, "t_14_2_1.rds"))

# Save as PDF
gt_table <-
  t_14_2_1 %>%
  as_gt() %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_labels(columns = everything())
  ) %>%
  cols_width(
    label ~ pct(30),
    everything() ~ pct(13)
  ) %>%
  cols_align_decimal(columns = -label) %>%
  tab_header(
    title = "Table 14-2.01",
    subtitle = "Summary of Demographic and Baseline Characteristics"
  )

# Add header and footer using docorator
gt_table %>%
  as_docorator(
    display_name = "t_14_2_1",
    display_loc  = path$table_output,
    header = fancyhead(
      fancyrow(left = "Protocol: CDISCPILOT01", center = NA, right = doc_pagenum()),
      fancyrow(left = "Population: Intent-to-Treat", center = NA, right = NA)
    ),
    footer = fancyfoot(
      fancyrow(
        left = paste0(
          "[1] P-values are results of ANOVA treatment group comparison for continuous",
          " variable and Pearson's chi-square test for"
        )
      ),
      fancyrow(
        "categorical variables."
      ),
      fancyrow(
        left = paste0(
          "NOTE: Duration of disease is computed as months between date of",
          " enrollment and date of onset of the first definite"
        )
      ),
      fancyrow(
        "symptoms of Alzheimer's disease"
      ),
      fancyrow(left = paste0("Source: ", doc_relative_path()), center = NA, right = doc_datetime())
    ),
    save_object =  FALSE
  ) %>%
  render_pdf(path$table_output)
