#************************************************************************
# Purpose:     Generate Table 14.1.3 - Summary of Number of Subjects By Site
# Input:       ADSL
# Output:      t_14_1_3.pdf
#************************************************************************

# Note to Reviewer
# To rerun the code below, please refer ADRG appendix.
# After required package are installed, the path variable needs to be defined
# in the .Rprofile file

# Load required packages
library(dplyr)
library(tidyr)
library(datasetjson)
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
  file.path(path$adam, "adsl.json")
) %>%
  convert_blanks_to_na_local() %>%
  rename_with(tolower)

# ----------------------------------------------------------------------------
# Prepare data
# ----------------------------------------------------------------------------

trt_levels <- c(
  "Placebo",
  "Xanomeline Low Dose",
  "Xanomeline High Dose",
  "Total"
)

adsl_total <- bind_rows(
  adsl,
  adsl %>% mutate(trt01p = "Total")
) %>%
  mutate(
    trt01p = factor(as.character(trt01p), levels = trt_levels),
    itt = as.integer(ittfl == "Y"),
    eff = as.integer(efffl == "Y"),
    com = as.integer(comp24fl == "Y"),
    pooled_id = as.character(sitegr1),
    site_id = as.character(siteid)
  )

site_counts <- adsl_total %>%
  filter(!is.na(pooled_id), !is.na(site_id)) %>%
  group_by(pooled_id, site_id, trt01p) %>%
  summarise(
    itt = sum(itt, na.rm = TRUE),
    eff = sum(eff, na.rm = TRUE),
    com = sum(com, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(
    nesting(pooled_id, site_id),
    trt01p,
    fill = list(itt = 0, eff = 0, com = 0)
  )

totals_by_trt <- adsl_total %>%
  group_by(trt01p) %>%
  summarise(
    itt = sum(itt, na.rm = TRUE),
    eff = sum(eff, na.rm = TRUE),
    com = sum(com, na.rm = TRUE),
    .groups = "drop"
  )

site_counts_wide <- site_counts %>%
  arrange(as.numeric(pooled_id), as.numeric(site_id), trt01p) %>%
  pivot_wider(
    names_from = trt01p,
    values_from = c(itt, eff, com),
    names_glue = "{trt01p}_{.value}"
  ) %>%
  mutate(
    pooled_id = if_else(duplicated(pooled_id), "", pooled_id)
  ) %>%
  select(
    pooled_id,
    site_id,
    `Placebo_itt`, `Placebo_eff`, `Placebo_com`,
    `Xanomeline Low Dose_itt`, `Xanomeline Low Dose_eff`, `Xanomeline Low Dose_com`,
    `Xanomeline High Dose_itt`, `Xanomeline High Dose_eff`, `Xanomeline High Dose_com`,
    `Total_itt`, `Total_eff`, `Total_com`
  )

total_row <- totals_by_trt %>%
  pivot_wider(names_from = trt01p, values_from = c(itt, eff, com), names_glue = "{trt01p}_{.value}") %>%
  mutate(
    pooled_id = "TOTAL",
    site_id = ""
  ) %>%
  select(colnames(site_counts_wide))

t_14_1_3_ard <- bind_rows(site_counts_wide, total_row)

# ----------------------------------------------------------------------------
# Output data
# ----------------------------------------------------------------------------

dir.create(path$table_ard, recursive = TRUE, showWarnings = FALSE)
dir.create(path$table_output, recursive = TRUE, showWarnings = FALSE)

saveRDS(t_14_1_3_ard, file.path(path$table_ard, "t_14_1_3.rds"))

n_placebo <- totals_by_trt %>% filter(trt01p == "Placebo") %>% pull(itt)
n_low <- totals_by_trt %>% filter(trt01p == "Xanomeline Low Dose") %>% pull(itt)
n_high <- totals_by_trt %>% filter(trt01p == "Xanomeline High Dose") %>% pull(itt)
n_total <- totals_by_trt %>% filter(trt01p == "Total") %>% pull(itt)

gt_table <- t_14_1_3_ard %>%
  gt() %>%
  cols_label(
    pooled_id = "Pooled Id",
    site_id = "Site Id",
    `Placebo_itt` = "ITT", `Placebo_eff` = "Eff", `Placebo_com` = "Com",
    `Xanomeline Low Dose_itt` = "ITT", `Xanomeline Low Dose_eff` = "Eff", `Xanomeline Low Dose_com` = "Com",
    `Xanomeline High Dose_itt` = "ITT", `Xanomeline High Dose_eff` = "Eff", `Xanomeline High Dose_com` = "Com",
    `Total_itt` = "ITT", `Total_eff` = "Eff", `Total_com` = "Com"
  ) %>%
  tab_spanner(
    label = md(paste0("Placebo<br>(N=", n_placebo, ")")),
    columns = c(`Placebo_itt`, `Placebo_eff`, `Placebo_com`)
  ) %>%
  tab_spanner(
    label = md(paste0("Xanomeline Low Dose<br>(N=", n_low, ")")),
    columns = c(`Xanomeline Low Dose_itt`, `Xanomeline Low Dose_eff`, `Xanomeline Low Dose_com`)
  ) %>%
  tab_spanner(
    label = md(paste0("Xanomeline High Dose<br>(N=", n_high, ")")),
    columns = c(`Xanomeline High Dose_itt`, `Xanomeline High Dose_eff`, `Xanomeline High Dose_com`)
  ) %>%
  tab_spanner(
    label = md(paste0("Total<br>(N=", n_total, ")")),
    columns = c(`Total_itt`, `Total_eff`, `Total_com`)
  ) %>%
  tab_header(
    title = "Table 14-1.03",
    subtitle = "Summary of Number of Subjects By Site"
  ) %>%
  cols_align(align = "center", columns = everything())

gt_table %>%
  as_docorator(
    display_name = "t_14_1_3",
    display_loc = path$table_output,
    header = fancyhead(
      fancyrow(left = "Protocol: CDISCPILOT01", center = NA, right = doc_pagenum()),
      fancyrow(left = "Population: All Subjects", center = NA, right = NA)
    ),
    footer = fancyfoot(
      fancyrow(
        left = paste0(
          "Note: ITT: Number of subjects in the ITT population, Eff: Number of subjects in the Efficacy population,"
        )
      ),
      fancyrow(
        left = "Com: Number of subjects completing Week 24."
      ),
      fancyrow(left = doc_path(), center = NA, right = doc_datetime())
    )
  ) %>%
  render_pdf(path$table_output)
