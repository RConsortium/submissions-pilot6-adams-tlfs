#************************************************************************
# Purpose:     Generate Table 14.3.02 - Primary Endpoint Analysis: CIBIC+ - Summary at Week 24 - LOCF
# Input:       ADCIBC
# Output:      t14_3_02.pdf
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
library(emmeans)

# Import utility functions
source(file.path("code", "utils", "doc_relative_path.r"))
source(file.path("code", "utils", "table_functions.r"))

# ----------------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------------

# Local helper to replace blank strings with NA
convert_blanks_to_na_local <- function(df) {
  df %>%
    mutate(across(where(is.character), ~ na_if(.x, "")))
}

# Insert a blank row before each new value of `group_var` (must be ordered before calling).
insert_blank_rows_by <- function(df, group_var) {
  group_var <- rlang::ensym(group_var)

  df_marked <- df %>%
    mutate(.insert = !!group_var != lag(!!group_var, default = (!!group_var)[1]) & row_number() > 1)

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

# Convert all variable names to lowercase for code style consistency
adsl <- read_dataset_json(file.path(path$adam, "adsl.json")) %>%
  convert_blanks_to_na_local() %>%
  rename_with(tolower) %>%
  filter(efffl == "Y") %>%
  mutate(
    trt01p = fct_reorder(factor(trt01p), trt01pn),
    trt01p_copy = trt01p
  )

adcibc <- read_dataset_json(file.path(path$adam, "adcibc.json")) %>%
  convert_blanks_to_na_local() %>%
  rename_with(tolower) %>%
  filter(efffl == "Y" & avisit == "Week 24" & anl01fl == "Y") %>%
  mutate(
    sitegr1 = factor(sitegr1),
    trtp = fct_reorder(factor(trtp), trtpn),
    trtp_rev = fct_reorder(factor(trtp), -trtpn),
  ) %>%
  select(trtp_rev, trtpn, trtp, aval, sitegr1, avisit)

# ----------------------------------------------------------------------------
# Get table summary data
# ----------------------------------------------------------------------------

# Population counts for the header
pop_count_ard <- adsl %>%
  rename(trtp = trt01p) %>%
  ard_tabulate(
    variables = trtp,
    statistic = trtp ~ "n"
  )

# Create continuous summary CIBIC+
continuous_stat_ard <- adcibc %>%
  ard_summary(
    by = trtp,
    variables = aval,
    statistic = ~ continuous_summary_fns(c(
      "N",
      "mean",
      "sd",
      "median",
      "min",
      "max"
    )),
    fmt_fun = list(aval ~ list(N ~ 0, mean ~ 1, sd ~ 2, median ~ 1, min ~ 0, max ~ 0))
  )

# Calculate dose-response statistics
dose_resp <- function(data, full_data, variable) {
  # Dose response is reported only in the last column
  if (data$trtp_rev[1] != "Xanomeline High Dose") {
    return(
      list(dose_resp = NA)
    )
  }
  model <- stats::lm(as.formula(paste0(variable, " ~ trtpn + sitegr1")), data = full_data)

  p_val <- stats::coef(summary(model))["trtpn", "Pr(>|t|)"]

  list(dose_resp = p_val)
}

dose_resp_ard <- adcibc %>%
  ard_mvsummary(
    by = trtp,
    variables = aval,
    statistic = everything() ~ list(dose_resp = \(data, full_data, variable, ...) dose_resp(data, full_data, variable)),
    stat_label = everything() ~ list(dose_resp = "p-value (Dose Response) [1][2]"),
    fmt_fun = list(aval ~ list(dose_resp ~ 3))
  )

# Pairwise comparisons of treatment groups using ANCOVA model with treatment and site group as factors.
form <- stats::as.formula(paste0("aval  ~ trtp_rev + sitegr1"))
model <- stats::lm(form, data = adcibc)

# Estimated marginal means and unadjusted pairwise contrasts
lsm <- emmeans::lsmeans(model, specs = "trtp_rev", weights = "proportional")

pairs_res <- pairs(lsm, adjust = "none")
pairs_ci <- confint(pairs_res)

# Xanomeline - Placebo
xan_pla_comp <- function(data, full_data, by, variable) {
  current_trt <- as.character(data$trtp_rev[1])

  # Nothing is reported in the Placebo column
  if (current_trt == "Placebo") {
    return(
      list(
        est = NA,
        se = NA,
        p_val = NA,
        lower = NA,
        upper = NA
      )
    )
  }

  # For Xanomeline Low Dose, we show Xanomeline Low - Placebo
  # For Xanomeline High Dose, we show Xanomeline High - Placebo
  match_row <- as.data.frame(pairs_res) %>%
    filter(contrast == paste0(current_trt, " - Placebo"))

  match_ci <- as.data.frame(pairs_ci) %>%
    filter(contrast == paste0(current_trt, " - Placebo"))

  # Parse the results
  if (nrow(match_row) > 0) {
    est <- match_row$estimate
    se <- match_row$SE
    p <- match_row$p.value
    lower <- match_ci$lower.CL
    upper <- match_ci$upper.CL

    list(
      est = est,
      se = se,
      p_val = p,
      lower = lower,
      upper = upper
    )
  } else {
    list(
      est = NA,
      se = NA,
      p_val = NA,
      lower = NA,
      upper = NA
    )
  }
}

xan_pla_comp_ard <- adcibc %>%
  ard_mvsummary(
    by = trtp,
    variables = aval,
    statistic = everything() ~ list(xan_pla_comp = \(data, full_data, by, variable, ...) {
      xan_pla_comp(data, full_data, by, variable)
    }),
    stat_label = everything() ~ list(
      p_val = "p-value (Xan - Placebo) [1][3]",
      est = "Diff of LS Means",
      se = "SE",
      lower = "95% CI",
      upper = "95% CI"
    ),
    fmt_fun = list(aval ~ list(est ~ 1, se ~ 2, p_val ~ 3, lower ~ 1, upper ~ 1))
  ) %>%
  mutate(
    context = "Xan - Placebo comparison"
  )

# Xanomeline High - Xanomeline Low
xan_xan_comp <- function(data, full_data, by, variable) {
  # The comparison is reported only in the Xanomeline High Dose column
  if (as.character(data$trtp_rev[1]) != "Xanomeline High Dose") {
    return(
      list(
        est_xan = NA,
        se_xan = NA,
        p_val_xan = NA,
        lower_xan = NA,
        upper_xan = NA
      )
    )
  }

  match_row <- as.data.frame(pairs_res) %>%
    filter(contrast == "Xanomeline High Dose - Xanomeline Low Dose")

  match_ci <- as.data.frame(pairs_ci) %>%
    filter(contrast == "Xanomeline High Dose - Xanomeline Low Dose")

  list(
    est_xan = match_row$estimate,
    se_xan = match_row$SE,
    p_val_xan = match_row$p.value,
    lower_xan = match_ci$lower.CL,
    upper_xan = match_ci$upper.CL
  )
}

xan_xan_comp_ard <- adcibc %>%
  ard_mvsummary(
    by = trtp,
    variables = aval,
    statistic = everything() ~ list(xan_xan_comp = \(data, full_data, by, variable, ...) {
      xan_xan_comp(data, full_data, by, variable)
    }),
    stat_label = everything() ~ list(
      p_val_xan = "p-value (Xan High - Xan Low) [1][3]",
      est_xan = "  Diff of LS Means (Xan High - Xan Low)",
      se_xan = "SE",
      lower_xan = "95% CI (Xan High - Xan Low)",
      upper_xan = "95% CI (Xan High - Xan Low)"
    ),
    fmt_fun = list(aval ~ list(est_xan ~ 1, se_xan ~ 2, p_val_xan ~ 3, lower_xan ~ 1, upper_xan ~ 1))
  ) %>%
  mutate(
    context = "Xan - Xan comparison"
  )


# Get final ARD
united_ard <- bind_ard(
  pop_count_ard,
  continuous_stat_ard,
  dose_resp_ard,
  xan_pla_comp_ard,
  xan_xan_comp_ard
)

t_14_3_2 <- united_ard %>%
  tbl_ard_summary(
    by = trtp,
    type = aval ~ "continuous2",
    include = aval,
    statistic = list(
      all_continuous() ~ c(
        "{N}", "{mean} ({sd})", "{median} ({min}; {max})",
        "{dose_resp}", "{p_val}", "{est} ({se})", "({lower}; {upper})",
        "{p_val_xan}", "{est_xan} ({se_xan})", "({lower_xan}; {upper_xan})"
      )
    ),
    label = list(
      aval = "Week 24"
    )
  ) %>%
  modify_header(!!!list(
    all_stat_cols() ~ "**{level}**  \n**(N={n})**",
    label = ""
  )) %>%
  # Various formatting modifications
  modify_table_body(
    ~ .x %>%
      # Group stats
      mutate(stat_type = case_when(
        row_number() %in% 1:4 ~ "cont",
        row_number() == 5 ~ "doseresp",
        row_number() %in% 6:8 ~ "pairwise_xan_pla",
        row_number() %in% 9:11 ~ "pairwise_xan_xan",
      )) %>%
      mutate(label = str_remove_all(label, fixed(" (Xan High - Xan Low)"))) %>%
      # Remove the contrast tag that keeps the Xan-Xan labels unique in the ARD
      mutate(label = str_remove_all(label, fixed(" (Xan High - Xan Low)"))) %>%
      # Manually set label to match shells
      mutate(label = case_when(
        (row_type == "level" & label == "(95% CI; 95% CI)") ~ "  95% CI",
        (row_type == "level" & label == "N") ~ "n",
        (row_type == "level" & label == "Median (Min; Max)") ~ "Median (Range)",
        (row_type == "level" & label == "Diff of LS Means (SE)") ~ paste0("  ", label),
        TRUE ~ label
      )) %>%
      # Show empty cells where a comparison is not applicable (statistic is NA).
      # Kept in the table only; the ARD still stores NA in the RDS.
      mutate(across(
        all_stat_cols(),
        ~ case_when(
          .x %in% c("NA", "NA (NA)", "(NA; NA)") ~ "",
          TRUE ~ .x
        )
      )) %>%
      insert_blank_rows_by(group_var = stat_type)
  )

# ----------------------------------------------------------------------------
# Output data
# ----------------------------------------------------------------------------

# Save the combined ARD as RDS
saveRDS(united_ard, file.path(path$table_ard, "t_14_3_02.rds"))

# Save as PDF
gt_table <-
  t_14_3_2 %>%
  as_gt(auto_align = FALSE) %>%
  cols_align(
    align = "center",
    columns = all_stat_cols()
  ) %>%
  cols_align(
    align = "left",
    columns = label
  ) %>%
  tab_style(
    style = cell_text(align = "center", v_align = "bottom"),
    locations = cells_column_labels(columns = everything())
  ) %>%
  cols_width(
    label ~ "40%",
    all_stat_cols() ~ "20%"
  ) %>%
  tab_header(
    title = "Table 14-3.02",
    subtitle = "Primary Endpoint Analysis: CIBIC+ - Summary at Week 24 - LOCF"
  )

# Add header and footer using docorator
gt_table %>%
  as_docorator(
    display_name = "t_14_3_02",
    display_loc = path$table_output,
    tbl_scale = FALSE,
    header = fancyhead(
      fancyrow(left = "Protocol: CDISCPILOT01", center = NA, right = doc_pagenum()),
      fancyrow(left = "Population: Efficacy", center = NA, right = NA)
    ),
    footer = fancyfoot(
      fancyrow(
        left = "[1] Based on Analysis of covariance (ANCOVA) model with treatment and site group as factors."
      ),
      fancyrow(
        left = "[2] Test for a non-zero coefficient for treatment (dose) as a continuous variable."
      ),
      fancyrow(
        left = paste0(
          "[3] Pairwise comparison with treatment as a categorical variable: p-values",
          " without adjustment for multiple comparisons."
        )
      ),
      fancyrow(left = paste0("Source: ", doc_relative_path()), center = NA, right = doc_datetime())
    ),
    save_object = FALSE
  ) %>%
  render_pdf(path$table_output)
