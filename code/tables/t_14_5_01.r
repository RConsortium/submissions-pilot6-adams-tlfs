#*************************************************************************
# Purpose:     Generate Table 14.5.01
#              Incidence of Treatment Emergent Adverse Events by Treatment Group
#
# Input:       ADSL, ADAE
# Output:      t14_5_01.pdf
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

fmt_pvalue <- function(x, flag_threshold = 0.15) {
  formatted <- ifelse(is.na(x), NA_character_,
               ifelse(x > 0.99, ">0.99",
                 formatC(x, digits = 3, format = "f")
               ))
  ifelse(!is.na(x) & x < flag_threshold,
    paste0(formatted, "*"),
    paste0(formatted, "\U2007")  # figure space — same width as *
  )
}

# ----------------------------------------------------------------------------
# Load datasets
# ----------------------------------------------------------------------------

adsl <- read_dataset_json(file.path(path$adam, "adsl.json")) %>%
  convert_blanks_to_na_local() %>%
  rename_with(tolower)

adae <- read_dataset_json(file.path(path$adam, "adae.json")) %>%
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

adae <- adae %>%
  filter(saffl == "Y", trtemfl == "Y") %>%
  left_join(adsl)


# ----------------------------------------------------------------------------
# Get table summary data
# ----------------------------------------------------------------------------

tbl1 <- adae %>%
  tbl_hierarchical(
    variables = c(aebodsys, aedecod),
    by = trt01p,
    id = usubjid,
    denominator = adsl,
    digits = everything() ~ list(p = 1),
    overall_row = TRUE,
    label = list(..ard_hierarchical_overall.. = "Any body system")
  ) %>%
  sort_hierarchical(
    sort = list(aebodsys = "alphanumeric", aedecod = "descending"),
    by_level = "Xanomeline High Dose"
  )

tbl2 <- adae %>%
  tbl_hierarchical_count(
    variables = c(aebodsys, aedecod),
    by = trt01p,
    denominator = adsl,
    digits = everything() ~ list(p = 1),
    overall_row = TRUE,
    label = list(..ard_hierarchical_overall.. = "Any body system")
  ) 

tbl <- tbl_merge(
  tbls = list(tbl1, tbl2),
  tab_spanner = FALSE
) %>% 
  modify_table_body(
  ~ .x %>% relocate(stat_1_2, .after = stat_1_1) %>%
            relocate(stat_2_2, .after = stat_2_1)
) %>%
  modify_header(
    label ~ "System Organ Class/\nPreferred Term",
    stat_1_1 ~ "\n\n\nn (%)",
    stat_1_2 ~ "\n\n\n[AEs]",
    stat_2_1 ~ "\n\n\nn (%)",
    stat_2_2 ~ "\n\n\n[AEs]",
    stat_3_1 ~ "\n\n\nn (%)",
    stat_3_2 ~ "\n\n\n[AEs]"
  ) %>%
  modify_spanning_header(
    c(stat_1_1, stat_1_2) ~ "Placebo\n(N={n})",
    c(stat_2_1, stat_2_2) ~ "Xanomeline Low Dose\n(N={n})",
    c(stat_3_1, stat_3_2) ~ "Xanomeline High Dose\n(N={n})"
  ) %>%
  remove_footnote_header(columns = everything()) %>%
  modify_post_fmt_fun(
    fmt_fun = ~ ifelse(. == "0 (0.0%)", "0", .),
    columns = c(stat_1_1, stat_2_1, stat_3_1)
  ) %>%
  modify_post_fmt_fun(
    fmt_fun = ~ ifelse(is.na(.) | . %in% c("0", "<br />", "NA", ""), "", paste0("[", ., "]")),
    columns = c(stat_1_2, stat_2_2, stat_3_2)
  )

# ----------------------------------------------------------------------------
# Custom derivations for columns with p-values from Fisher's exact test
# ----------------------------------------------------------------------------

compute_fisher <- function(adae, adsl, arm1, arm2) {
  sub_sl <- adsl %>% filter(trt01p %in% c(arm1, arm2))
  sub_ae <- adae %>% filter(trt01p %in% c(arm1, arm2))

  # Helper: Fisher p-value for a subject set vs denominator
  pval <- function(subj_with_event) {
    mat <- table(
      sub_sl$trt01p,
      sub_sl$usubjid %in% subj_with_event
    )
    tryCatch(fisher.test(mat)$p.value, error = function(e) NA_real_)
  }

  # Overall (any TEAE)
  overall_p <- pval(unique(sub_ae$usubjid))

  # By body system
  bodsys_p <- sub_ae %>%
    group_by(aebodsys) %>%
    summarise(p.value = pval(unique(usubjid)), .groups = "drop") %>%
    rename(label = aebodsys)

  # By preferred term
  pt_p <- sub_ae %>%
    group_by(aedecod) %>%
    summarise(p.value = pval(unique(usubjid)), .groups = "drop") %>%
    rename(label = aedecod)

  bind_rows(
    tibble(label = "Any body system", p.value = overall_p),
    bodsys_p,
    pt_p
  )
}

pval_low  <- compute_fisher(adae, adsl, "Placebo", "Xanomeline Low Dose")
pval_high <- compute_fisher(adae, adsl, "Placebo", "Xanomeline High Dose")

tbl_final <- tbl %>%
  modify_table_body(
    ~ .x %>%
      left_join(pval_low  %>% rename(p_low  = p.value), by = "label") %>%
      left_join(pval_high %>% rename(p_high = p.value), by = "label")
  ) %>%
  modify_table_body(~ insert_blank_rows_by(.x,
     group_var = group1_level, insert_before_first = TRUE)
  ) %>%
  modify_fmt_fun(
  p_low  ~ fmt_pvalue,
  p_high ~ fmt_pvalue
) %>%
  modify_header(
    p_low  ~ "Placebo\nvs.\nLow Dose",
    p_high ~ "Placebo\nvs.\nHigh Dose"
  ) %>%
  modify_spanning_header(
    c(p_low, p_high) ~ "Fisher's Exact\np-values"
  ) 

tbl <- tbl_final
# ----------------------------------------------------------------------------
# Output data
# ----------------------------------------------------------------------------
# Convert gtsummary object to ARD format
t_14_5_01_ard <- gather_ard(tbl)

pval_ard <- bind_rows(
  pval_low  %>% mutate(comparison = "Placebo vs Low Dose"),
  pval_high %>% mutate(comparison = "Placebo vs High Dose")
) %>%
  transmute(
    group1         = "trt01p",
    group1_level   = as.list(comparison),
    variable       = "p.value",
    variable_level = as.list(label),
    stat_name      = "p.value",
    stat_label     = "Fisher's Exact p-value",
    stat           = as.list(p.value),
    fmt_fun        = replicate(n(), list()),
    warning        = replicate(n(), list()),
    error          = replicate(n(), list())
  ) %>%
  as_card()

t_14_5_01_ard <- c(t_14_5_01_ard, list(pval_ard))
# Save ARD as RDS
saveRDS(t_14_5_01_ard, file.path(path$table_ard, "t_14_5_01.rds"))

# Save as PDF
gt_table <-
  tbl %>%
  as_gt() %>%
 cols_width(
   label ~ pct(20),
   all_stat_cols() ~ pct(80 / 8),
   p_low ~ pct(80 / 8),
   p_high ~ pct(80 / 8)
 ) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = list(
      cells_column_spanners(spanners = everything()),
      cells_column_labels(columns = all_stat_cols())
    )
  ) %>%
  tab_header(
    title = "Table 14-5.01",
    subtitle = "Incidence of Treatment Emergent Adverse Events by Treatment Group"
  )

# Add header and footer using docorator
gt_table %>%
  as_docorator(
    display_name = "t_14_5_01",
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
