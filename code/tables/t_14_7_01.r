# Load required packages
library(dplyr)
library(rlang)
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

replace_with_mapping <- function(df, var, mapping, new_var = NULL) {
  sym_var <- ensym(var)
  var_name <- as_name(sym_var)

  new_expr <- enexpr(new_var)
  if (!is.null(new_expr)) {
    sym_new <- ensym(new_var)
    new_name <- as_name(sym_new)
    if (new_name %in% names(df)) {
      stop("Argument `new_var` cannot be an existing column name.")
    }
  } else {
    new_name <- NULL
  }

  out_name <- new_name %||% var_name
  df %>%
    mutate("{out_name}" := replace_values(!!sym_var, from = names(mapping), to = mapping))
}


suppress_repeats <- function(x) {
  changed <- x != lag(x)
  changed[1] <- TRUE
  if (is.numeric(x)) {
    ifelse(changed, x, NA)
  } else {
    ifelse(changed, x, "")
  }
}

insert_blank_rows <- function(df, col) {
  change_idx <- which(
    df[[col]] != lead(df[[col]], default = df[[col]][nrow(df)])
  )

  blank <- df[1, ] %>% mutate(across(everything(), ~NA))
  for (i in rev(change_idx)) {
    df <- bind_rows(df[seq_len(i), ], blank, df[(i + 1):nrow(df), ])
  }
  df
}

# ----------------------------------------------------------------------------
# Load data
# ----------------------------------------------------------------------------
adsl <- read_dataset_json(file.path(path$adam, "adsl.json")) %>%
  convert_blanks_to_na_local() %>%
  rename_with(tolower)

advs <- read_dataset_json(file.path(path$adam, "advs.json")) %>%
  convert_blanks_to_na_local() %>%
  rename_with(tolower)

# ----------------------------------------------------------------------------
# Prepare data and data manipulation steps
# ----------------------------------------------------------------------------
tableid <- "t_14_7_01"
# --- 1. Pre-aggregate variables used in script ---
population <- exprs(saffl)
trtvar <- exprs(trta)
trtvar_n <- exprs(trtan)
trtvar_adsl <- exprs(trt01a)
trtvar_adsl_n <- exprs(trt01an)
analvar <- exprs(aval)

groupbyvars <- exprs(param, atpt, trt_disp, avisit)
groupbyvars_n <- exprs(paramn, atptn, trt_dispn, avisitn)

# for usage with replace_with_mapping
trt_relabel <- list(
  "Xanomeline Low Dose" = "Xan.Low",
  "Xanomeline High Dose" = "Xan.High"
)

visit_relabel <- list(
  "End of Treatment" = "End of Trt."
)

# data manipulation step
adsl <- adsl %>%
  filter(!!population[[1]] == "Y") %>%
  replace_with_mapping(!!trtvar_adsl[[1]], trt_relabel, new_var = trt_disp) %>%
  mutate(trt_dispn = !!trtvar_adsl_n[[1]])

# modify treatment group + visit names
advs <- advs %>%
  filter(!!population[[1]] == "Y") %>%
  filter(
    paramcd %in% c("SYSBP", "DIABP", "PULSE"),
    ablfl == "Y" | avisit %in% c("Week 24", "End of Treatment")
  ) %>%
  replace_with_mapping(!!trtvar[[1]], trt_relabel, new_var = trt_disp) %>%
  mutate(trt_dispn = !!trtvar_n[[1]]) %>%
  replace_with_mapping(avisit, visit_relabel) %>%
  select(
    usubjid, !!analvar[[1]], !!!groupbyvars, !!trtvar[[1]], trt_disp,
    !!!groupbyvars_n, !!trtvar_n[[1]]
  )

# ----------------------------------------------------------------------------
# Perform statistical analysis
# ----------------------------------------------------------------------------
# get N at treatment level
n_trt <- adsl %>%
  group_by(trt_disp, trt_dispn) %>%
  summarise(n_trt = n_distinct(usubjid), .groups = "drop")

# get statistics at param/pos/trt/visit level
vs_summary <- advs %>%
  group_by(!!!groupbyvars, !!!groupbyvars_n) %>%
  summarise(
    n = sum(!is.na(!!analvar[[1]])),
    mean = mean(!!analvar[[1]], na.rm = TRUE),
    sd = sd(!!analvar[[1]], na.rm = TRUE),
    median = median(!!analvar[[1]], na.rm = TRUE),
    min = min(!!analvar[[1]], na.rm = TRUE),
    max = max(!!analvar[[1]], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(!!!groupbyvars_n)

# combine statistics
vs_summary <- vs_summary %>%
  left_join(n_trt)

# prepare vs_summary in ard format ----
t_ard <- vs_summary %>%
  rename(variable = param, group1 = trt_disp, group1_val = trt_disp) %>%
  tidyr::pivot_longer(
    cols = c(n_trt, n, mean, sd, median, min, max),
    names_to = "stat_name",
    values_to = "stat"
  ) %>%
  mutate(
    stat = as.list(stat),
    stat_label = stat_name,
    fmt_fun = vector("list", n()),
    warning = vector("list", n()),
    error = vector("list", n())
  ) %>%
  as_card(check = TRUE)

# prepare vs_summary for output format ----

display_groupvars <- exprs(param, atpt, trt_disp)

vs_display <- vs_summary %>%
  # --- 1. suppress repeats in grouping variables
  mutate(
    across(c(!!!display_groupvars), suppress_repeats),
    n_trt = ifelse(trt_disp == "", NA_integer_, n_trt)
  ) %>%
  # --- 2. add blank row after last value of column trt_disp
  # use the original var (trt_dispn), not the blanked out one
  insert_blank_rows(col = as_label(exprs(trt_dispn)[[1]])) %>%
  # restrict to columns for table
  # be careful: order in select does reflect columns on output
  select(!!!display_groupvars, n_trt, avisit, n, mean, sd, median, min, max)

# --- 3. Render with gt ---
# define column widths (percentages)
pcw <- c(20, 15, 8, 15)
remainderw <- (100 - sum(pcw[1:4])) / (ncol(vs_display) - 4)

gt_table <- vs_display %>%
  gt() %>%
  cols_label(
    param = "Measure",
    atpt = "Position",
    trt_disp = "Treatment",
    avisit = "Planned Relative Time",
    # labels for statistics
    n_trt = "N",
    n = "n",
    mean = "Mean",
    sd = "SD",
    median = "Median",
    min = "Min.",
    max = "Max."
  ) %>%
  cols_width(
    param ~ pct(pcw[1]),
    atpt ~ pct(pcw[2]),
    trt_disp ~ pct(pcw[3]),
    avisit ~ pct(pcw[4]),
    everything() ~ pct(remainderw)
  ) %>%
  fmt_number(columns = c(mean, median, min, max), decimals = 1) %>%
  fmt_number(columns = c(sd), decimals = 2) %>%
  sub_missing(columns = everything(), missing_text = "") %>%
  cols_align(align = "left", columns = c(param, atpt, trt_disp)) %>%
  cols_align(align = "right", columns = c(n_trt, n, mean, sd, median, min, max)) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_labels(columns = everything())
  ) %>%
  tab_header(
    title    = "Table 14-7.01",
    subtitle = "Summary of Vital Signs at End of Treatment"
  )

# ----------------------------------------------------------------------------
# Output data
# ----------------------------------------------------------------------------

# Save ARD as RDS
saveRDS(t_ard, file.path(path$table_ard, paste0(tableid, ".rds")))

# Save as PDF
# Add header and footer using docorator
gt_table %>%
  as_docorator(
    display_name = tableid,
    tbl_scale = FALSE,
    display_loc = path$table_output,
    header = fancyhead(
      fancyrow(left = "Protocol: CDISCPILOT01", center = NA, right = doc_pagenum()),
      fancyrow(left = "Population: All Subjects", center = NA, right = NA)
    ),
    footer = fancyfoot(
      fancyrow(
        left = paste0(
          "End of treatment is the last on-treatment assessment of the specified vital sign (on or before the Week 24
visit). "
        )
      ),
      fancyrow(left = paste0("Source: ", doc_relative_path()), center = NA, right = doc_datetime())
    ),
    save_object = FALSE
  ) %>%
  render_pdf(path$table_output)
