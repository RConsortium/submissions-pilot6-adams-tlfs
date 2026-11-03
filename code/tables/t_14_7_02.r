# Load required packages
library(dplyr)
library(purrr)
library(datasetjson)
library(cards)
library(gtsummary)
library(gt)
library(docorator)

# Import utility functions
source(file.path("code", "utils", "doc_relative_path.r"))

# ----------------------------------------------------------------------------
# Local script functions
# ----------------------------------------------------------------------------

convert_blanks_to_na_local <- function(df) {
  df %>% mutate(across(where(is.character), ~ na_if(.x, "")))
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

# Suppress repeating values in `suppress_vars` by replacing with blank string.
suppress_repeats <- function(df, suppress_vars, mirror_vars = NULL) {
  df <- df %>%
    mutate(across(
      any_of(suppress_vars),
      ~ {
        x <- as.character(.x)
        last_val <- ""
        result <- character(length(x))
        for (j in seq_along(x)) {
          if (is.na(x[j])) {
            result[j] <- x[j]
          } else if (x[j] == last_val) {
            result[j] <- ""
          } else {
            result[j] <- x[j]
            last_val <- x[j]
          }
        }
        result
      }
    ))
  if (!is.null(mirror_vars) && mirror_vars %in% names(df)) {
    mirror_col <- mirror_vars
    ref_col <- suppress_vars[length(suppress_vars)]
    df <- df %>% mutate(across(
      all_of(mirror_col),
      ~ if_else(.data[[ref_col]] == "", "", as.character(.x))
    ))
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
# Prepare data
# ----------------------------------------------------------------------------
tableid <- "t_14_7_02"
analysis_var <- "chg"
group_vars <- c("param", "atpt", "trt_disp", "avisit")
group_var_labels <- c("Measure", "Position", "Treatment", "Planned Relative Time")
blank_row_var <- "trt_disp"

trt_relabel <- list(
  "Xanomeline Low Dose" ~ "Xan.Low",
  "Xanomeline High Dose" ~ "Xan.High"
)

visit_relabel <- list("End of Treatment" ~ "End of Trt.")

adsl <- adsl %>%
  filter(saffl == "Y") %>%
  mutate(trt_disp = replace_values(trt01a, !!!trt_relabel)) %>%
  mutate(trt_disp = factor(trt_disp, levels = unique(trt_disp[order(trt01an)])))

advs <- advs %>%
  filter(
    saffl == "Y",
    paramcd %in% c("SYSBP", "DIABP", "PULSE"),
    avisit %in% c("Week 24", "End of Treatment")
  ) %>%
  mutate(
    trt_disp = replace_values(trta, !!!trt_relabel),
    trt_disp = factor(trt_disp, levels = unique(trt_disp[order(trtan)])),
    param = factor(param, levels = unique(param[order(paramn)])),
    atpt = factor(atpt, levels = unique(atpt[order(atptn)])),
    avisit = replace_values(avisit, !!!visit_relabel),
    avisit = factor(avisit, levels = unique(avisit[order(avisitn)]))
  )

# ----------------------------------------------------------------------------
# Treatment-level N for header
# ----------------------------------------------------------------------------
n_trt <- adsl %>%
  count(trt_disp, name = "n_trt") %>%
  arrange(trt_disp)

# ----------------------------------------------------------------------------
# Build table: one tbl_wide_summary per param x atpt x trt_disp x avisit
# ----------------------------------------------------------------------------
group_keys <- advs %>%
  group_by(across(all_of(group_vars))) %>%
  group_keys() %>%
  arrange(across(all_of(group_vars))) %>%
  left_join(n_trt, by = "trt_disp")

statlist <- c("{N_nonmiss}", "{mean}", "{sd}", "{median}", "{min}", "{max}")
statlist_labels <- c("n", "Mean", "SD", "Median", "Min.", "Max.")

tbls <- advs %>%
  arrange(across(all_of(group_vars))) %>%
  group_by(across(all_of(group_vars))) %>%
  group_split() %>%
  imap(
    function(slice, i) {
      tbl_wide_summary(
        slice,
        include = all_of(analysis_var),
        statistic = statlist,
        digits = setNames(list(list(N_nonmiss = 0, mean = 1, sd = 2, median = 1, min = 1, max = 1)), analysis_var),
        label = setNames(list(" "), analysis_var)
      ) %>%
        modify_table_body(
          ~ .x %>%
            bind_cols(group_keys[i, c(group_vars, "n_trt"), drop = FALSE]) %>%
            relocate(all_of(group_vars), .before = label) %>%
            relocate(n_trt, .after = trt_disp)
        ) %>%
        modify_column_unhide(all_of(c(group_vars, "n_trt"))) %>%
        modify_column_hide(label) %>%
        modify_header(
          !!!c(
            setNames(as.list(group_var_labels), group_vars),
            list(n_trt = "N"),
            setNames(as.list(statlist_labels), paste0("stat_", seq_along(statlist_labels)))
          )
        )
    }
  )

# define column widths (percentages)
# ensure pcw has same length as group_vars
pcw <- c(20, 15, 10, 15)
# extra 1 after length(statlist) is due to column n_trt
remainderw <- (100 - sum(pcw[1:4])) / (length(statlist) + 1)

gtsum_table <- tbl_stack(tbls,
  tbl_ids = apply(group_keys[, group_vars], 1, paste, collapse = " | ")
) %>%
  modify_table_body(~ insert_blank_rows_by(.x, group_var = !!sym(blank_row_var))) %>%
  modify_table_body(~ suppress_repeats(.x, suppress_vars = head(group_vars, -1), mirror_vars = "n_trt")) %>%
  modify_column_alignment(columns = all_of(group_vars), align = "left") %>%
  modify_column_alignment(columns = "n_trt", align = "right") %>%
  modify_column_alignment(
    columns = all_of(paste0("stat_", seq_along(statlist_labels))),
    align   = "right"
  )

gt_table <- gtsum_table %>%
  as_gt() %>%
  cols_width(
    .list = c(
      map2(group_vars, pcw, ~ as.formula(paste0(.x, " ~ pct(", .y, ")"))),
      list(n_trt ~ pct(5)),
      list(everything() ~ pct(remainderw))
    )
  ) %>%
  tab_header(
    title = "Table 14-7.02",
    subtitle = "Summary of Vital Signs Change from Baseline at End of Treatment"
  )

# ----------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------
ard_out <- gather_ard(gtsum_table)

saveRDS(ard_out, file.path(path$table_ard, paste0(tableid, ".rds")))

# Save as PDF via docorator
gt_table %>%
  as_docorator(
    display_name = tableid,
    tbl_scale = FALSE,
    display_loc = path$table_output,
    header = fancyhead(
      fancyrow(left = "Protocol: CDISCPILOT01", center = NA, right = doc_pagenum()),
      fancyrow(left = "Population: Safety", center = NA, right = NA)
    ),
    footer = fancyfoot(
      fancyrow(
        left = paste0(
          "End of treatment is the last on-treatment assessment of the specified ",
          "vital sign (on or before the Week 24 visit). "
        )
      ),
      fancyrow(left = paste0("Source: ", doc_relative_path()), center = NA, right = doc_datetime())
    ),
    save_object = FALSE
  ) %>%
  render_pdf(path$table_output)
