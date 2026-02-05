#************************************************************************
# Purpose:     Generate ADLBHY dataset
# Input:       LB (from datasetjson)
# Output:      adlbhy.json
#************************************************************************

# Note to Reviewer
# To rerun the code below, please refer ADRG appendix.
# After required package are installed, the path variable needs to be defined
# in the .Rprofile file

# Libraries ---------------------------------------------------------------

library(haven)
library(dplyr)
library(admiral)
library(purrr)
library(tidyr)
library(datasetjson)
library(metacore)
library(metatools)
library(xportr)

# Load Data ---------------------------------------------------------------

dat_to_load <- list(
  adlbc = file.path(path$adam_reference, "adlbc.json")
)

datasets <- map(
  dat_to_load,
  ~ convert_blanks_to_na(read_dataset_json(.x, decimals_as_floats = TRUE))
)

list2env(datasets, envir = .GlobalEnv)

# Metadata ----------------------------------------------------------------
## Very noisy function - remove suppress if you want to see warnings
adlbhy_spec <- suppressWarnings(
  spec_to_metacore(
    file.path(path$adam_reference, "pilot6-specs.xlsx"),
    where_sep_sheet = FALSE,
    quiet = TRUE
    )
  ) %>%
  select_dataset("ADLBHY")

# Select Parameters from LB -----------------------------------------------
allowed_avisit <- c("Baseline", "Week 2", "Week 4", "Week 6", "Week 8",
                    "Week 12", "Week 16", "Week 20", "Week 24")


adlbhy <- adlbc %>%
  filter(PARAMCD %in% c("ALT", "AST", "BILI"),
         AVISIT %in% allowed_avisit) %>%
  mutate(
    PARAMN = case_when(
      PARAMCD == "ALT" ~ 1,
      PARAMCD == "AST" ~ 2,
      PARAMCD == "BILI" ~ 3
    ),
    PARAMTYP = NA_character_,
    CRIT1 = NA_character_,
    CRIT1FL = NA_character_,
    SHIFT1 = NA_character_,
    SHIFT1N = NA_integer_
  )

# Hys Law -----------------------------------------------------------------
# BILIHY: BILI > 2*ULN
bilihy <- adlbhy %>%
  filter(PARAMCD == "BILI") %>%
  mutate(
    PARAMCD = "BILIHY",
    PARAM = "Bilirubin 1.5 x ULN",
    PARAMN = 4,
    AVAL = if_else(R2A1HI > 1.5, 1, 0),
    PARAMTYP = "DERIVED",
    CRIT1 = "R2A1HI > 1.5",
    CRIT1FL = if_else(R2A1HI > 1.5, "Y", "N"),
    PARCAT1 = "HYLAW"
  )

# TRANSHY: ALT or AST > 3*ULN
transhy <- adlbhy %>%
  filter(PARAMCD %in% c("ALT", "AST")) %>%
  group_by(STUDYID, USUBJID, VISIT, VISITNUM, ADT) %>%
  summarise(
    has_elevated = any(R2A1HI > 1.5),
    .groups = "drop"
  ) %>%
  left_join(
    adlbhy %>% 
      filter(PARAMCD %in% c("ALT", "AST")) %>%
      select(-PARAMCD, -PARAM, -PARAMN, -AVAL) %>%
      distinct(),
    by = c("STUDYID", "USUBJID", "VISIT", "VISITNUM", "ADT")
  ) %>%
  mutate(
    PARAMCD = "TRANSHY",
    PARAM = "Transaminase 1.5 x ULN",
    PARAMN = 6,
    AVAL = if_else(has_elevated, 1, 0),
    PARAMTYP = "DERIVED",
    CRIT1 = "R2A1HI > 1.5",
    CRIT1FL = if_else(has_elevated, "Y", "N"),
    PARCAT1 = "HYLAW"
  ) %>%
  select(-has_elevated)

# HYLAW: Both conditions met  (0/1 for all visits)
hylaw_visits <- adlbhy %>%
  filter(PARAMCD %in% c("ALT", "AST", "BILI")) %>%
  group_by(STUDYID, USUBJID, AVISIT, AVISITN, ADT) %>%
  summarise(
    has_bili = any(PARAMCD == "BILI" & R2A1HI > 1.5),
    has_trans = any(PARAMCD %in% c("ALT", "AST") & R2A1HI > 1.5),
    .groups = "drop"
  ) %>%
  left_join(
    adlbhy %>% 
      select(STUDYID, USUBJID, TRTP, TRTPN, TRTA, TRTAN, TRTSDT, TRTEDT,
             AGE, AGEGR1, AGEGR1N, RACE, RACEN, SEX, COMP24FL, DSRAEFL, SAFFL) %>% 
      distinct(),
    by = c("STUDYID", "USUBJID")
  ) %>%
  mutate(
    PARAMCD = "HYLAW",
    PARAM = "Total Bili 1.5 x ULN and Transaminase 1.5 x ULN",
    PARAMN = 5,
    AVAL = if_else(has_bili & has_trans, 1, 0),
    PARCAT1 = "HYLAW",
    PARAMTYP = "DERIVED",
    CRIT1 = "R2A1HI > 1.5",
    CRIT1FL = if_else(has_bili & has_trans, "Y", "N"),
    ADY = NA_integer_,
    A1LO = NA_real_,
    A1HI = NA_real_,
    R2A1LO = NA_real_,
    R2A1HI = NA_real_,
    BASE = NA_real_,
    BR2A1LO = NA_real_,
    BR2A1HI = NA_real_,
    ABLFL = NA_character_,
    SHIFT1 = NA_character_,
    SHIFT1N = NA_integer_,
    LBSEQ = NA_integer_,
    SUBJID = sub(".*-", "", USUBJID)
  ) %>%
  select(-has_bili, -has_trans)

# Combine all
adlbhy <- bind_rows(adlbhy, bilihy, transhy, hylaw_visits) %>%
  arrange(STUDYID, USUBJID, PARAMN, ADT)

# Final clean up ----------------------------------------------------------

adlbhy <- adlbhy %>%
  select(STUDYID, SUBJID, USUBJID, TRTP, TRTPN, TRTA, TRTAN, TRTSDT, TRTEDT,
         AGE, AGEGR1, AGEGR1N, RACE, RACEN, SEX, COMP24FL, DSRAEFL, SAFFL,
         AVISIT, AVISITN, ADY, ADT, VISIT, VISITNUM, PARAMTYP, PARAM, PARAMCD, PARAMN,
         PARCAT1, AVAL, BASE, A1LO, A1HI, R2A1LO, R2A1HI, BR2A1LO, BR2A1HI,
         ABLFL, SHIFT1, SHIFT1N, CRIT1, CRIT1FL, LBSEQ)

# Output ------------------------------------------------------------------

adlbhy <- adlbhy %>%
  drop_unspec_vars(adlbhy_spec) %>%
  check_ct_data(adlbhy_spec, na_acceptable = TRUE) %>%
  order_cols(adlbhy_spec) %>%
  sort_by_key(adlbhy_spec) %>%
  set_variable_labels(adlbhy_spec) %>%
  xportr_label(adlbhy_spec) %>%
  xportr_df_label(adlbhy_spec, domain = "adlbhy") %>%
  xportr_format(
    adlbhy_spec$var_spec %>% mutate_at(c("format"), ~ replace_na(., "")),
    "ADLBHY"
  ) %>%
  convert_na_to_blanks()

# FIX: attribute issues where sas.format attributes set to DATE9. are changed to DATE9,
# and missing formats are set to NULL (instead of an empty character vector)
# when reading original xpt file
for (col in colnames(adlbhy)) {
  if (attr(adlbhy[[col]], "format.sas") == "") {
    attr(adlbhy[[col]], "format.sas") <- NULL
  } else if (attr(adlbhy[[col]], "format.sas") == "DATE9.") {
    attr(adlbhy[[col]], "format.sas") <- "DATE9"
  }
}

write_dataset_json_with_metadata(adlbhy, adlbhy_spec, "adlbhy", path$adam_json)
