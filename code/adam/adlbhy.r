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
library(datasetjson)
library(metacore)
library(metatools)

# Load Data ---------------------------------------------------------------

dat_to_load <- list(
  lb = file.path(path$sdtm, "lb.json"),
  supplb = file.path(path$sdtm, "supplb.json"),
  adsl = file.path(path$adam_reference, "adsl.json")
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

lb_hy <- lb %>%
  filter(LBTESTCD %in% c("ALT", "AST", "BILI"))

# Merge with ADSL ---------------------------------------------------------

adlbhy <- lb_hy %>%
  derive_vars_merged(
    dataset_add = adsl,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(TRTSDT, TRTEDT, TRT01P, TRT01PN, TRT01A, TRT01AN, 
                     AGE, AGEGR1, AGEGR1N, RACE, RACEN, SEX,
                     SAFFL, COMP24FL, DSRAEFL)
  ) %>%
  mutate(
    TRTP = TRT01P,
    TRTPN = TRT01PN,
    TRTA = TRT01A,
    TRTAN = TRT01AN,
    PARAMCD = LBTESTCD,
    PARAM = LBTEST,
    PARAMN = case_when(
      PARAMCD == "ALT" ~ 1,
      PARAMCD == "AST" ~ 2,
      PARAMCD == "BILI" ~ 3
    ),
    AVAL = LBSTRESN,
    ADT = as.numeric(as.Date(LBDTC)),
    ADY = LBDY,
    AVISIT = VISIT,
    AVISITN = VISITNUM,
    A1LO = LBSTNRLO,
    A1HI = LBSTNRHI,
    R2A1LO = AVAL / A1LO,
    R2A1HI = AVAL / A1HI,
    PARCAT1 = "CHEM"
  )


# Baseline/Base -----------------------------------------------------------

adlbhy <- adlbhy %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(STUDYID, USUBJID, PARAMCD),
      order = exprs(ADT, LBSEQ),
      new_var = ABLFL,
      mode = "last"
    ),
    filter = ADT <= TRTSDT
  )

adlbhy <- adlbhy %>%
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD),
    source_var = AVAL,
    new_var = BASE,
    filter = ABLFL == "Y"
  ) %>%
  mutate(
    BR2A1LO = BASE / A1LO,
    BR2A1HI = BASE / A1HI
  )


# Hys Law -----------------------------------------------------------------
# BILIHY: BILI > 2*ULN
bilihy <- adlbhy %>%
  filter(PARAMCD == "BILI", R2A1HI > 2) %>%
  mutate(
    PARAMCD = "BILIHY",
    PARAM = "Bilirubin > 2*ULN",
    PARAMN = 4,
    AVAL = 1
  )

# TRANSHY: ALT or AST > 3*ULN
transhy <- adlbhy %>%
  filter(PARAMCD %in% c("ALT", "AST"), R2A1HI > 3) %>%
  group_by(STUDYID, USUBJID, VISIT, VISITNUM, ADT) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    PARAMCD = "TRANSHY",
    PARAM = "Transaminase > 3*ULN",
    PARAMN = 6,
    AVAL = 1
  )

# HYLAW: Both conditions met
hylaw <- adlbhy %>%
  filter(PARAMCD %in% c("ALT", "AST", "BILI")) %>%
  group_by(STUDYID, USUBJID, VISIT, VISITNUM, ADT) %>%
  summarise(
    has_bili = any(PARAMCD == "BILI" & R2A1HI > 2),
    has_trans = any(PARAMCD %in% c("ALT", "AST") & R2A1HI > 3),
    .groups = "drop"
  ) %>%
  filter(has_bili & has_trans) %>%
  left_join(
    adlbhy %>% select(STUDYID, USUBJID, TRTSDT, TRTEDT, TRTP, TRTPN, TRTA, TRTAN,
                      AGE, AGEGR1, AGEGR1N, RACE, RACEN, SEX, SAFFL, COMP24FL, DSRAEFL) %>%
      distinct(),
    by = c("STUDYID", "USUBJID")
  ) %>%
  mutate(
    PARAMCD = "HYLAW",
    PARAM = "Hy's Law",
    PARAMN = 5,
    AVAL = 1,
    PARCAT1 = "CHEM"
  )

# Combine all
adlbhy_final <- bind_rows(adlbhy, bilihy, transhy, hylaw) %>%
  mutate(
    PARAMTYP = if_else(PARAMCD %in% c("BILIHY", "TRANSHY", "HYLAW"), "DERIVED", NA_character_),
    CRIT1FL = if_else(AVAL == 1 & PARAMCD %in% c("BILIHY", "TRANSHY", "HYLAW"), "Y", NA_character_)
  ) %>%
  arrange(STUDYID, USUBJID, PARAMN, ADT)

# Final clean up ----------------------------------------------------------

adlbhy <- adlbhy %>%
  select(STUDYID, SUBJID, USUBJID, TRTP, TRTPN, TRTA, TRTAN, TRTSDT, TRTEDT,
         AGE, AGEGR1, AGEGR1N, RACE, RACEN, SEX, COMP24FL, DSRAEFL, SAFFL,
         AVISIT, AVISITN, ADY, ADT, VISIT, VISITNUM, PARAMTYP, PARAM, PARAMCD, PARAMN,
         PARCAT1, AVAL, BASE, A1LO, A1HI, R2A1LO, R2A1HI, BR2A1LO, BR2A1HI,
         ABLFL, CRIT1FL)

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
  if (attr(adsl[[col]], "format.sas") == "") {
    attr(adsl[[col]], "format.sas") <- NULL
  } else if (attr(adsl[[col]], "format.sas") == "DATE9.") {
    attr(adsl[[col]], "format.sas") <- "DATE9"
  }
}

write_dataset_json_with_metadata(adlbhy, adlbhy_spec, "adlbhy", path$adam_json)
