# ============================================================================
# Program: advs.R
# Purpose: Create ADaM Vital Signs Analysis Dataset (ADVS)
# Description: Derives analysis variables for vital signs data including
#              change from baseline, analysis flags, and parameter mappings.
# Input: SDTM domains (VS), ADaM (ADSL)
# Output: advs.xpt
# ============================================================================

# ----------------------------------------------------------------------------
# SETUP AND LIBRARY LOADING
# ----------------------------------------------------------------------------

# Load required packages
library(admiral)      # ADaM derivations
library(dplyr)        # Data manipulation
library(lubridate)    # Date handling
library(haven)        # Reading/writing SAS datasets
library(stringr)      # String manipulation
library(purrr)        # Functional programming
library(tibble)       # Creating tibbles
library(datasetjson)  # Dataset JSON handling
library(metacore)     # Metadata handling

# ----------------------------------------------------------------------------
# LOAD DATASETS
# ----------------------------------------------------------------------------

# Define data to load
data_to_load_xpt <- list(
  vs = file.path(path$sdtm, "vs.xpt")
)

data_to_load_json <- list(
  adsl = file.path(path$adam, "adsl.json")
)

# Load datasets using map and convert blanks to NA
datasets_xpt <- map(
  data_to_load_xpt,
  ~ convert_blanks_to_na(read_xpt(.x))
)

datasets_json <- map(
  data_to_load_json,
  ~ convert_blanks_to_na(read_dataset_json(.x))
)

# Put datasets into the global environment
list2env(c(datasets_xpt, datasets_json), envir = .GlobalEnv)

# ----------------------------------------------------------------------------
# DERIVATIONS
# ----------------------------------------------------------------------------

## Merge ADSL variables
# Select ADSL variables to merge
adsl_vars <- exprs(
  SITEID, TRTSDT, TRTEDT, TRT01P, TRT01PN, TRT01A, TRT01AN,
  AGE, AGEGR1, AGEGR1N, RACE, RACEN, SEX, SAFFL
)

advs_adsl <- vs %>%
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  mutate(
    TRTP = TRT01P,
    TRTPN = TRT01PN,
    TRTA = TRT01A,
    TRTAN = TRT01AN
  )

## Derive Analysis Date (ADT) and Relative Day (ADY)
advs_dt <- advs_adsl %>%
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = VSDTC
  ) %>%
  derive_vars_dy(
    reference_date = TRTSDT,
    source_vars = exprs(ADT)
  )

## Derive Parameters (PARAMCD, PARAM, PARAMN)
# Create a lookup table for parameters
param_lookup <- tribble(
  ~VSTESTCD, ~PARAMCD, ~PARAM,                               ~PARAMN,
  "SYSBP",   "SYSBP",  "Systolic Blood Pressure (mmHg)",     1,
  "DIABP",   "DIABP",  "Diastolic Blood Pressure (mmHg)",    2,
  "PULSE",   "PULSE",  "Pulse Rate (beats/min)",             3,
  "WEIGHT",  "WEIGHT", "Weight (kg)",                        4,
  "HEIGHT",  "HEIGHT", "Height (cm)",                        5,
  "TEMP",    "TEMP",   "Temperature (C)",                    6,
)

advs_param <- advs_dt %>%
  inner_join(param_lookup, by = "VSTESTCD")

## Derive Analysis Value
advs_aval <- advs_param %>%
  mutate(AVAL = VSSTRESN)

## Derive Analysis Timepoints (ATPT, ATPTN)
atpt_lookup <- tribble(
  ~ATPT,                              ~ATPTN,
  "AFTER LYING DOWN FOR 5 MINUTES",   815,
  "AFTER STANDING FOR 1 MINUTE",      816,
  "AFTER STANDING FOR 3 MINUTES",     817,
  NA,                                NA_real_
)

advs_atpt <- advs_aval %>%
  mutate(ATPT = VSTPT) %>%
  left_join(atpt_lookup, by = "ATPT")

## Derive Analysis Visits (AVISIT, AVISITN)
# Lookup table for visits
avisit_lookup <- tribble(
  ~AVISIT,            ~AVISITN,
  "Baseline",         0,
  "Week 2",           2,
  "Week 4",           4,
  "Week 6",           6,
  "Week 8",           8,
  "Week 12",          12,
  "Week 16",          16,
  "Week 20",          20,
  "Week 24",          24,
  "Week 26",          26,
  "End of Treatment", 99,
  NA,                NA_real_
)

advs_avisit <- advs_atpt %>%
  mutate(
    AVISIT = case_when(
      str_detect(VISIT, "SCREENING|UNSCH|AMBUL ECG|RETRIEVAL") ~ NA_character_,
      TRUE ~ str_to_title(VISIT)
    )
  ) %>%
  left_join(avisit_lookup, by = "AVISIT")

# Derive end of treatment visit
advs_eot <- advs_avisit %>%
  derive_extreme_records(
    dataset_add = advs_avisit,
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, ATPT),
    order = exprs(ADT, AVISITN, ATPTN, AVAL),
    mode = "last",
    filter_add = (AVISITN > 2),
    set_values_to = exprs(
      AVISIT = "End of Treatment",
      AVISITN = 99
    )
  )

## Derive Baseline Variables (ABLFL, BASE, BASETYPE, CHG, PCHG)
advs_baseline <- advs_eot %>%
  mutate(
    ABLFL = VSBLFL,
    BASETYPE = VSTPT
  ) %>%
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVAL,
    new_var = BASE,
    filter = ABLFL == "Y"
  ) %>%
  derive_var_chg() %>%
  derive_var_pchg()

## Derive Analysis Flags (ANL01FL)
advs_anl01fl <- advs_baseline %>%
  mutate(
    ANL01FL = case_when(
      !is.na(AVISIT) ~ "Y",
      TRUE ~ NA_character_
    )
  )

# Sort the data
sort_order <- advs_anl01fl %>%
  select(USUBJID, PARAMCD, ATPTN, VISITNUM, AVISIT, VSSEQ) %>%
  colnames()

advs_sorted <- advs_anl01fl %>%
  arrange(across(all_of(sort_order)))

# ----------------------------------------------------------------------------
# EXPORT
# ----------------------------------------------------------------------------

# Load define.xml metadata;
advs_spec <- define_to_metacore(
  file.path(path$adam, "define.xml"),
  quiet = TRUE
) %>%
  select_dataset("ADVS")

# Select variables as per spec
advs_final <- advs_sorted %>%
  select(
    advs_spec$ds_vars$variable
  )

# Prepare column metadata
oid_cols <- advs_spec$ds_vars %>%
  select(dataset, variable, key_seq) %>%
  left_join(advs_spec$var_spec, by = c("variable")) %>%
  rename(name = variable, dataType = type, keySequence = key_seq, displayFormat = format) %>%
  mutate(itemOID = paste0("IT.", dataset, ".", name)) %>%
  select(itemOID, name, label, dataType, length, keySequence, displayFormat) %>%
  mutate(
    dataType =
      case_when(
        displayFormat == "DATE9." ~ "date",
        displayFormat == "DATETIME20." ~ "datetime",
        substr(name, nchar(name) - 3 + 1, nchar(name)) == "DTC" & length == "8" ~ "date",
        substr(name, nchar(name) - 3 + 1, nchar(name)) == "DTC" & length == "20" ~ "datetime",
        dataType == "text" ~ "string",
        .default = as.character(dataType)
      ),
    targetDataType =
      case_when(
        displayFormat == "DATE9." ~ "integer",
        displayFormat == "DATETIME20." ~ "integer",
        .default = NA
      ),
    length = case_when(
      dataType == "string" ~ length,
      .default = NA
    )
  ) %>%
  data.frame()

# Create and write dataset JSON
dataset_json(advs_final,
  last_modified = strftime(as.POSIXlt(Sys.time(), "UTC"), "%Y-%m-%dT%H:%M"),
  originator = "R Submission Pilot 6",
  sys = paste0("R on ", R.Version()$os, " ", unname(Sys.info())[[2]]),
  sys_version = R.Version()$version.string,
  version = "1.1.0",
  study = "Pilot 6",
  metadata_version = "MDV.TDF_ADaM.ADaM-IG.1.1", # from define
  metadata_ref = file.path(path$adam, "define.xml"),
  item_oid = paste0("IG.ADVS"),
  name = "ADVS",
  dataset_label = advs_spec$ds_spec[["label"]],
  file_oid = file.path(path$adam, "advs.json"),
  columns = oid_cols
) %>%
  write_dataset_json(file = file.path(path$adam, "advs.json"), float_as_decimals = FALSE)