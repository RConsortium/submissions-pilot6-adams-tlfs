# ============================================================================
# Program: adlbh.R
# Purpose: Create ADaM Laboratory Data - Hematology Dataset (ADLBH)
# Description: Derives analysis variables for hematology laboratory data
#              including change from baseline, analysis flags, and parameter
#              mappings. Includes normal range indicators and shift analysis.
# Input: SDTM domains (LB), ADaM (ADSL)
# Output: adlbh.json
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
# LOAD METADATA
# ----------------------------------------------------------------------------

# Load define.xml metadata
adlbh_spec <- define_to_metacore(
  file.path(path$adam, "define.xml"),
  quiet = TRUE
) %>%
  select_dataset("ADLBH")

# Tibbles with codelists
avisit_codelist <- adlbh_spec$codelist %>%
  filter(code_id == "CL.ADLBH.AVISITN") %>%
  pull(codes) %>%
  pluck(1) %>%
  select(code, decode) %>%
  mutate(
    code = str_trim(code),
    code = as.integer((code))
  ) %>%
  rename(AVISITN = code, AVISIT = decode)

paramcd_codelist <- adlbh_spec$codelist %>%
  filter(code_id == "CL.PARAMCD_ADLBH") %>%
  pull(codes) %>%
  pluck(1) %>%
  select(code, decode) %>%
  rename(PARAMCD = code, PARAM = decode)

paramn_codelist <- adlbh_spec$codelist %>%
  filter(code_id == "CL.PARAMN_ADLBH") %>%
  pull(codes) %>%
  pluck(1) %>%
  select(code, decode) %>%
  mutate(
    code = str_trim(code),
    code = as.integer((code))
  ) %>%
  rename(PARAMN = code, PARAM = decode)

param_lookup <- paramcd_codelist %>%
  inner_join(paramn_codelist, by = "PARAM") %>%
  select(PARAMCD, PARAM, PARAMN) %>%
  mutate(LBTESTCD = PARAMCD)

# ----------------------------------------------------------------------------
# LOAD DATASETS
# ----------------------------------------------------------------------------

# Define data to load
data_to_load_xpt <- list(
  lb = file.path(path$sdtm, "lb.xpt")
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

# Filter to hematology tests only
lb_hematology <- lb %>%
  filter(LBCAT == "HEMATOLOGY")

## Merge ADSL variables
# Select ADSL variables to merge
adsl_vars <- exprs(
  SITEID, TRTSDT, TRTEDT, TRT01P, TRT01PN, TRT01A, TRT01AN,
  AGE, AGEGR1, AGEGR1N, RACE, RACEN, SEX, SAFFL
)

adlbh_adsl <- lb_hematology %>%
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
adlbh_dt <- adlbh_adsl %>%
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = LBDTC
  ) %>%
  derive_vars_dy(
    reference_date = TRTSDT,
    source_vars = exprs(ADT)
  )

## Derive Parameters (PARAMCD, PARAM, PARAMN)
adlbh_param <- adlbh_dt %>%
  inner_join(param_lookup, by = "LBTESTCD")

## Derive Analysis Value (AVAL)
adlbh_aval <- adlbh_param %>%
  mutate(
    AVAL = LBSTRESN,
    AVALC = LBSTRESC
  )

## Derive Analysis Visits (AVISIT, AVISITN)
# Map visits to analysis visits
adlbh_avisit <- adlbh_aval %>%
  mutate(
    AVISIT = case_when(
      str_detect(VISIT, "SCREENING|UNSCH|AMBUL|RETRIEVAL") ~ NA_character_,
      TRUE ~ str_to_title(VISIT)
    )
  ) %>%
  left_join(avisit_codelist, by = "AVISIT")

# Derive end of treatment visit
adlbh_eot <- adlbh_avisit %>%
  derive_extreme_records(
    dataset_add = adlbh_avisit,
    by_vars = exprs(STUDYID, USUBJID, PARAMCD),
    order = exprs(ADT, AVISITN, AVAL),
    mode = "last",
    filter_add = (AVISITN > 2),
    set_values_to = exprs(
      AVISIT = "End of Treatment",
      AVISITN = 99
    )
  )

## Derive Baseline Variables (ABLFL, BASE, CHG, PCHG)
adlbh_baseline <- adlbh_eot %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(STUDYID, USUBJID, PARAMCD),
      order = exprs(ADT, LBSEQ),
      new_var = ABLFL,
      mode = "last"
    ),
    filter = (!is.na(AVAL) & ADT <= TRTSDT & !is.na(AVISITN))
  ) %>%
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD),
    source_var = AVAL,
    new_var = BASE,
    filter = ABLFL == "Y"
  ) %>%
  derive_var_chg() %>%
  derive_var_pchg()

## Derive Normal Range Indicators (ANRLO, ANRHI, ANRIND, BNRIND)
adlbh_nrange <- adlbh_baseline %>%
  mutate(
    ANRLO = LBSTNRLO,
    ANRHI = LBSTNRHI,
    ANRIND = case_when(
      !is.na(AVAL) & !is.na(ANRLO) & !is.na(ANRHI) & AVAL < ANRLO ~ "LOW",
      !is.na(AVAL) & !is.na(ANRLO) & !is.na(ANRHI) & AVAL > ANRHI ~ "HIGH",
      !is.na(AVAL) & !is.na(ANRLO) & !is.na(ANRHI) ~ "NORMAL",
      TRUE ~ NA_character_
    )
  ) %>%
  derive_vars_merged(
    dataset_add = .,
    by_vars = exprs(STUDYID, USUBJID, PARAMCD),
    filter_add = ABLFL == "Y",
    new_vars = exprs(BNRIND = ANRIND)
  )

## Derive Shift Variables
adlbh_shift <- adlbh_nrange %>%
  mutate(
    SHIFT1 = case_when(
      !is.na(BNRIND) & !is.na(ANRIND) ~ paste0(BNRIND, " to ", ANRIND),
      TRUE ~ NA_character_
    )
  )

## Derive Analysis Flags (ANL01FL)
adlbh_anl01fl <- adlbh_shift %>%
  mutate(
    ANL01FL = case_when(
      !is.na(AVISIT) & !is.na(AVAL) ~ "Y",
      TRUE ~ NA_character_
    )
  )

## Derive ATOXGR and BTOXGR (Toxicity Grade)
# Common Toxicity Criteria for Hematology Parameters
# These are example grades and should be customized based on specific study requirements
adlbh_toxicity <- adlbh_anl01fl %>%
  mutate(
    ATOXGR = case_when(
      # Hemoglobin (g/dL)
      PARAMCD == "HGB" & !is.na(AVAL) & AVAL >= ANRLO ~ "0",
      PARAMCD == "HGB" & !is.na(AVAL) & AVAL < ANRLO & AVAL >= 10.0 ~ "1",
      PARAMCD == "HGB" & !is.na(AVAL) & AVAL < 10.0 & AVAL >= 8.0 ~ "2",
      PARAMCD == "HGB" & !is.na(AVAL) & AVAL < 8.0 & AVAL >= 6.5 ~ "3",
      PARAMCD == "HGB" & !is.na(AVAL) & AVAL < 6.5 ~ "4",
      # Leukocytes (10^9/L)
      PARAMCD == "WBC" & !is.na(AVAL) & AVAL >= ANRLO ~ "0",
      PARAMCD == "WBC" & !is.na(AVAL) & AVAL < ANRLO & AVAL >= 3.0 ~ "1",
      PARAMCD == "WBC" & !is.na(AVAL) & AVAL < 3.0 & AVAL >= 2.0 ~ "2",
      PARAMCD == "WBC" & !is.na(AVAL) & AVAL < 2.0 & AVAL >= 1.0 ~ "3",
      PARAMCD == "WBC" & !is.na(AVAL) & AVAL < 1.0 ~ "4",
      # Platelets (10^9/L)
      PARAMCD == "PLAT" & !is.na(AVAL) & AVAL >= ANRLO ~ "0",
      PARAMCD == "PLAT" & !is.na(AVAL) & AVAL < ANRLO & AVAL >= 75.0 ~ "1",
      PARAMCD == "PLAT" & !is.na(AVAL) & AVAL < 75.0 & AVAL >= 50.0 ~ "2",
      PARAMCD == "PLAT" & !is.na(AVAL) & AVAL < 50.0 & AVAL >= 25.0 ~ "3",
      PARAMCD == "PLAT" & !is.na(AVAL) & AVAL < 25.0 ~ "4",
      # Neutrophils (10^9/L)
      PARAMCD == "NEUT" & !is.na(AVAL) & AVAL >= ANRLO ~ "0",
      PARAMCD == "NEUT" & !is.na(AVAL) & AVAL < ANRLO & AVAL >= 1.5 ~ "1",
      PARAMCD == "NEUT" & !is.na(AVAL) & AVAL < 1.5 & AVAL >= 1.0 ~ "2",
      PARAMCD == "NEUT" & !is.na(AVAL) & AVAL < 1.0 & AVAL >= 0.5 ~ "3",
      PARAMCD == "NEUT" & !is.na(AVAL) & AVAL < 0.5 ~ "4",
      # Lymphocytes (10^9/L)
      PARAMCD == "LYMPH" & !is.na(AVAL) & AVAL >= ANRLO ~ "0",
      PARAMCD == "LYMPH" & !is.na(AVAL) & AVAL < ANRLO & AVAL >= 0.8 ~ "1",
      PARAMCD == "LYMPH" & !is.na(AVAL) & AVAL < 0.8 & AVAL >= 0.5 ~ "2",
      PARAMCD == "LYMPH" & !is.na(AVAL) & AVAL < 0.5 & AVAL >= 0.2 ~ "3",
      PARAMCD == "LYMPH" & !is.na(AVAL) & AVAL < 0.2 ~ "4",
      TRUE ~ NA_character_
    )
  ) %>%
  derive_vars_merged(
    dataset_add = .,
    by_vars = exprs(STUDYID, USUBJID, PARAMCD),
    filter_add = ABLFL == "Y",
    new_vars = exprs(BTOXGR = ATOXGR)
  )

## Derive Criterion Flags (CRIT1, CRIT1FL)
adlbh_criterion <- adlbh_toxicity %>%
  mutate(
    CRIT1 = "Abnormal",
    CRIT1FL = case_when(
      ANRIND %in% c("HIGH", "LOW") ~ "Y",
      TRUE ~ NA_character_
    )
  )

# Sort the data
sort_order <- adlbh_criterion %>%
  select(USUBJID, PARAMCD, AVISITN, VISITNUM, ADT, LBSEQ) %>%
  colnames()

adlbh_sorted <- adlbh_criterion %>%
  arrange(across(all_of(sort_order)))

# ----------------------------------------------------------------------------
# EXPORT
# ----------------------------------------------------------------------------

# Select variables as per spec
adlbh_final <- adlbh_sorted %>%
  select(
    adlbh_spec$ds_vars$variable
  )

# Prepare column metadata
oid_cols <- adlbh_spec$ds_vars %>%
  select(dataset, variable, key_seq) %>%
  left_join(adlbh_spec$var_spec, by = c("variable")) %>%
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
dataset_json(adlbh_final,
  last_modified = strftime(as.POSIXlt(Sys.time(), "UTC"), "%Y-%m-%dT%H:%M"),
  originator = "R Submission Pilot 6",
  sys = paste0("R on ", R.Version()$os, " ", unname(Sys.info())[[2]]),
  sys_version = R.Version()$version.string,
  version = "1.1.0",
  study = "Pilot 6",
  metadata_version = "MDV.TDF_ADaM.ADaM-IG.1.1", # from define
  metadata_ref = file.path(path$adam, "define.xml"),
  item_oid = paste0("IG.ADLBH"),
  name = "ADLBH",
  dataset_label = adlbh_spec$ds_spec[["label"]],
  file_oid = file.path(path$adam, "adlbh.json"),
  columns = oid_cols
) %>%
  write_dataset_json(file = file.path(path$adam, "adlbh.json"), float_as_decimals = FALSE)

# Print summary
cat("\n============================================================================\n")
cat("ADLBH Dataset Creation Complete\n")
cat("============================================================================\n")
cat("Output file:", file.path(path$adam, "adlbh.json"), "\n")
cat("Number of records:", nrow(adlbh_final), "\n")
cat("Number of subjects:", length(unique(adlbh_final$USUBJID)), "\n")
cat("Number of parameters:", length(unique(adlbh_final$PARAMCD)), "\n")
cat("============================================================================\n")
