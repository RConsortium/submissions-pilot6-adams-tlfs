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
library(metatools)

# ----------------------------------------------------------------------------
# LOAD METADATA
# ----------------------------------------------------------------------------

#googlesheets4::read_sheet("https://docs.google.com/spreadsheets/u/1/d/e/2PACX-1vR5B7pDzgiHFCjFnmmBVNxytwrH5A06iahdIfkN_WnXf7eoRvdeoUWgJvhImvnn3eJE1DfUq9S2CadT/pubhtml#gid=826399612")
 
## Load dataset specs -------------
# Very noisy function - remove suppress if you want to see warning(s)
metacore <- suppressWarnings(
  spec_to_metacore(
    file.path(path$adam_reference, "pilot6-specs.xlsx"),
    where_sep_sheet = FALSE,
    quiet = TRUE
  )
)

# Get the specifications for the dataset we are currently building
adlbh_spec <- metacore %>%
  select_dataset("ADLBH")

# Tibbles with codelists
avisit_codelist <- adlbh_spec$codelist %>%
  filter(code_id == "ADLBH.AVISITN") %>%
  pull(codes) %>%
  pluck(1) %>%
  select(code, decode) %>%
  mutate(
    code = str_trim(code),
    code = as.integer((code))
  ) %>%
  rename(AVISITN = code, AVISIT = decode)

paramcd_codelist <- adlbh_spec$codelist %>%
  filter(code_id == "PARAMCD_ADLBH") %>%
  pull(codes) %>%
  pluck(1) %>%
  select(code, decode) %>%
  rename(PARAMCD = code, PARAM = decode)

paramn_codelist <- adlbh_spec$codelist %>%
  filter(code_id == "PARAMN_ADLBH") %>%
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

# Load SDTM datasets
lb <- read_dataset_json(file.path(path$sdtm, "lb.json"))
adsl <- read_dataset_json(file.path(path$adam, "adsl.json"))

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
  derive_vars_merged_lookup(
    dataset_add = param_lookup, 
    by_vars = exprs(LBTESTCD)
    ) %>% 
  mutate(PARCAT1 = str_to_title(LBCAT))

## Derive Analysis Value (AVAL)
adlbh_aval <- adlbh_param %>%
  mutate(
    AVAL = as.numeric(LBSTRESN), # this needs to be looked at
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
  derive_vars_merged_lookup(
    dataset_add = avisit_codelist, 
    by_vars = exprs(AVISIT)
    )

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
  mutate(ABLFL = LBBLFL) %>% 
  # restrict_derivation(
  #   derivation = derive_var_extreme_flag,
  #   args = params(
  #     by_vars = exprs(STUDYID, USUBJID, PARAMCD),
  #     order = exprs(ADT, LBSEQ),
  #     new_var = ABLFL,
  #     mode = "last"
  #   ),
  #   filter = (!is.na(AVAL) & ADT <= TRTSDT & !is.na(AVISITN))
  # ) %>%
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
    A1LO = if_else(!is.na(LBSTNRLO), as.numeric(LBSTNRLO), NA),
    A1HI = if_else(!is.na(LBSTNRHI), as.numeric(LBSTNRHI), NA),
    R2A1LO = if_else(!is.na(LBSTNRLO), AVAL / A1LO, NA),
    R2A1HI = if_else(!is.na(LBSTNRHI), AVAL / A1HI, NA),
    BR2A1LO = if_else(!is.na(LBSTNRLO) & ABLFL == "Y", AVAL / A1HI, NA),
    BR2A1HI = if_else(!is.na(LBSTNRHI) & ABLFL == "Y", AVAL / A1HI, NA),
    ANRIND = case_when(
      !is.na(AVAL) & !is.na(A1LO) & !is.na(A1HI) & AVAL < A1LO ~ "LOW",
      !is.na(AVAL) & !is.na(A1LO) & !is.na(A1HI) & AVAL > A1HI ~ "HIGH",
      !is.na(AVAL) & !is.na(A1LO) & !is.na(A1HI) ~ "NORMAL",
      TRUE ~ NA_character_
    )
  )  


adlbh_bnrind <- derive_var_base(
  adlbh_nrange,
  by_vars = exprs(STUDYID, USUBJID, PARAMCD),
  source_var = ANRIND,
  new_var = BNRIND
) %>%
  mutate(VISITNUM = as.numeric(VISITNUM))

adlbh_max_aval <- adlbh_bnrind %>%
  mutate(VISITNUM = as.numeric(VISITNUM)) %>% 
  filter(VISITNUM >= 4 & VISITNUM <= 12) %>% 
  group_by(USUBJID, PARAMCD) %>%         
  summarise(
    max_aval_in_range = max(AVAL, na.rm = TRUE),
    .groups = 'drop'                            
  ) %>% 
  mutate(ANL01FL = "Y",
         AVAL = max_aval_in_range) 

adlbh_anl01fl <- left_join(adlbh_bnrind, adlbh_max_aval)

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
