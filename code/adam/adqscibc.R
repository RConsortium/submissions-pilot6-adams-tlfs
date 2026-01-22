# ============================================================================
# Program: adqscibc.R
# Purpose: Create ADaM CIBIC+ Analysis Dataset (ADQSCIBC)
# Description: Derives analysis variables for CIBIC+ Dataset.
# Input: SDTM domains (QS), ADaM (ADSL)
# Output: adqscibc.json
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
library(xportr)

# ----------------------------------------------------------------------------
# LOAD METADATA
# ----------------------------------------------------------------------------
adcibc_spec <- define_to_metacore(
  #path$define_path,
  file.path("~/Downloads", "define.xml"),
  quiet = TRUE
) %>%
  select_dataset("ADCIBC")

# ----------------------------------------------------------------------------
# LOAD DATASETS
# ----------------------------------------------------------------------------
# Define data to load
data_to_load_rds <- list(
  adsl = file.path(path$adam, "adsl.rds")
)

data_to_load_json <- list(
  qs = file.path(path$sdtm, "qs.json")
)

# Load datasets using map and convert blanks to NA
datasets_rds <- map(
  data_to_load_rds,
  ~ convert_blanks_to_na(readRDS(.x))
)

datasets_json <- map(
  data_to_load_json,
  ~ convert_blanks_to_na(read_dataset_json(.x))
)

# Put datasets into the global environment
list2env(c(datasets_rds, datasets_json), envir = .GlobalEnv)

# ----------------------------------------------------------------------------
# DERIVATIONS
# ----------------------------------------------------------------------------

# filter QS domain for qstestcd = CIBIC
adcibc00 <- qs %>%
  filter(QSTESTCD == "CIBIC") %>%
  select(STUDYID, USUBJID, VISIT, VISITNUM, QSDTC, QSSTRESN,
         QSSEQ) %>%
  # TODO: these are read as character from read_dataset_json
  mutate(VISITNUM = as.numeric(VISITNUM),
         QSSTRESN = as.numeric(QSSTRESN),
         QSSEQ = as.numeric(QSSEQ))

## ADSL information ----------------------------------------------
adsl_vars <- exprs(
  STUDYID,
  SITEID,
  SITEGR1,
  USUBJID,
  TRTSDT,
  TRTEDT,
  TRT01P,
  TRT01PN,
  AGE,
  AGEGR1,
  AGEGR1N,
  RACE,
  RACEN,
  SEX,
  ITTFL,
  #TODO: spec says this is called FASFL in ADSL??
  EFFFL,
  COMP24FL
)

adcibc1 <- adcibc00 %>%
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by = exprs(STUDYID, USUBJID)
  ) %>%
  rename(TRTP = TRT01P,
         TRTPN = TRT01PN)

# Derive dates -----------------------------------------------
# derive ADT and ADY
adcibc2 <- adcibc1 %>%
  derive_vars_dt(new_vars_prefix = "A",
                 dtc = QSDTC) %>%
  derive_vars_dy(reference_date = TRTSDT,
                 source_vars = exprs(ADT))

## Derive AVISIT, AVAL, PARAM, AVISITN, PARAMN -------------------
adcibc3 <- adcibc2 %>%
  mutate(
    AVISIT = case_when(
      ADY <= 1 ~ "Baseline",
      ADY >= 2 & ADY <= 84 ~ "Week 8",
      ADY >= 85 & ADY <= 140 ~ "Week 16",
      ADY > 140 ~ "Week 24",
      TRUE ~ NA_character_
    ),
    AVAL = QSSTRESN,
    PARAM = "CIBIC Score"
  ) %>%
  create_var_from_codelist(adcibc_spec, AVISIT, AVISITN) %>%
  create_var_from_codelist(adcibc_spec, PARAM, PARAMN) %>%
  create_var_from_codelist(adcibc_spec, PARAM, PARAMCD)

## Derive AWRANGE, AWTARGET, AWTDIFF, AWLO, AWHI, AWU ----------------
aw_lookup <- tribble(
  ~AVISIT, ~AWRANGE, ~AWTARGET, ~AWLO, ~AWHI,
  "Baseline", "<=1", 1, NA_integer_, 1,
  "Week 8", "2-84", 56, 2, 84,
  "Week 16", "85-140", 112, 85, 140,
  "Week 24", ">140", 168, 141, NA_integer_
)

adcibc4 <- derive_vars_merged(
  adcibc3,
  dataset_add = aw_lookup,
  by_vars = exprs(AVISIT)
) %>%
  mutate(
    AWTDIFF = abs(AWTARGET - ADY),
    AWU = "DAYS"
  )

## Derive ANL01FL ----------------------------------------
adcibc5 <- adcibc4 %>%
  mutate(diff = AWTARGET - ADY) %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(USUBJID, PARAMCD, AVISIT),
      order = exprs(AWTDIFF, diff),
      new_var = ANL01FL,
      mode = "first"
    ),
    filter = !is.na(AVISIT)
  )

## Derive DTYPE=LOCF -----------------------------------------
# A dataset with combinations of PARAMCD, AVISIT which are expected.
cibic_expected_obsv <- tibble::tribble(
  ~PARAMCD, ~AVISITN, ~AVISIT,
  #"CIBICVAL", 0, "Baseline",
  "CIBICVAL", 8, "Week 8",
  "CIBICVAL", 16, "Week 16",
  "CIBICVAL", 24, "Week 24"
)

adcibc_locf <- adcibc5 %>%
  restrict_derivation(
    derivation = derive_locf_records,
    args = params(
      dataset_ref = cibic_expected_obsv,
      by_vars = exprs(
        STUDYID, SITEID, SITEGR1, USUBJID, TRTSDT, TRTEDT,
        TRTP, TRTPN, AGE, AGEGR1, AGEGR1N, RACE, RACEN, SEX,
        ITTFL, EFFFL, COMP24FL, PARAMCD
      ),
      order = exprs(AVISITN, AVISIT),
      keep_vars = exprs(VISIT, VISITNUM, ADY, ADT, PARAM, PARAMN, QSSEQ)
    ),
    filter = TRUE #!is.na(ANL01FL)
  ) %>%
  # assign ANL01FL for new records
  mutate(ANL01FL = if_else(is.na(DTYPE), ANL01FL, "Y")) %>%
  # re-derive AWRANGE/AWTARGET/AWTDIFF/AWLO/AWHI/AWU
  select(-c("AWRANGE", "AWTARGET", "AWLO", "AWHI")) %>%
  derive_vars_merged(
    dataset_add = aw_lookup,
    by_vars = exprs(AVISIT)
  ) %>%
  mutate(
    AWTDIFF = abs(AWTARGET - ADY),
    AWU = "DAYS"
  ) %>%
  filter(!is.na(ADT))

adcibc <- adcibc_locf %>%
  mutate_if(is.numeric, as.integer) %>%
  drop_unspec_vars(adcibc_spec) %>%
  check_ct_data(adcibc_spec, na_acceptable = TRUE) %>%
  order_cols(adcibc_spec) %>%
  sort_by_key(adcibc_spec) %>%
  set_variable_labels(adcibc_spec) %>%
  xportr_df_label(adcibc_spec, domain = "ADAE") %>%
  xportr_label(adcibc_spec) %>%
  xportr_format(adcibc_spec$var_spec, "ADAE") %>%
  convert_na_to_blanks()

# Prepare column metadata
oid_cols <- adcibc_spec$ds_vars %>%
  select(dataset, variable, key_seq) %>%
  left_join(adcibc_spec$var_spec, by = c("variable")) %>%
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
dataset_json(adcibc,
             last_modified = strftime(as.POSIXlt(Sys.time(), "UTC"), "%Y-%m-%dT%H:%M"),
             originator = "R Submission Pilot 6",
             sys = paste0("R on ", R.Version()$os, " ", unname(Sys.info())[[2]]),
             sys_version = R.Version()$version.string,
             version = "1.1.0",
             study = "Pilot 6",
             metadata_version = "MDV.TDF_ADaM.ADaM-IG.1.1", # from define
             metadata_ref = file.path(path$adam, "define.xml"),
             item_oid = paste0("IG.ADCIBC"),
             name = "ADCIBC",
             dataset_label = adcibc_spec$ds_spec[["label"]],
             file_oid = file.path(path$adam, "adqscibc.json"),
             columns = oid_cols
) %>%
  write_dataset_json(file = file.path(path$adam, "adqscibc.json"), float_as_decimals = FALSE)


# Print summary
cat("\n============================================================================\n")
cat("ADCIBC Dataset Creation Complete\n")
cat("============================================================================\n")
cat("Output file:", file.path(path$adam, "adqscibc.json"), "\n")
cat("Number of records:", nrow(adcibc), "\n")
cat("Number of subjects:", length(unique(adcibc$USUBJID)), "\n")
cat("Number of parameters:", length(unique(adcibc$PARAMCD)), "\n")
cat("============================================================================\n")

# Temp: compare -----------------------
test <- read_dataset_json("~/Downloads/adqscibc.json")

diffdf(adcibc, test, keys = c("USUBJID", "PARAMCD", "AVISIT", "ADT"))

adcibc %>%
  filter(USUBJID == "01-705-1310") %>%
  select(USUBJID, AVISIT, AVISITN, ADT, AWTARGET, AWTDIFF, AWRANGE, DTYPE, ANL01FL, AVAL, QSSEQ) %>%
  arrange(AVISITN, ADT)

test %>%
  filter(USUBJID == "01-705-1310") %>%
  select(USUBJID, AVISIT, AVISITN, ADT, AWTARGET, AWTDIFF, AWRANGE, ADY, DTYPE, ANL01FL, AVAL, QSSEQ)

qs %>%
  filter(USUBJID == "01-705-1310",
         QSTESTCD == "CIBIC") %>%
  select(USUBJID, QSSEQ, QSTESTCD, QSSTRESN)