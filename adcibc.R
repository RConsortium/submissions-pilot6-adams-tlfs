# Setup -----------------
## Load libraries -------
library(dplyr)
library(tidyr)
library(admiral)
library(metacore)
library(metatools)
library(stringr)
library(purrr)
library(metacore)

path <- list(sdtm = "~/submissions-pilot5-datasetjson/pilot5-submission/pilot5-input/sdtmdata/",
             adam = "~/submissions-pilot5-datasetjson/pilot5-submission/pilot5-input/adamdata/")

## Load datasets ------------
dat_to_load <- list(
  qs = file.path(path$sdtm, "qs.rds"),
  adsl = file.path(path$adam, "adsl.rds")
)

datasets <- map(
  dat_to_load,
  ~ convert_blanks_to_na(readRDS(.x))
)

list2env(datasets, envir = .GlobalEnv)

## Load dataset specs -----------
adcibc_spec <- define_to_metacore(
  file.path("~/Downloads", "define.xml"),
  quiet = TRUE
) %>%
  select_dataset("ADCIBC")

# filter QS domain for qstestcd = CIBIC
adcibc00 <- qs %>%
  filter(QSTESTCD == "CIBIC") %>%
  select(STUDYID, USUBJID, VISIT, VISITNUM, QSDTC, QSSTRESN,
         QSSEQ)

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

# derive dates ------------
# derive ADT and ADY
adcibc2 <- adcibc1 %>%
  derive_vars_dt(new_vars_prefix = "A",
                 dtc = QSDTC) %>%
  derive_vars_dy(reference_date = TRTSDT,
                 source_vars = exprs(ADT))

## Derive AVISIT, AVAL, PARAM, AVISITN, PARAMN ----------
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

## Derive AWRANGE, AWTARGET, AWTDIFF, AWLO, AWHI, AWU -----------
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

## Derive ANL01FL -----------
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

## Derive DTYPE=LOCF -------------
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
  drop_unspec_vars(adcibc_spec) %>%
  check_ct_data(adcibc_spec, na_acceptable = TRUE) %>%
  order_cols(adcibc_spec) %>%
  sort_by_key(adcibc_spec) %>%
  set_variable_labels(adcibc_spec) %>%
  # xportr_df_label(adcibc_spec, domain = "ADAE") %>%
  # xportr_label(adcibc_spec) %>%
  # xportr_format(adae_spec$var_spec, "ADAE") %>%
  convert_na_to_blanks() %>%
  mutate_if(is.numeric, as.integer)

# compare
test <- read_dataset_json("~/Downloads/adqscibc.json")

adcibc5 %>% select(USUBJID, AVISIT, AVISITN) %>% slice(1:20)
test %>% select(USUBJID, AVISIT, AVISITN) %>% slice(1:20) %>% View()
test %>% filter(USUBJID == "01-701-1111") %>% select(USUBJID, AVISIT, AVISITN)
adcibc5 %>% filter(USUBJID == "01-701-1111") %>% select(USUBJID, AVISIT, AVISITN)

adcibc %>%
  filter(USUBJID == "01-701-1294") %>%
  select(USUBJID, AVISIT, AVISITN, ADT, AWTARGET, AWTDIFF, AWRANGE, DTYPE, ANL01FL, AVAL) %>%
  arrange(AVISITN, ADT)

test %>%
  filter(USUBJID == "01-701-1294") %>%
  select(USUBJID, AVISIT, AVISITN, ADT, AWTARGET, AWTDIFF, AWRANGE, DTYPE, ANL01FL, AVAL)

adcibc00 %>% filter(USUBJID == "01-701-1294")


diffdf(adcibc, test, keys = c("USUBJID", "PARAMCD", "AVISIT", "ADT"))

