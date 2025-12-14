library(dplyr)
library(datasetjson)  # Dataset JSON handling
library(metacore)     # Metadata handling

save_dataset_json <- function(output_dir, dataset, ds_spec) {

  # Check parameters
  stopifnot(is.character(output_dir), length(output_dir) == 1, nzchar(output_dir))
  stopifnot(is.data.frame(dataset))
  stopifnot(
    !is.null(ds_spec$ds_vars),
    !is.null(ds_spec$var_spec),
    !is.null(ds_spec$ds_spec)
  )


  ds_dataset <- ds_spec$ds_spec[["dataset"]]
  ds_name <- as.character(ds_dataset)[1]

  # Select variables as per spec
  dataset_final <- dataset %>%
    dplyr::select(
      ds_spec$ds_vars$variable
    )

  oid_cols <- ds_spec$ds_vars %>%
    select(.data$dataset, .data$variable, .data$key_seq) %>%
    left_join(ds_spec$var_spec, by = c("variable")) %>%
    rename(name = .data$variable, dataType = .data$type, keySequence = .data$key_seq, displayFormat = .data$format) %>%
    mutate(itemOID = paste0("IT.", .data$dataset, ".", .data$name)) %>%
    select(.data$itemOID, .data$name, .data$label, .data$dataType,
           .data$length, .data$keySequence, .data$displayFormat) %>%
    mutate(
      dataType =
        case_when(
          displayFormat == "DATE9." ~ "date",
          displayFormat == "DATETIME20." ~ "datetime",
          substr(name, nchar(name) - 3 + 1, nchar(name)) == "DTC" & length == "8" ~ "date",
          substr(name, nchar(name) - 3 + 1, nchar(name)) == "DTC" & length == "20" ~ "datetime",
          dataType == "text" ~ "string",
          .default = as.character(.data$dataType)
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

  dataset_json(dataset_final,
    last_modified = strftime(as.POSIXlt(Sys.time(), "UTC"), "%Y-%m-%dT%H:%M"),
    originator = study$originator,
    sys = paste0("R on ", R.Version()$os, " ", unname(Sys.info())[[2]]),
    sys_version = R.Version()$version.string,
    version = "1.1.0",
    study = study$name,
    metadata_version = "MDV.TDF_ADaM.ADaM-IG.1.1", # from define
    metadata_ref = file.path(path$adam, "define.xml"),
    item_oid = paste0("IG.", toupper(ds_name)),
    name = toupper(ds_name),
    dataset_label = ds_spec$ds_spec[["label"]],
    file_oid = file.path(path$adam, ds_name),
    columns = oid_cols
  ) %>%
    write_dataset_json(file = file.path(output_dir, paste0(tolower(ds_name), ".json")))
}
