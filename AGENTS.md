# AGENTS.md

## Purpose
This repository builds ADaM analysis datasets from SDTM Dataset-JSON inputs using R, tidyverse workflows, and admiral/metacore tooling.

## Repository conventions discovered from current files

### Project layout
- `code/adam/`: one dataset-building script per ADaM domain (e.g., `adsl.r`, `adae.r`, `advs.r`).
- `code/utils/`: shared utilities (currently `save_dataset_json.r`).
- `data/sdtm/`: SDTM Dataset-JSON inputs (`*.json`, tracked with DVC).
- `data/adam/`: generated ADaM Dataset-JSON outputs (`*.json`, tracked with DVC).
- `data/adam_reference/`: specs and define metadata inputs.
- `.Rprofile`: global `path` and `study` lists used by scripts.

### Style and linting
- Lint config is in `.lintr`.
- Line length target is **120** characters.
- Pipe style is enforced as `%>%` (not base `|>`).

### Script structure
- Scripts commonly follow:
  1. header comments (purpose/input/output),
  2. library loading,
  3. data/spec loading,
  4. derivations in staged objects (`adsl00`, `adsl01`, etc.),
  5. export/save.
- Shared paths come from `.Rprofile` via the global `path` list.
- Metadata workflow typically uses `spec_to_metacore()` + `select_dataset()`.
- Finalization commonly uses:
  `drop_unspec_vars()` → `check_ct_data()` → `order_cols()` → `sort_by_key()` → labeling/formatting → `convert_na_to_blanks()`.

### Data I/O conventions
- Dataset-JSON is the canonical read/write format.
- Inputs are usually read with `read_dataset_json(..., decimals_as_floats = TRUE)`.
- Outputs are written via `save_dataset_json()` utility in `code/utils/save_dataset_json.r`.

## Tidyverse style expectations for future changes

- Use `<-` for assignment.
- Prefer `%>%` pipelines for multi-step transformations.
- Prefer `snake_case` for local object names and helper functions.
- Keep CDISC variable names in canonical uppercase where required by standards.
- Use `dplyr::mutate()`, `dplyr::select()`, `dplyr::left_join()`, etc. in clear staged steps (avoid large monolithic mutations).
- Prefer `case_when()` over deeply nested `ifelse()` logic.
- Keep one logical derivation block per object to preserve traceability.
- Keep library imports at top of script.

## Additional practical R conventions for this repo

- Preserve deterministic ordering before export (explicit `arrange()` or `sort_by_key()`).
- Avoid hidden side effects; scripts should only create/modify expected output datasets.
- Validate assumptions with explicit checks for key edge cases (missing dates, mixed treatments, absent joins).
- Reuse `code/utils` helpers when logic is shared across domains.
- Keep comments focused on ADaM derivation intent and reviewer traceability.

## Validation notes

- Existing lint expectations are defined in `.lintr`.
- In this execution environment, `Rscript` may be unavailable; run lint/tests/checks in an R-enabled environment when needed.
