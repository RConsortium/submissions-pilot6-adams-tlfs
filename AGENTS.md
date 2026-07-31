# AGENTS.md

## Purpose
This repository builds ADaM analysis datasets from SDTM Dataset-JSON inputs using R, tidyverse workflows, and admiral/metacore tooling.

## Repository conventions discovered from current files

### Project layout
- `code/adam/`: one dataset-building script per ADaM domain (e.g., `adsl.r`, `adae.r`, `advs.r`).
- `code/utils/`: shared utilities (currently `save_dataset_json.r`).
- `code/tables/`: one dataset-building script per ADaM domain (e.g., `adsl.r`, `adae.r`, `advs.r`).
- `data/sdtm/`: SDTM Dataset-JSON inputs (`*.json`, tracked with DVC).
- `data/adam/`: generated ADaM Dataset-JSON outputs (`*.json`, tracked with DVC).
- `data/adam_reference/`: specs and define metadata inputs.
- `.Rprofile`: global `path` and `study` lists used by scripts.

### Style and linting
- Lint config is in `.lintr`.
- Line length target is **120** characters.
- Pipe style is enforced as `%>%` (not base `|>`).
- If a migration to base `|>` is desired, update `.lintr` and affected scripts together in the same change.

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

### TLF output conventions (docorator)
- For table/listing/figure outputs, use `gt` + `docorator` for PDF rendering.
- Use `gtsummary` for calculating table summaries (descriptive statistics, frequencies, etc.).
- TLF program naming convention:
  - `t_<name>.r` for tables,
  - `l_<name>.r` for listings,
  - `f_<name>.r` for figures.
  - Use zero-padding in the last numeric segment of the name (e.g., `t_14_6_02.r`, not `t_14_6_2.r`).
- Output artifact naming should mirror program naming (for example, `t_<name>.rds` and `t_<name>.pdf`).
- Keep TLF programs under `code/tables/` and produce both:
  - ARD-style RDS output (for reproducible intermediate data),
  - PDF output (final review artifact).
- If `code/tables/` is absent in your branch, create it manually before adding TLF programs.
- Use configured path entries (for example `path$table_ard` and `path$table_output`) for output locations.
- Suggested table program flow:
  1. header comments (purpose/input/output),
  2. library loading,
  3. load ADaM source dataset(s) (typically from `data/adam/`),
  4. staged derivations for ARD/table-ready data,
  5. formatting/labels/footnotes,
  6. render with `gt` + `docorator`,
  7. save ARD `.rds` and final `.pdf`.
- Before rendering, validate input assumptions (required variables, join keys, and unexpected missingness patterns) so issues are flagged early.
- Keep header/footer/footnote formatting in the table script so reviewer-facing output remains traceable and reproducible.
- Use define metadata (`data/adam_reference/define.xml`) and ADaM variable lineage to keep table columns traceable to source domains.
- DVC workflow for TLF work:
  - run `dvc pull` before table generation when SDTM/ADaM inputs may be missing or stale,
  - run `dvc add` for generated table artifacts that are intended to be tracked,
  - commit corresponding `.dvc` files with program changes.

## Tidyverse style expectations for future changes

- Use `<-` for assignment.
- Prefer `|>` pipelines for multi-step transformations.
- Prefer `snake_case` for local object names and helper functions.
- Use lowercase for variable names in programming scripts (local variables, variables from source datasets).
- Keep CDISC variable names in canonical uppercase where required by standards (variables in the final SDTM/ADaM datasets).
- Use `dplyr::mutate()`, `dplyr::select()`, `dplyr::left_join()`, etc. in clear staged steps (avoid large monolithic mutations).
- Prefer `case_when()` over deeply nested `ifelse()` logic.
- Keep one logical derivation block per object to preserve traceability.
- Keep library imports at top of script.

## Additional practical R conventions for this repo

- Preserve deterministic ordering before export (explicit `arrange()` or `sort_by_key()`).
- Avoid hidden side effects; scripts should only create/modify expected output datasets.
- Validate assumptions with explicit checks for key edge cases (missing dates, mixed treatments, absent joins).
- Reuse `code/utils` helpers when logic is shared across domains.
- When adding new helper functions, include a short high-level overview of created functions (purpose and expected inputs/outputs) in script documentation or PR summary.
- Before introducing a new helper, scan existing scripts/utilities for overlapping logic that can be re-used or consolidated.
- If overlap is identified, propose a reusable utility in `code/utils/` and update call sites when in scope; otherwise, explicitly alert developers with the consolidation opportunity so they can decide.
- Keep comments focused on ADaM derivation intent and reviewer traceability.

## Validation notes

- Existing lint expectations are defined in `.lintr`.
- In this execution environment, `Rscript` may be unavailable; run lint/tests/checks in an R-enabled environment when needed.
