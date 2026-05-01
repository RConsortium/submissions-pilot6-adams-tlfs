# submissions-pilot6-expansion
Repo to build out more ADaM and TLF programs for future submissions

## DVC Comparison

Add the following alias to your local git configuration to enable easy comparison of DVC-tracked files:

```
git config --local alias.compare '!cd "${GIT_PREFIX:-.}" && "$(git rev-parse --show-toplevel)/code/utils/dvc_compare.sh"'
```

You can modify the alias to add default options. For example, to use the text comparison mode by default, update the alias as follows:
```
git config --local alias.compare '!cd "${GIT_PREFIX:-.}" && "$(git rev-parse --show-toplevel)/code/utils/dvc_compare.sh" --text'
```

This allows you to run `git compare <file>` to compare the specified file against the current HEAD revision.

Compare can be done in text and visual modes. For the visual mode you will need to have the [vde-dataset-viewer](https://github.com/defineEditor/vde-dataset-viewer). For the text mode you will need Rscript in your PATH and both `diffdf` and `datasetjson` packages installed.

Usage:

```
git compare [options] <path/to/file>
```

Options:

`--rev <rev>`: Git revision to compare against (default: HEAD)
`--vde-path <path>`: Path to vde-dataset-viewer executable (default: vde-dataset-viewer in PATH)
`--r-path <path>`: Path to Rscript executable (default: Rscript in PATH)
`--text, -t`: Use text-based comparison mode
`--visual, -v`: Use visual comparison mode (default)