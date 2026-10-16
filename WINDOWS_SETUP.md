# Running Pilot 6 on Windows with Positron

This guide gets you from a fresh Windows machine to running the R programs in this repo from Positron.
The SDTM/ADaM data lives in an AWS S3 bucket managed by DVC, so you need Python (for DVC) in addition to R.

## 1. Prerequisites

You need all three of these installed:

| Tool | Check | Get it |
|------|-------|--------|
| R | `R --version` in a terminal, or R starts in Positron | <https://cran.r-project.org/bin/windows/base/> |
| Positron | — | <https://positron.posit.co/download> |
| Python 3.10+ | `py --list` (see below) | <https://www.python.org/downloads/windows/> |

If you can't install software yourself, ask IT to install **Python 3 (64-bit) from python.org, with the
"py launcher" option and "Add python.exe to PATH" checked**. The Microsoft Store version also works, but the
python.org installer is the more reliable option on managed machines.

## 2. Check whether Python is already installed

Open a terminal in Positron (**Terminal > New Terminal**; it's PowerShell by default) and run:

```powershell
py --list
```

- **You see a version like `-V:3.13 *  Python 3.13 (64-bit)`:** Python is installed. Continue to step 3.
- **`py` is not recognized:** Python isn't installed. Install it, or ask IT to (see step 1).

> **Why `py` and not `python`?** On Windows, typing `python` often runs a Microsoft Store placeholder that just
> prints *"Python was not found"*, even when Python is installed. The `py` launcher finds the real installation
> either way, so this guide uses `py` throughout.

## 3. Clone the repo and open it in Positron

```powershell
git clone https://github.com/RConsortium/submissions-pilot6-adams-tlfs.git
```

In Positron, use **File > Open Folder…** and select the `submissions-pilot6-adams-tlfs` folder. Opening the
project root matters: `.Rprofile` builds the shared `path` list from the working directory.

## 4. Install DVC with S3 support

In the Positron terminal (from the project root):

```powershell
py -m pip install --upgrade pip
py -m pip install "dvc[s3]"
```

- Keep the quotes around `"dvc[s3]"` so PowerShell passes the brackets through unchanged.
- `py -m pip` runs pip through Python, so it works even when pip isn't on your PATH.

You may see a warning that scripts were installed to a `...\Python313\Scripts` folder that is **not on PATH**.
That's fine. Run DVC through Python instead:

```powershell
py -m dvc --version
```

<details>
<summary>Optional: make plain <code>dvc</code> work</summary>

1. Open the Start menu and search for **"Edit environment variables for your account"** (no admin needed).
2. Under *User variables*, edit **Path** and add the `Scripts` folder from the pip warning, for example:
   `C:\Users\<you>\AppData\Local\Programs\Python\Python313\Scripts`
3. Fully close and reopen Positron so it picks up the new PATH.
4. Check with `dvc --version`.

</details>

## 5. Pull the data

The DVC remote is already configured in the repo and allows anonymous read access, so you don't need AWS
credentials to pull:

```powershell
py -m dvc status
py -m dvc pull
```

This fills `data/sdtm/` and `data/adam/` with the Dataset-JSON files.

> To **push** data back to the bucket you do need credentials. See *Contributor Setup* in [README.md](README.md).

## 6. Run the R programs

Start an R console in Positron and run the scripts you need, for example:

- ADaM datasets: `code/adam/adsl.r`, `code/adam/adae.r`, …
- Tables: `code/tables/t_14_1_01.r`, …

Install any missing R packages that the scripts ask for (for example admiral, metacore, gt, gtsummary,
docorator, datasetjson).

## Troubleshooting

| Symptom | Fix |
|---------|-----|
| `Python was not found; run without arguments to install from the Microsoft Store` | You typed `python`. Use `py`. |
| `py` is not recognized | Python isn't installed, or was installed without the launcher. Reinstall from python.org with the **py launcher** option checked. |
| `dvc` is not recognized | Use `py -m dvc …`, or add the Scripts folder to PATH (step 4). |
| `ERROR: No matching distribution found for install` | You typed `install` twice: `pip install install …`. |
| `pip install` is blocked by a proxy or SSL error | Ask IT for your organization's pip proxy/index settings. |
| `dvc pull` fails with an S3 error | Run it from the project root, and check `py -m dvc remote list` shows `myremote`. |
