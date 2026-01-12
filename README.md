# submissions-pilot6-expansion
Repo to build out more ADaM and TLF programs for future submissions


## DVC (Data Versioning Control) Project Setup Instructions
1. Install DVC from these instructions https://doc.dvc.org/start, also install dvc-s3 with `pip install dvc-s3`
2. Initialize project with `dvc init`
3. Track 'data/' directory with `dvc add data`
4. Push up all data with `dvc push`

## DVC Update instructions
1. Install DVC from these instructions https://doc.dvc.org/start, also install dvc-s3 with `pip install dvc-s3`
2. Pull current data down with `dvc pull`
3. Push new data up with `dvc update` and `dvc push`
