# submissions-pilot6-expansion

Repo to build out more ADaM and TLF R programs for future submissions.

Check out the [CDISC SAP for Datasets and TLFs](https://github.com/RConsortium/submissions-pilot6-adams-tlfs/blob/827ca9ae0025df78451ae6824c1fb8ada051c9ef/cdiscpilot01.pdf)

Welcome to the **R Consortium Pilot6-Expansion Project** of the Submissions Working Group!

Please feel welcome to try out the code on your device or even on posit cloud.
This pilot uses both a code repository (using git) and a separate data repository
(using dvc) with the stdm/adam datasets stored on AWS s3 bucket. In order to run
the R programs for the production datasets a little setup is required so your sandbox
has the proper coordinates and package support. 

# Sandbox Setup

*Pre-requisite Steps*

1. Sign-on or Sign-up to access the remote R Consortium dvc configuration for Pilot 6.

    https://dvc-callback-095d6d03.s3.amazonaws.com/index.html

2. Upon logging in the dvc repository you will receive the remote configuration with 
   the access/secret keys and session token. These are temporary and must be refreshed 
   periodically.
   
3. Once you can authenticate to access the dvc configuration, keys, and session token,
   please follow the 'Sandbox Steps' below to run any/all of the production adam datasets.

*Sandbox Steps*

1. Clone the Pilot6 repository to your device or workbench

   $git clone https://github.com/RConsortium/submissions-pilot6-adams-tlfs.git
   

2. After the repository is successfully cloned and upnpacked in your sandbox, change
   directory (cd) into the 'submissions-pilot6-adams-tlfs' folder.
   
   $cd pilot6-dvc/submissions-pilot6-adams-tlfs
   
   
3. Checkout a new git branch 'pilot6-dvc-test' to run your ADaM peoduction tests.

   $git checkout -b pilot6-dvc-test
   
4. Confirm installtion of dvc in your workspace (device, cloud). Check version of dvc.

   $dvc --version
   
   If you do not have an environment that supports dvc (i.e. dvc command not found result)
   then you will need to install dvc, in particular, the dvc[s3] module since Pilot 6 uses
   AWS s3 bucket storage.
   
   4(a). The Pilot 6 Team has been using posit Workbench so Python and the pip Python
         package manager are supported by default in the environment. 
         
         $pip --version
         
    4(b). If your environment does not support the Python package manager then please 
          reach out to your organization sysadmin to install Python and the Python 
          Package manager (pip).

 5. Install dvc with AWS s3 support module.
 
    $pip install dvc[s3]
    
 6. With dvc[s3] module installed be sure to configure your environment by running the 
    commands to provision/refresh the access/secret keys and session token. Note the
    dvc remote configuration comes out-of-the-box with your clone of the  Pilot 6 repo.
    
 7. Once the dev remote is configured and your session keys and token pull the sdtm and adam 
    data sets from the dvc repote.
    
    $dvc pull
    
  8. Run any/all of the production datasets using the R programs found under code/adam relative
     path from your project root. The R production code will run the dataset is output to the
     data/adam relative path from root directory in your project sanbox.
     
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
