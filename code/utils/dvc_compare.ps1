#!/usr/bin/env pwsh

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$vdePath = 'vde-dataset-viewer'
$rPath = 'Rscript'
$textMode = $false
$rev = 'HEAD'
$file = $null

function Show-Usage {
    Write-Error "Usage: $($MyInvocation.MyCommand.Name) [--rev <rev>] [--vde-path <path/to/vde>] [--r-path <path/to/Rscript>] [--text] <path/to/file>"
    Write-Error '  --rev <rev>           Git revision to compare against (default: HEAD)'
    Write-Error '  --vde-path <path>     Path to vde-dataset-viewer executable (default: vde-dataset-viewer in PATH)'
    Write-Error '  --r-path <path>       Path to Rscript executable (default: Rscript in PATH)'
    Write-Error '  --text, -t            Use text-based comparison mode'
    Write-Error '  --visual, -v          Use visual comparison mode (default)'
    Write-Error '  --help, -h            Show this help message and exit'
}

function Fail-WithUsage {
    param(
        [int] $ExitCode
    )

    Show-Usage
    exit $ExitCode
}

$index = 0
while ($index -lt $args.Count) {
    $arg = [string] $args[$index]

    switch -Regex ($arg) {
        '^--rev$' {
            $index++
            if ($index -ge $args.Count) {
                Fail-WithUsage -ExitCode 2
            }
            $rev = [string] $args[$index]
        }
        '^--rev=(.+)$' {
            $rev = $Matches[1]
        }
        '^--vde-path$' {
            $index++
            if ($index -ge $args.Count) {
                Fail-WithUsage -ExitCode 2
            }
            $vdePath = [string] $args[$index]
        }
        '^--vde-path=(.+)$' {
            $vdePath = $Matches[1]
        }
        '^--r-path$' {
            $index++
            if ($index -ge $args.Count) {
                Fail-WithUsage -ExitCode 2
            }
            $rPath = [string] $args[$index]
        }
        '^--r-path=(.+)$' {
            $rPath = $Matches[1]
        }
        '^(--text|-t)$' {
            $textMode = $true
        }
        '^(--visual|-v)$' {
            $textMode = $false
        }
        '^(--help|-h)$' {
            Show-Usage
            exit 0
        }
        '^--' {
            Write-Error "Unknown option: $arg"
            Fail-WithUsage -ExitCode 2
        }
        default {
            $file = $arg
        }
    }

    $index++
}

if (-not $file) {
    Fail-WithUsage -ExitCode 2
}

if (-not (Test-Path -LiteralPath "$file.dvc" -PathType Leaf)) {
    Write-Error "not in dvc: $file"
    exit 1
}

$resolvedFile = Resolve-Path -LiteralPath $file
$filePath = $resolvedFile.Path
$directory = Split-Path -Path $filePath -Parent
$baseName = Split-Path -Path $filePath -Leaf
$tmp = Join-Path -Path $directory -ChildPath ("~compare.{0}" -f $baseName)
$scriptDir = Split-Path -Path $PSCommandPath -Parent

try {
    $gitRootResult = & git rev-parse --show-toplevel 2>&1
    if ($LASTEXITCODE -ne 0) {
        throw "git rev-parse failed: $($gitRootResult -join [Environment]::NewLine)"
    }

    $gitRoot = (($gitRootResult | Select-Object -First 1).ToString()).Trim()
    $normalizedRoot = [System.IO.Path]::GetFullPath($gitRoot).TrimEnd('\', '/')
    $normalizedFile = [System.IO.Path]::GetFullPath($filePath)
    $rootPrefix = $normalizedRoot + [System.IO.Path]::DirectorySeparatorChar

    if (-not $normalizedFile.StartsWith($rootPrefix, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "File is not inside git root: $filePath"
    }

    $fileRelativePath = $normalizedFile.Substring($rootPrefix.Length) -replace '\\', '/'

    $statusOutput = & dvc status $filePath --json 2>&1
    $statusJson = ($statusOutput | ForEach-Object { $_.ToString() }) -join [Environment]::NewLine
    if ($LASTEXITCODE -ne 0) {
        Write-Error "dvc status failed for $fileRelativePath: $statusJson"
        exit 3
    }

    if ($statusJson.Trim() -eq '{}') {
        Write-Output "File not modified relative to committed DVC version: $file"
        if ($rev -eq 'HEAD') {
            exit 0
        }
    }
    else {
        Write-Output "Workspace file modified; fetching committed DVC version to: $tmp"
    }

    & dvc get . $fileRelativePath -o $tmp --rev $rev *> $null
    if ($LASTEXITCODE -ne 0) {
        Write-Error "Failed to retrieve committed DVC version for $file (rev=$rev)"
        exit 5
    }

    if ($textMode) {
        if (-not (Get-Command $rPath -ErrorAction SilentlyContinue)) {
            Write-Error 'Rscript not found in PATH'
            exit 6
        }

        Write-Output ("Running {0} {1} {2} {3}" -f $rPath, (Join-Path $scriptDir 'compare_data.r'), $filePath, $tmp)
        & $rPath (Join-Path $scriptDir 'compare_data.r') $filePath $tmp
        exit $LASTEXITCODE
    }

    Write-Output "Launching $vdePath --compare $filePath $tmp"
    & $vdePath --compare $filePath $tmp *> $null
    exit $LASTEXITCODE
}
finally {
    if (Test-Path -LiteralPath $tmp -PathType Leaf) {
        Remove-Item -LiteralPath $tmp -Force
    }
}