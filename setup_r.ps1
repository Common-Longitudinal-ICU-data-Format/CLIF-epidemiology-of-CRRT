# setup_r.ps1 — one-time R provisioning for the CRRT pipeline on Windows.
#
# WHY THIS EXISTS
#   On managed Windows workstations, R packages auto-installed from CRAN land in a
#   user folder as freshly-downloaded binaries carrying Windows' "Mark of the Web"
#   (a Zone.Identifier stream). Windows Defender Application Control / Smart App
#   Control then BLOCKS loading those DLLs — e.g.:
#       unable to load shared object '...\utf8\libs\x64\utf8.dll':
#       LoadLibrary failure: An Application Control policy has blocked this file.
#   Installing the packages once and then stripping Mark-of-the-Web with Unblock-File
#   resolves it, so the automated pipeline (run_pipeline.bat) loads them cleanly.
#
# USAGE  (run ONCE, from the project root, before run_pipeline.bat):
#   powershell -ExecutionPolicy Bypass -File .\setup_r.ps1
#
# NOTES
#   - Requires R (Rscript) on PATH.
#   - No admin needed for the default user library. If your org's policy is
#     signature-based (nothing here helps), fall back to an IT allow-rule for the
#     R library path, or run the R stack on the HPC.

$ErrorActionPreference = "Stop"

# Package list mirrors code/05_PSM_IPTW_CRRT_dose.R + 05b.
$packages = @(
  "tidyverse","readr","arrow","gtsummary","cmprsk","survival","jsonlite",
  "MatchIt","WeightIt","broom","cobalt","EValue","SuperLearner","randomForest",
  "xgboost","gam","survminer","survey","mice"
)

if (-not (Get-Command Rscript -ErrorAction SilentlyContinue)) {
  Write-Error "Rscript not found on PATH. Install R (>= 4.5) and re-run."
}

Write-Host "== Installing R packages (first run can take several minutes) =="
$pkgArg = ($packages | ForEach-Object { "'$_'" }) -join ","
Rscript -e "options(repos=c(CRAN='https://cloud.r-project.org')); pk<-c($pkgArg); new<-pk[!(pk %in% rownames(installed.packages()))]; if(length(new)){cat('Installing:', paste(new, collapse=', '), '\n'); install.packages(new)} else {cat('All packages already installed.\n')}"

Write-Host ""
Write-Host "== Stripping Mark-of-the-Web from R libraries =="
$libs = @("$env:LOCALAPPDATA\R\win-library", "$env:R_LIBS_USER")
try {
  $sysLib = (Rscript -e "cat(.Library)").Trim()
  if ($sysLib) { $libs += $sysLib }
} catch {}

foreach ($lib in ($libs | Where-Object { $_ } | Select-Object -Unique)) {
  if (Test-Path $lib) {
    Write-Host "  Unblocking: $lib"
    Get-ChildItem $lib -Recurse -File -ErrorAction SilentlyContinue | Unblock-File -ErrorAction SilentlyContinue
  }
}

Write-Host ""
Write-Host "== Done. Now run the pipeline: =="
Write-Host "   .\run_pipeline.bat                 (full: descriptive + causal)"
Write-Host "   .\run_pipeline.bat --descriptive-only   (skip the R causal stack)"
