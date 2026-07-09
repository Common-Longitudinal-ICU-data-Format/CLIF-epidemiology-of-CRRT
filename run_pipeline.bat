@echo off
setlocal enabledelayedexpansion

set PYTHONIOENCODING=utf-8
set SCRIPT_DIR=%~dp0

:: ── Auto-log (parity with run_pipeline.sh): on first entry, re-run through
::    PowerShell Tee-Object so output goes to BOTH the console AND a timestamped
::    log file (wmic-free). The "__teed" sentinel prevents infinite recursion. ──
if /i not "%~1"=="__teed" (
    if not exist "%SCRIPT_DIR%output\final_no_phi" mkdir "%SCRIPT_DIR%output\final_no_phi"
    set "SITE_NAME=unknown"
    for /f "usebackq delims=" %%i in (`uv run python -c "import json; print(json.load(open(r'%SCRIPT_DIR%config\config.json'))['site_name'])" 2^>nul`) do set "SITE_NAME=%%i"
    for /f %%t in ('powershell -NoProfile -Command "Get-Date -Format yyyyMMdd_HHmmss"') do set "TS=%%t"
    set "LOG=%SCRIPT_DIR%output\final_no_phi\!SITE_NAME!_pipeline_!TS!.log"
    echo === Logging to: !LOG! ===
    powershell -NoProfile -Command "& '%~f0' __teed %* 2>&1 | Tee-Object -FilePath '!LOG!'"
    exit /b !errorlevel!
)
shift

:: ── Mode parsing: --descriptive-only runs the descriptive + SMR deliverable
::    (00 01 02 03 03b 06) and skips the causal stack (04 + the R scripts). ──
set DESCRIPTIVE_ONLY=false
if /i "%~1"=="--descriptive-only" set DESCRIPTIVE_ONLY=true
if /i "%~1"=="-h"      ( echo Usage: run_pipeline.bat [--descriptive-only] & exit /b 0 )
if /i "%~1"=="--help"  ( echo Usage: run_pipeline.bat [--descriptive-only] & exit /b 0 )
if not "%~1"=="" if /i not "%~1"=="--descriptive-only" (
    echo Unknown argument: %~1
    echo Usage: run_pipeline.bat [--descriptive-only]
    exit /b 1
)

:: Check config exists
if not exist "%SCRIPT_DIR%config\config.json" (
    echo ERROR: config\config.json not found.
    echo Copy the template and edit it:
    echo   copy config\config_template.json config\config.json
    exit /b 1
)

:: Site name (via the project's uv env, not bare python)
set SITE_NAME=unknown
for /f "usebackq delims=" %%i in (`uv run python -c "import json; print(json.load(open(r'%SCRIPT_DIR%config\config.json'))['site_name'])" 2^>nul`) do set "SITE_NAME=%%i"

echo === CRRT Epidemiology Pipeline ===
echo   Started: %date% %time%
echo   Site: %SITE_NAME%
if "%DESCRIPTIVE_ONLY%"=="true" (
    echo   Mode: DESCRIPTIVE-ONLY ^(cohort + Table 1 + descriptive epi + SMR + low-dose; no causal R stack^)
) else (
    echo   Mode: FULL ^(descriptive + causal data prep + R stack^)
)
echo.

:: Scripts use relative paths like ..\config\config.json, so run from code\
cd /d "%SCRIPT_DIR%code"

echo === Descriptive + SMR data prep (Python) ===
echo.

:: Common descriptive steps (both modes) — output streams live to the console
for %%S in (
    00_cohort.py
    01_create_wide_df.py
    02_construct_crrt_tableone.py
    03_crrt_epidemiology.py
    03b_crrt_epi_smr.py
    06_low_dose_characterization.py
) do (
    echo --- Running %%~nS ---
    uv run python "%%S"
    if errorlevel 1 (
        echo ERROR: %%~nS failed.
        exit /b 1
    )
    echo --- %%~nS complete ---
    echo.
)

if "%DESCRIPTIVE_ONLY%"=="true" goto :done

:: Causal data prep (full mode only)
echo --- Running 04_build_causal_df ---
uv run python "04_build_causal_df.py"
if errorlevel 1 (
    echo ERROR: 04_build_causal_df failed.
    exit /b 1
)
echo --- 04_build_causal_df complete ---
echo.

echo === Causal inference (R) ===
echo.

set R_FAILED=
where Rscript >nul 2>nul
if errorlevel 1 (
    echo   Rscript not found on PATH -- the R ^(causal^) steps did NOT run.
    echo   This is expected when Python ^(uv^) and R live in different environments.
    echo   Run the causal step manually in RStudio ^(see summary below^).
    set "R_FAILED=05_PSM_IPTW_CRRT_dose.R 05b_dose_response_analysis.R"
    goto :after_r
)
:: Run R from the project root (NOT code\) so the project .Rprofile's
:: source("renv/activate.R") resolves and the pinned renv library is used.
:: Do NOT pass --no-init-file (it would skip .Rprofile and thus renv).
cd /d "%SCRIPT_DIR%"
echo --- Ensuring pinned R packages (renv::restore) ---
Rscript -e "renv::restore(prompt = FALSE)"
echo.
for %%S in (
    05_PSM_IPTW_CRRT_dose.R
    05b_dose_response_analysis.R
) do (
    echo --- Running %%~nS ---
    Rscript "%SCRIPT_DIR%code\%%S"
    if errorlevel 1 (
        echo --- %%~nS FAILED ---
        set "R_FAILED=!R_FAILED! %%S"
    ) else (
        echo --- %%~nS complete ---
    )
    echo.
)
:after_r

:done
echo === Pipeline complete ===
echo   Finished: %date% %time%
echo   Results: output\final_no_phi\

if defined R_FAILED (
    echo.
    echo === The R ^(causal^) steps did not complete ===
    echo   Run them manually in RStudio ^(the scripts self-locate config.json^):
    echo     1. Open this project folder in RStudio, then run:  renv::restore^(^)
    echo     2. Open  code\05_PSM_IPTW_CRRT_dose.R  and click Source.
    echo     3. Open  code\05b_dose_response_analysis.R  and click Source.
    echo.
    echo   R outputs: output\final_no_phi\psm_iptw\
)
