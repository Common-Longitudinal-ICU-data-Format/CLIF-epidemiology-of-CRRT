@echo off
setlocal

set PYTHONIOENCODING=utf-8
set SCRIPT_DIR=%~dp0

:: Check config exists
if not exist "%SCRIPT_DIR%config\config.json" (
    echo ERROR: config\config.json not found.
    echo Copy the template and edit it:
    echo   copy config\config_template.json config\config.json
    exit /b 1
)

:: Scripts use relative paths like ..\config\config.json, so run from code\
cd /d "%SCRIPT_DIR%code"

echo === CRRT Epidemiology Pipeline ===
echo.

echo === Descriptive + MSM data prep (Python) ===
echo.
for %%S in (
    00_cohort.py
    01_create_wide_df.py
    02_construct_crrt_tableone.py
    03_crrt_visualizations.py
    04_build_msm_competing_risk_df.py
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

echo === Causal inference (R) ===
echo.
for %%S in (
    05_PSM_IPTW_CRRT_dose.R
    05b_dose_response_analysis.R
    06_time_varying_MSM.R
    06b_time_varying_MSM_sensitivity.R
) do (
    echo --- Running %%~nS ---
    Rscript "%%S"
    if errorlevel 1 (
        echo ERROR: %%~nS failed.
        exit /b 1
    )
    echo --- %%~nS complete ---
    echo.
)

echo === Pipeline complete ===
echo Results are in output\final\
