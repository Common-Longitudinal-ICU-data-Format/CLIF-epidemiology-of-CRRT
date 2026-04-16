@echo off
setlocal enabledelayedexpansion

set PYTHONIOENCODING=utf-8
set SCRIPT_DIR=%~dp0

:: Check config exists
if not exist "%SCRIPT_DIR%config\config.json" (
    echo ERROR: config\config.json not found.
    echo Copy the template and edit it:
    echo   copy config\config_template.json config\config.json
    exit /b 1
)

:: Set up logging
if not exist "%SCRIPT_DIR%output\final" mkdir "%SCRIPT_DIR%output\final"
for /f "tokens=*" %%i in ('python -c "import json; print(json.load(open(r'%SCRIPT_DIR%config\config.json'))['site_name'])" 2^>nul') do set SITE_NAME=%%i
if "%SITE_NAME%"=="" set SITE_NAME=unknown
for /f "tokens=2 delims==" %%a in ('wmic os get localdatetime /value') do set dt=%%a
set TIMESTAMP=%dt:~0,8%_%dt:~8,6%
set LOG_FILE=%SCRIPT_DIR%output\final\%SITE_NAME%_pipeline_%TIMESTAMP%.log

:: Tee all output to log file + console via PowerShell
echo === CRRT Epidemiology Pipeline === > "%LOG_FILE%"
echo   Started: %date% %time% >> "%LOG_FILE%"
echo   Site: %SITE_NAME% >> "%LOG_FILE%"
echo   Log: %LOG_FILE% >> "%LOG_FILE%"

echo === CRRT Epidemiology Pipeline ===
echo   Started: %date% %time%
echo   Site: %SITE_NAME%
echo   Log: %LOG_FILE%
echo.

:: Scripts use relative paths like ..\config\config.json, so run from code\
cd /d "%SCRIPT_DIR%code"

echo === Descriptive + MSM data prep (Python) ===
echo === Descriptive + MSM data prep (Python) === >> "%LOG_FILE%"
echo.
for %%S in (
    00_cohort.py
    01_create_wide_df.py
    02_construct_crrt_tableone.py
    03_crrt_visualizations.py
    04_build_msm_competing_risk_df.py
) do (
    echo --- Running %%~nS ---
    echo --- Running %%~nS --- >> "%LOG_FILE%"
    uv run python "%%S" >> "%LOG_FILE%" 2>&1
    if errorlevel 1 (
        echo ERROR: %%~nS failed. See log: %LOG_FILE%
        echo ERROR: %%~nS failed. >> "%LOG_FILE%"
        exit /b 1
    )
    echo --- %%~nS complete ---
    echo --- %%~nS complete --- >> "%LOG_FILE%"
    echo.
)

echo === Causal inference (R) ===
echo === Causal inference (R) === >> "%LOG_FILE%"
echo.

set R_LIBS_USER=%USERPROFILE%\R\win-library\4.5
mkdir "%R_LIBS_USER%" 2>nul

set R_FAILED=
for %%S in (
    05_PSM_IPTW_CRRT_dose.R
    05b_dose_response_analysis.R
    06_time_varying_MSM.R
    06b_time_varying_MSM_sensitivity.R
) do (
    echo --- Running %%~nS ---
    echo --- Running %%~nS --- >> "%LOG_FILE%"
    Rscript --no-init-file "%SCRIPT_DIR%code\%%S" >> "%LOG_FILE%" 2>&1
    if errorlevel 1 (
        echo --- %%~nS FAILED ---
        echo --- %%~nS FAILED --- >> "%LOG_FILE%"
        set "R_FAILED=!R_FAILED! %%S"
    ) else (
        echo --- %%~nS complete ---
        echo --- %%~nS complete --- >> "%LOG_FILE%"
    )
    echo.
)

echo === Pipeline complete ===
echo   Finished: %date% %time%
echo   Results: output\final\
echo   Log: %LOG_FILE%
echo === Pipeline complete === >> "%LOG_FILE%"
echo   Finished: %date% %time% >> "%LOG_FILE%"

if defined R_FAILED (
    echo.
    echo === WARNING: Some R scripts failed ===
    echo   Failed scripts:%R_FAILED%
    echo.
    echo   To re-run manually from the project root directory:
    echo.
    for %%F in (%R_FAILED%) do (
        echo     Rscript --no-init-file code\%%F
    )
    echo.
    echo   Required R scripts and their outputs:
    echo     05_PSM_IPTW_CRRT_dose.R            -^> output\final\psm_iptw\
    echo     05b_dose_response_analysis.R        -^> output\final\psm_iptw\ (dose-response)
    echo     06_time_varying_MSM.R               -^> output\final\time_varying\ (primary 12h)
    echo     06b_time_varying_MSM_sensitivity.R  -^> output\final\time_varying_sensitivity\ (24h)
    echo.
    echo   Check the log for error details: %LOG_FILE%
)
