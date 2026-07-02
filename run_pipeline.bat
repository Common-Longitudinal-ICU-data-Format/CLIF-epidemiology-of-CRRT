@echo off
setlocal enabledelayedexpansion

set PYTHONIOENCODING=utf-8
set SCRIPT_DIR=%~dp0

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

:: Set up logging
if not exist "%SCRIPT_DIR%output\final" mkdir "%SCRIPT_DIR%output\final"
for /f "tokens=*" %%i in ('python -c "import json; print(json.load(open(r'%SCRIPT_DIR%config\config.json'))['site_name'])" 2^>nul') do set SITE_NAME=%%i
if "%SITE_NAME%"=="" set SITE_NAME=unknown
for /f "tokens=2 delims==" %%a in ('wmic os get localdatetime /value') do set dt=%%a
set TIMESTAMP=%dt:~0,8%_%dt:~8,6%
set LOG_FILE=%SCRIPT_DIR%output\final\%SITE_NAME%_pipeline_%TIMESTAMP%.log

echo === CRRT Epidemiology Pipeline === > "%LOG_FILE%"
echo   Started: %date% %time% >> "%LOG_FILE%"
echo   Site: %SITE_NAME% >> "%LOG_FILE%"
echo   Log: %LOG_FILE% >> "%LOG_FILE%"

echo === CRRT Epidemiology Pipeline ===
echo   Started: %date% %time%
echo   Site: %SITE_NAME%
if "%DESCRIPTIVE_ONLY%"=="true" (
    echo   Mode: DESCRIPTIVE-ONLY ^(cohort + Table 1 + descriptive epi + SMR + low-dose; no causal R stack^)
) else (
    echo   Mode: FULL ^(descriptive + causal data prep + R stack^)
)
echo   Log: %LOG_FILE%
echo.

:: Scripts use relative paths like ..\config\config.json, so run from code\
cd /d "%SCRIPT_DIR%code"

echo === Descriptive + SMR data prep (Python) ===
echo === Descriptive + SMR data prep (Python) === >> "%LOG_FILE%"
echo.

:: Common descriptive steps (both modes)
for %%S in (
    00_cohort.py
    01_create_wide_df.py
    02_construct_crrt_tableone.py
    03_crrt_epidemiology.py
    03b_crrt_epi_smr.py
    06_low_dose_characterization.py
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

if "%DESCRIPTIVE_ONLY%"=="true" goto :done

:: Causal data prep (full mode only)
echo --- Running 04_build_causal_df ---
echo --- Running 04_build_causal_df --- >> "%LOG_FILE%"
uv run python "04_build_causal_df.py" >> "%LOG_FILE%" 2>&1
if errorlevel 1 (
    echo ERROR: 04_build_causal_df failed. See log: %LOG_FILE%
    echo ERROR: 04_build_causal_df failed. >> "%LOG_FILE%"
    exit /b 1
)
echo --- 04_build_causal_df complete ---
echo.

echo === Causal inference (R) ===
echo === Causal inference (R) === >> "%LOG_FILE%"
echo.

set R_LIBS_USER=%USERPROFILE%\R\win-library\4.5
mkdir "%R_LIBS_USER%" 2>nul

set R_FAILED=
for %%S in (
    05_PSM_IPTW_CRRT_dose.R
    05b_dose_response_analysis.R
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

:done
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
    echo     05_PSM_IPTW_CRRT_dose.R        -^> output\final\psm_iptw\
    echo     05b_dose_response_analysis.R    -^> output\final\psm_iptw\ (dose-response)
    echo.
    echo   Check the log for error details: %LOG_FILE%
)
