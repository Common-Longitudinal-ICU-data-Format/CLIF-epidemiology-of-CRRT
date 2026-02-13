@echo off
setlocal

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

for %%S in (
    00_cohort.py
    01_create_wide_df.py
    02_construct_crrt_tableone.py
    03_crrt_visualizations.py
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

echo === Pipeline complete ===
echo Results are in output\final\
