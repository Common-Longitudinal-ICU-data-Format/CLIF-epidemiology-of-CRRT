@echo off
SETLOCAL ENABLEDELAYEDEXPANSION

REM ── Step 0: Go to script directory ──
cd /d %~dp0

:create_venv
REM ── Step 2: Create virtual environment if missing ──
if not exist ".crrt\" (
    echo Creating virtual environment...
    python -m venv .crrt
) else (
    echo Virtual environment already exists.
)

REM ── Step 3: Activate virtual environment ──
call .crrt\Scripts\activate.bat

REM ── Step 4: Install required packages ──
echo Installing dependencies...
pip install --quiet -r requirements.txt
pip install --quiet jupyter ipykernel papermill

REM ── Step 5: Register kernel ──
python -m ipykernel install --user --name=.crrt --display-name="Python (CRRT)"

REM ── Step 6: Set environment variables ──
set PYTHONWARNINGS=ignore
set PYTHONPATH=%cd%\code;%PYTHONPATH%

REM ── Step 7: Change to code directory ──
cd code

REM ── Step 8: Create logs folder ──
if not exist logs (
    mkdir logs
)

REM ── Step 9: Run analysis notebooks using papermill ──
echo.
echo Running 01_cohort_identification.ipynb ...
papermill 01_cohort_identification.ipynb 01_cohort_identification.ipynb > logs\01_cohort_identification.log
echo Finished 01_cohort_identification.ipynb

echo.
echo Running 02_analysis_summary.ipynb ...
papermill 02_analysis_summary.ipynb 02_analysis_summary.ipynb > logs\02_analysis_summary.log
echo Finished 02_analysis_summary.ipynb

REM ── Step 10: Done ──
echo.
echo ✅ All CRRT epidemiology analysis steps completed successfully!
pause

REM ── Step 11: Ask to launch dashboard (if available) ──
echo.
echo Would you like to launch the visualization dashboard?
echo The dashboard provides interactive CRRT epidemiology analysis.
echo.

exit /b 0