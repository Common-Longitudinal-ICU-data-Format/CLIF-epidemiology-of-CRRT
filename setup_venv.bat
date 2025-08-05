@echo off
REM CLIF CRRT Epidemiology - Virtual Environment Setup Script for Windows
REM This script creates a Python virtual environment and installs required packages

setlocal enabledelayedexpansion

REM Check if Python is installed
python --version >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Python is not installed. Please install Python 3.8 or higher.
    exit /b 1
)

REM Get Python version
for /f "tokens=2" %%i in ('python --version 2^>^&1') do set PYTHON_VERSION=%%i
echo [INFO] Found Python version: %PYTHON_VERSION%

REM Check Python version is 3.8 or higher
python -c "import sys; exit(0 if sys.version_info >= (3, 8) else 1)"
if errorlevel 1 (
    echo [ERROR] Python 3.8 or higher is required. Current version: %PYTHON_VERSION%
    exit /b 1
)

REM Set virtual environment name
set VENV_NAME=venv

REM Check if virtual environment already exists
if exist "%VENV_NAME%" (
    echo [WARNING] Virtual environment '%VENV_NAME%' already exists.
    set /p REPLY="Do you want to delete and recreate it? (y/n): "
    if /i "!REPLY!"=="y" (
        echo [INFO] Removing existing virtual environment...
        rmdir /s /q "%VENV_NAME%"
    ) else (
        echo [INFO] Keeping existing virtual environment.
        echo [INFO] Activating virtual environment...
        call "%VENV_NAME%\Scripts\activate.bat"
        echo [INFO] Upgrading pip...
        python -m pip install --upgrade pip
        echo [INFO] Installing/updating requirements...
        pip install -r requirements.txt
        echo [INFO] Setup complete!
        goto :end
    )
)

REM Create virtual environment
echo [INFO] Creating virtual environment '%VENV_NAME%'...
python -m venv "%VENV_NAME%"

REM Check if virtual environment was created successfully
if not exist "%VENV_NAME%" (
    echo [ERROR] Failed to create virtual environment.
    exit /b 1
)

REM Activate virtual environment
echo [INFO] Activating virtual environment...
call "%VENV_NAME%\Scripts\activate.bat"

REM Upgrade pip
echo [INFO] Upgrading pip...
python -m pip install --upgrade pip

REM Install wheel for faster installations
echo [INFO] Installing wheel for faster package installations...
pip install wheel

REM Install requirements
if exist "requirements.txt" (
    echo [INFO] Installing requirements from requirements.txt...
    pip install -r requirements.txt
    
    if !errorlevel! equ 0 (
        echo [INFO] All requirements installed successfully!
    ) else (
        echo [ERROR] Some packages failed to install. Please check the error messages above.
        exit /b 1
    )
) else (
    echo [ERROR] requirements.txt not found in the current directory.
    exit /b 1
)

REM Install Jupyter kernel for the virtual environment
echo [INFO] Installing Jupyter kernel for this environment...
python -m ipykernel install --user --name="crrt-epidemiology" --display-name="CRRT Epidemiology (Python %PYTHON_VERSION%)"

REM Create directories if they don't exist
echo [INFO] Creating project directories...
if not exist "output\final\graphs" mkdir "output\final\graphs"
if not exist "output\intermediate" mkdir "output\intermediate"
if not exist "logs" mkdir "logs"

REM Summary
echo.
echo [INFO] ====================================
echo [INFO] Setup completed successfully!
echo [INFO] ====================================
echo.
echo To activate the virtual environment, run:
echo     %VENV_NAME%\Scripts\activate.bat
echo.
echo To deactivate the virtual environment, run:
echo     deactivate
echo.
echo To use this environment in Jupyter, select the kernel:
echo     'CRRT Epidemiology (Python %PYTHON_VERSION%)'

:end
endlocal