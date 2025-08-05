#!/bin/bash

# CLIF CRRT Epidemiology - Virtual Environment Setup Script
# This script creates a Python virtual environment and installs required packages

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Check if Python 3 is installed
if ! command -v python3 &> /dev/null; then
    print_error "Python 3 is not installed. Please install Python 3.8 or higher."
    exit 1
fi

# Get Python version
PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}')
print_status "Found Python version: $PYTHON_VERSION"

# Check Python version is 3.8 or higher
REQUIRED_VERSION="3.8"
if ! python3 -c "import sys; exit(0 if sys.version_info >= (3, 8) else 1)"; then
    print_error "Python 3.8 or higher is required. Current version: $PYTHON_VERSION"
    exit 1
fi

# Set virtual environment name
VENV_NAME="venv"

# Check if virtual environment already exists
if [ -d "$VENV_NAME" ]; then
    print_warning "Virtual environment '$VENV_NAME' already exists."
    read -p "Do you want to delete and recreate it? (y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        print_status "Removing existing virtual environment..."
        rm -rf "$VENV_NAME"
    else
        print_status "Keeping existing virtual environment."
        print_status "Activating virtual environment..."
        source "$VENV_NAME/bin/activate"
        print_status "Upgrading pip..."
        pip install --upgrade pip
        print_status "Installing/updating requirements..."
        pip install -r requirements.txt
        print_status "Setup complete!"
        exit 0
    fi
fi

# Create virtual environment
print_status "Creating virtual environment '$VENV_NAME'..."
python3 -m venv "$VENV_NAME"

# Check if virtual environment was created successfully
if [ ! -d "$VENV_NAME" ]; then
    print_error "Failed to create virtual environment."
    exit 1
fi

# Activate virtual environment
print_status "Activating virtual environment..."
source "$VENV_NAME/bin/activate"

# Upgrade pip
print_status "Upgrading pip..."
pip install --upgrade pip

# Install wheel for faster installations
print_status "Installing wheel for faster package installations..."
pip install wheel

# Install requirements
if [ -f "requirements.txt" ]; then
    print_status "Installing requirements from requirements.txt..."
    pip install -r requirements.txt
    
    # Check if installation was successful
    if [ $? -eq 0 ]; then
        print_status "All requirements installed successfully!"
    else
        print_error "Some packages failed to install. Please check the error messages above."
        exit 1
    fi
else
    print_error "requirements.txt not found in the current directory."
    exit 1
fi

# Install Jupyter kernel for the virtual environment
print_status "Installing Jupyter kernel for this environment..."
python -m ipykernel install --user --name="crrt-epidemiology" --display-name="CRRT Epidemiology (Python $PYTHON_VERSION)"

# Create directories if they don't exist
print_status "Creating project directories..."
mkdir -p output/final/graphs
mkdir -p output/intermediate
mkdir -p logs

# Summary
echo
print_status "===================================="
print_status "Setup completed successfully!"
print_status "===================================="
echo
echo "To activate the virtual environment, run:"
echo "    source $VENV_NAME/bin/activate"
echo
echo "To deactivate the virtual environment, run:"
echo "    deactivate"
echo
echo "To use this environment in Jupyter, select the kernel:"
echo "    'CRRT Epidemiology (Python $PYTHON_VERSION)'"