#!/usr/bin/env bash

# Enhanced setup_and_run.sh — Interactive CLIF project execution script (Mac/Linux)

# set -o pipefail will still be used for individual commands where needed

# ── ANSI colours ─────────────────────────────────────────────
YELLOW="\033[33m"
CYAN="\033[36m"
GREEN="\033[32m"
RED="\033[31m"
BLUE="\033[34m"
BOLD="\033[1m"
RESET="\033[0m"

# ── Setup logging ──────────────────────────────────────────────────────────────
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
LOG_FILE="${SCRIPT_DIR}/logs/execution_log_${TIMESTAMP}.log"
mkdir -p "${SCRIPT_DIR}/logs"

# Function to log and display
log_echo() {
    echo -e "$1" | tee -a "$LOG_FILE"
}

# Initialize log file with header
echo "CLIF Epidemiology of CRRT Project - Execution Log" > "$LOG_FILE"
echo "Started at: $(date)" >> "$LOG_FILE"
echo "========================================" >> "$LOG_FILE"

show_banner() {
    clear
    log_echo "${CYAN}${BOLD}"
    log_echo "╔══════════════════════════════════════════════════════════════════════════════╗"
    log_echo "║                                                                              ║"
    log_echo "║                              ██████ ██████  ██████  ████████                 ║"
    log_echo "║                             ██      ██   ██ ██   ██    ██                    ║"
    log_echo "║                             ██      ██████  ██████     ██                    ║"
    log_echo "║                             ██      ██   ██ ██   ██    ██                    ║"
    log_echo "║                              ██████ ██   ██ ██   ██    ██                    ║"
    log_echo "║                                                                              ║"
    log_echo "║                                       CLIF PROJECT                           ║"
    log_echo "║                                                                              ║"
    log_echo "╚══════════════════════════════════════════════════════════════════════════════╝"
    log_echo "${RESET}"
    log_echo ""
    log_echo "${GREEN}Welcome to the CLIF Epidemiology of CRRT Project!${RESET}"
    log_echo ""
}

# ── Progress separator ─────────────────────────────────────────────────────────
separator() {
    log_echo "${YELLOW}==================================================${RESET}"
}

# ── Progress bar function ──────────────────────────────────────────────────────
show_progress() {
    local step=$1
    local total=$2
    local description=$3
    
    separator
    log_echo "${CYAN}${BOLD}Step ${step}/${total}: ${description}${RESET}"
    log_echo "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')] Starting: ${description}${RESET}"
    separator
}

# ── Error handler ──────────────────────────────────────────────────────────────
# Global array to track failed steps
FAILED_STEPS=()

handle_error() {
    local exit_code=$?
    local step_name=$1
    
    log_echo ""
    log_echo "${RED}${BOLD}❌ ERROR OCCURRED!${RESET}"
    log_echo "${RED}Step failed: ${step_name}${RESET}"
    log_echo "${RED}Exit code: ${exit_code}${RESET}"
    log_echo "${RED}Check the log file for full details: ${LOG_FILE}${RESET}"
    log_echo "${YELLOW}Continuing with next step...${RESET}"
    log_echo ""
    
    # Add to failed steps array
    FAILED_STEPS+=("$step_name")
    
    return 0  # Continue execution
}

# Export environment variables for unbuffered output
export PYTHONUNBUFFERED=1
export MPLBACKEND=Agg

# ── Main execution flow ────────────────────────────────────────────────────────
show_banner

# Step 1: Create virtual environment
show_progress 1 7 "Create Virtual Environment"
if [ ! -d ".crrt" ]; then
    log_echo "Creating virtual environment (.crrt)..."
    python3 -m venv .crrt 2>&1 | tee -a "$LOG_FILE" || handle_error "Create Virtual Environment"
else
    log_echo "Virtual environment already exists."
fi
log_echo "${GREEN}✅ Completed: Create Virtual Environment${RESET}"

# Step 2: Activate virtual environment
show_progress 2 7 "Activate Virtual Environment"
log_echo "Activating virtual environment..."
source .crrt/bin/activate || handle_error "Activate Virtual Environment"
log_echo "${GREEN}✅ Completed: Activate Virtual Environment${RESET}"

# Step 3: Install dependencies
show_progress 3 7 "Install Dependencies"
log_echo "Upgrading pip..."
python -m pip install --upgrade pip 2>&1 | tee -a "$LOG_FILE" || handle_error "Upgrade pip"

log_echo "Installing dependencies..."
pip install -r requirements.txt 2>&1 | tee -a "$LOG_FILE" || handle_error "Install requirements"
pip install jupyter ipykernel 2>&1 | tee -a "$LOG_FILE" || handle_error "Install jupyter"
log_echo "${GREEN}✅ Completed: Install Dependencies${RESET}"

# Step 4: Register Jupyter kernel
show_progress 4 7 "Register Jupyter Kernel"
log_echo "Registering Jupyter kernel..."
python -m ipykernel install --user --name=.crrt --display-name="Python (CRRT)" 2>&1 | tee -a "$LOG_FILE" || handle_error "Register Jupyter Kernel"
log_echo "${GREEN}✅ Completed: Register Jupyter Kernel${RESET}"

# Step 5: Change to code directory and validate data
show_progress 5 7 "Setup Working Directory & Validate Data"
log_echo "Changing to code directory..."
cd code || handle_error "Change to code directory"

# Check if data path exists
log_echo "Checking data configuration..."
CONFIG_FILE="../config/config.json"
if [ -f "$CONFIG_FILE" ]; then
    DATA_PATH=$(python -c "import json; print(json.load(open('$CONFIG_FILE'))['tables_path'])" 2>/dev/null)
    if [ -n "$DATA_PATH" ]; then
        if [ -d "$DATA_PATH" ]; then
            log_echo "${GREEN}✅ Data path found: $DATA_PATH${RESET}"
        else
            log_echo "${YELLOW}⚠️  Data path not found: $DATA_PATH${RESET}"
            log_echo "${YELLOW}Please ensure the CLIF data is available at this location${RESET}"
            log_echo "${YELLOW}Or update config/config.json with the correct path${RESET}"
            log_echo ""
            read -p "Continue anyway? (y/n): " continue_choice
            if [[ ! "$continue_choice" =~ ^[Yy]$ ]]; then
                log_echo "${BLUE}Exiting. Please set up the data path and try again.${RESET}"
                exit 0
            fi
        fi
    fi
else
    log_echo "${YELLOW}⚠️  Config file not found: $CONFIG_FILE${RESET}"
fi

log_echo "${GREEN}✅ Completed: Setup Working Directory${RESET}"

# Step 6: Execute notebooks
show_progress 6 7 "Execute Analysis Notebooks"

log_echo "Executing 01_cohort_identification.ipynb..."
# Ensure buffer is flushed before executing
sync
# Convert notebook to script, suppress nbconvert messages, then run with Python
log_echo "Converting and executing notebook..."
set -o pipefail  # Make pipes fail if any command fails
jupyter nbconvert --to script --stdout --log-level ERROR 01_cohort_identification.ipynb 2>/dev/null | python -u 2>&1 | tee ../logs/01_cohort_identification.log | tee -a "$LOG_FILE"
# Check if the pipeline failed
if [ $? -ne 0 ]; then
    handle_error "Execute 01_cohort_identification.ipynb"
    log_echo "${RED}❌ Failed: 01_cohort_identification.ipynb${RESET}"
else
    log_echo "${GREEN}✅ Completed: 01_cohort_identification.ipynb${RESET}"
fi

log_echo ""
log_echo "Executing 02_analysis_summary.ipynb..."
# Convert notebook to script, suppress nbconvert messages, then run with Python
log_echo "Converting and executing notebook..."
jupyter nbconvert --to script --stdout --log-level ERROR 02_analysis_summary.ipynb 2>/dev/null | python -u 2>&1 | tee ../logs/02_analysis_summary.log | tee -a "$LOG_FILE"
# Check if the pipeline failed
if [ $? -ne 0 ]; then
    handle_error "Execute 02_analysis_summary.ipynb"
    log_echo "${RED}❌ Failed: 02_analysis_summary.ipynb${RESET}"
else
    log_echo "${GREEN}✅ Completed: 02_analysis_summary.ipynb${RESET}"
fi

# Step 7: Analysis Complete
show_progress 7 7 "Analysis Complete"
log_echo "${GREEN}✅ All notebooks executed successfully!${RESET}"

# Final summary
separator
log_echo "${CYAN}${BOLD}📋 EXECUTION SUMMARY${RESET}"
separator

# Display success/failure summary
if [ ${#FAILED_STEPS[@]} -eq 0 ]; then
    log_echo "${GREEN}${BOLD}🎉 SUCCESS! All CRRT epidemiology analysis steps completed successfully!${RESET}"
else
    log_echo "${YELLOW}${BOLD}⚠️  PARTIAL SUCCESS: Some steps failed${RESET}"
    log_echo ""
    log_echo "${RED}${BOLD}Failed steps:${RESET}"
    for step in "${FAILED_STEPS[@]}"; do
        log_echo "${RED}  ❌ $step${RESET}"
    done
    log_echo ""
    log_echo "${YELLOW}Please check the individual log files for error details${RESET}"
fi

log_echo ""
log_echo "${BLUE}📊 Results saved to: output/${RESET}"
log_echo "${BLUE}📝 Full log saved to: ${LOG_FILE}${RESET}"
log_echo "${BLUE}📄 Individual logs in: logs/${RESET}"
separator


log_echo ""
log_echo "${GREEN}Thank you for running the CLIF Epidemiology of CRRT Project!${RESET}"
read -rp "Press [Enter] to exit..."