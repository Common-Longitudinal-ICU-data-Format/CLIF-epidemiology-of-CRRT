#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONIOENCODING=utf-8

# Check config exists
if [ ! -f "$SCRIPT_DIR/config/config.json" ]; then
    echo "ERROR: config/config.json not found."
    echo "Copy the template and edit it:"
    echo "  cp config/config_template.json config/config.json"
    exit 1
fi

# Set up logging — capture all output to a timestamped log file + console
LOG_DIR="$SCRIPT_DIR/output/final"
mkdir -p "$LOG_DIR"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
SITE_NAME=$(python3 -c "import json; print(json.load(open('$SCRIPT_DIR/config/config.json'))['site_name'])" 2>/dev/null || echo "unknown")
LOG_FILE="$LOG_DIR/${SITE_NAME}_pipeline_${TIMESTAMP}.log"

exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== CRRT Epidemiology Pipeline ==="
echo "  Started: $(date)"
echo "  Site: $SITE_NAME"
echo "  Log: $LOG_FILE"
echo ""

PIPELINE_START=$SECONDS

# Scripts use relative paths like ../config/config.json, so run from code/
cd "$SCRIPT_DIR/code"

PYTHON_STEPS=(
    "00_cohort.py"
    "01_create_wide_df.py"
    "02_construct_crrt_tableone.py"
    "03_crrt_visualizations.py"
    "04_build_msm_competing_risk_df.py"
)

R_STEPS=(
    "05_PSM_IPTW_CRRT_dose.R"
    "05b_dose_response_analysis.R"
    "06_time_varying_MSM.R"
    "06b_time_varying_MSM_sensitivity.R"
)

echo "=== Descriptive + MSM data prep (Python) ==="
echo ""
for step in "${PYTHON_STEPS[@]}"; do
    name=$(basename "$step" .py)
    STEP_START=$SECONDS
    echo "--- Running $name ---"
    uv run python "$step"
    elapsed=$(( SECONDS - STEP_START ))
    echo "--- $name complete (${elapsed}s) ---"
    echo ""
done

# Load R module on HPC systems (no-op if not available)
if command -v module &>/dev/null; then
    module load R 2>/dev/null || true
fi

echo "=== Causal inference (R) ==="
echo ""
R_FAILED=()
for step in "${R_STEPS[@]}"; do
    name=$(basename "$step" .R)
    STEP_START=$SECONDS
    echo "--- Running $name ---"
    if Rscript --no-init-file "$SCRIPT_DIR/code/$step"; then
        elapsed=$(( SECONDS - STEP_START ))
        echo "--- $name complete (${elapsed}s) ---"
    else
        elapsed=$(( SECONDS - STEP_START ))
        echo "--- $name FAILED (${elapsed}s) ---"
        R_FAILED+=("$step")
    fi
    echo ""
done

total_elapsed=$(( SECONDS - PIPELINE_START ))
total_min=$(( total_elapsed / 60 ))
total_sec=$(( total_elapsed % 60 ))

echo "=== Pipeline complete ==="
echo "  Finished: $(date)"
echo "  Total time: ${total_min}m ${total_sec}s"
echo "  Results: output/final/"
echo "  Log: $LOG_FILE"

if [ ${#R_FAILED[@]} -gt 0 ]; then
    echo ""
    echo "=== WARNING: ${#R_FAILED[@]} R script(s) failed ==="
    echo "  The following R scripts did not complete:"
    for f in "${R_FAILED[@]}"; do
        echo "    - $f"
    done
    echo ""
    echo "  To re-run manually from the project root directory:"
    echo ""
    for f in "${R_FAILED[@]}"; do
        echo "    Rscript --no-init-file code/$f"
    done
    echo ""
    echo "  Required R scripts and their outputs:"
    echo "    05_PSM_IPTW_CRRT_dose.R        -> output/final/psm_iptw/"
    echo "    05b_dose_response_analysis.R    -> output/final/psm_iptw/ (dose-response)"
    echo "    06_time_varying_MSM.R           -> output/final/time_varying/ (primary 12h)"
    echo "    06b_time_varying_MSM_sensitivity.R -> output/final/time_varying_sensitivity/ (24h)"
    echo ""
    echo "  Check the log for error details: $LOG_FILE"
fi
