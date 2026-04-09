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

echo "=== CRRT Epidemiology Pipeline ==="
echo ""

echo "=== Descriptive + MSM data prep (Python) ==="
echo ""
for step in "${PYTHON_STEPS[@]}"; do
    name=$(basename "$step" .py)
    echo "--- Running $name ---"
    uv run python "$step"
    echo "--- $name complete ---"
    echo ""
done

echo "=== Causal inference (R) ==="
echo ""
for step in "${R_STEPS[@]}"; do
    name=$(basename "$step" .R)
    echo "--- Running $name ---"
    Rscript "$SCRIPT_DIR/code/$step"
    echo "--- $name complete ---"
    echo ""
done

echo "=== Pipeline complete ==="
echo "Results are in output/final/"
