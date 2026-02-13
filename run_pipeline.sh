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

STEPS=(
    "00_cohort.py"
    "01_create_wide_df.py"
    "02_construct_crrt_tableone.py"
    "03_crrt_visualizations.py"
)

echo "=== CRRT Epidemiology Pipeline ==="
echo ""

for step in "${STEPS[@]}"; do
    name=$(basename "$step" .py)
    echo "--- Running $name ---"
    uv run python "$step"
    echo "--- $name complete ---"
    echo ""
done

echo "=== Pipeline complete ==="
echo "Results are in output/final/"
