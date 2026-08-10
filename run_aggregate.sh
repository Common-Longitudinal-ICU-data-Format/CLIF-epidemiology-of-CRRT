#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Coordinator-only: build the pooled manuscript + supplement artifacts (08)
# from every per-site run collected under manuscript_run_results/output_<site>/,
# writing to manuscript_run_results/aggregate/{pooled,figures,tables}.
#
# Each output_<site>/ is symlinked into manuscript_run_results/sites/<KEY> with
# the canonical site key report_core expects (matches SITE_LABELS), then 08 runs
# with CRRT_SITES_ROOT / CRRT_AGG_OUT pointed at those dirs.
#
#   ./run_aggregate.sh            # build the pooled artifacts
#   ./run_aggregate.sh dashboard  # also build the 07 dashboard
# ---------------------------------------------------------------------------
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MRR="$ROOT_DIR/manuscript_run_results"
SITES="$MRR/sites"
AGG="$MRR/aggregate"
mkdir -p "$SITES" "$AGG"

# output_<slug> -> canonical site key (report_core SITE_LABELS). Unknown slugs
# fall back to uppercase; add mixed-case sites to the case below as they join.
# (case, not an associative array, so this works on macOS's default bash 3.2.)
site_key() {
    case "$1" in
        nu) echo NU ;;       ucmc) echo UCMC ;;   ucsf) echo UCSF ;;
        umn) echo UMN ;;     ohsu) echo OHSU ;;   rush) echo Rush ;;
        penn) echo Penn ;;   emory) echo Emory ;; michigan) echo Michigan ;;
        hopkins) echo Hopkins ;;
        *) printf '%s' "$1" | tr '[:lower:]' '[:upper:]' ;;
    esac
}

shopt -s nullglob
for d in "$MRR"/output_*/; do
    base=$(basename "$d"); slug=${base#output_}
    ln -sfn "../$base" "$SITES/$(site_key "$slug")"
done
echo "=== sites linked for aggregation ==="
ls -l "$SITES"
echo

export CRRT_SITES_ROOT="$SITES" CRRT_AGG_OUT="$AGG"
echo "=== building pooled manuscript artifacts (08) -> $AGG ==="
uv run python "$ROOT_DIR/code/08_manuscript_artifacts.py"

if [ "${1:-}" = "dashboard" ]; then
    echo; echo "=== building dashboard (07) ==="
    uv run python "$ROOT_DIR/code/07_build_dashboard.py"
fi
echo "Done -> $AGG"
