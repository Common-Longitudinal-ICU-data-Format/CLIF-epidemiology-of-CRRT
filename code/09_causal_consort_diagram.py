"""
Generate a CONSORT-style flow diagram showing cohort narrowing from
descriptive analysis (N=2,136) through causal analysis (N=1,555).

Reads strobe_counts.csv and msm_competing_risk_df.parquet to compute
exact exclusion counts dynamically.

Usage:
    uv run python code/09_causal_consort_diagram.py
"""

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch

# ---------------------------------------------------------------------------
# 0. Configuration
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = REPO_ROOT / "config" / "config.json"
with open(CONFIG_PATH) as f:
    config = json.load(f)

site = config["site_name"]
OUTPUT_DIR = REPO_ROOT / "output" / "final" / "psm_iptw"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Load data and compute counts
# ---------------------------------------------------------------------------
strobe = pd.read_csv(REPO_ROOT / "output" / "final" / "crrt_epi" / "strobe_counts.csv")
strobe_dict = dict(zip(strobe["counter"], strobe["value"]))

outcomes_df = pd.read_parquet(REPO_ROOT / "output" / "intermediate" / "outcomes_df.parquet")
index_crrt = pd.read_parquet(REPO_ROOT / "output" / "intermediate" / "index_crrt_df.parquet")
msm_df = pd.read_parquet(REPO_ROOT / "output" / "intermediate" / "msm_competing_risk_df.parquet")

# --- Descriptive cohort counts (from STROBE) ---
n_total_hosp = int(strobe_dict.get("1b_after_stitching", strobe_dict.get("1_adult_hospitalizations", 0)))
n_crrt = int(strobe_dict["2_crrt_blocks"])
n_no_esrd = int(strobe_dict["3_encounter_blocks_without_esrd"])
n_with_weight = int(strobe_dict["4_encounter_blocks_with_weight"])
n_with_settings = int(strobe_dict.get("5_encounter_blocks_with_crrt_settings", n_with_weight))
n_with_labs = int(strobe_dict["6_encounter_blocks_with_required_labs"])
n_descriptive = n_with_labs  # Final descriptive cohort

# --- Step 04 exclusion breakdown ---
# Reconstruct the short_mask logic from 04_build_msm_competing_risk_df.py
merged = outcomes_df.merge(
    index_crrt[["hospitalization_id", "encounter_block", "duration_days", "crrt_initiation_time"]],
    on=["hospitalization_id", "encounter_block"],
    how="left",
)
merged["crrt_initiation_time"] = pd.to_datetime(merged["crrt_initiation_time"])
merged["death_dttm"] = pd.to_datetime(merged["death_dttm"])
merged["time_crrt_to_death_days"] = (
    (merged["death_dttm"] - merged["crrt_initiation_time"]).dt.total_seconds() / 86400
)

died_24h = (merged["in_hosp_death"] == 1) & (merged["time_crrt_to_death_days"] <= 1.0)
crrt_short = merged["duration_days"] < 1.0
short_mask = died_24h | crrt_short
n_excluded_short = int(short_mask.sum())

remaining_after_short = merged[~short_mask]

# SCUF-only exclusion (check if any were excluded beyond short_mask)
n_after_short = len(remaining_after_short)
n_msm = len(msm_df)
n_excluded_scuf = n_after_short - n_msm  # 0 at UCMC

# --- Step 05: MICE imputes remaining NAs, so no further exclusion ---
# oxygenation_index_0 has NAs in the raw parquet, but MICE (pmm, m=5)
# fills them in the R script before analysis. The drop_na() is a safety
# net that drops 0 rows after imputation.
n_missing_oi_pre_mice = int(msm_df["oxygenation_index_0"].isna().sum())
n_causal = n_msm  # No further exclusion after MICE

# Dose group breakdown (on full msm_df since MICE fills NAs)
n_high_dose = int((msm_df["crrt_dose_ml_kg_hr_0"] >= 30).sum())
n_low_dose = int((msm_df["crrt_dose_ml_kg_hr_0"] < 30).sum())

print(f"=== {site} CONSORT Flow ===")
print(f"Adult hospitalizations:     {n_total_hosp:,}")
print(f"  -> No CRRT:               -{n_total_hosp - n_crrt:,}")
print(f"CRRT hospitalizations:      {n_crrt:,}")
print(f"  -> ESRD:                  -{n_crrt - n_no_esrd:,}")
print(f"  -> Missing weight:        -{n_no_esrd - n_with_weight:,}")
if n_with_weight != n_with_settings:
    print(f"  -> Missing CRRT settings: -{n_with_weight - n_with_settings:,}")
print(f"  -> Missing labs:          -{n_with_settings - n_with_labs:,}")
print(f"Descriptive cohort:         {n_descriptive:,}")
print(f"  -> Died/off CRRT <=24h:   -{n_excluded_short:,}")
if n_excluded_scuf > 0:
    print(f"  -> SCUF-only:             -{n_excluded_scuf:,}")
print(f"Causal analysis eligible:   {n_msm:,}")
print(f"  (oxygenation_index NAs before MICE: {n_missing_oi_pre_mice})")
print(f"Final causal cohort:        {n_causal:,}  (MICE imputes remaining NAs)")
print(f"  High dose (>=30):         {n_high_dose:,}")
print(f"  Low dose (<30):           {n_low_dose:,}")

# ---------------------------------------------------------------------------
# 2. Build the CONSORT figure (matching 00_cohort.py style)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 8))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

# Geometry — identical to 00_cohort.py
box_h = 0.08
box_w = 0.40
x_main_start = 0.05
x_main_center = x_main_start + box_w / 2  # 0.25
x_excl_start = 0.55
excl_arrow_gap = 0.015
v_spacing = 0.14

arrow_props = dict(arrowstyle="->", lw=2, color="black")


def draw_box(x, y, w, h, text, fontsize=11, weight="normal"):
    rect = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.01",
        linewidth=2,
        edgecolor="black",
        facecolor="white",
    )
    ax.add_patch(rect)
    ax.text(x + w / 2, y + h / 2, text,
            ha="center", va="center", fontsize=fontsize, fontweight=weight,
            wrap=True)
    return x + w / 2, y


# ---------------------------------------------------------------------------
# 3. Draw the diagram
# ---------------------------------------------------------------------------

# Title
ax.text(0.5, 0.98, "CRRT Causal Cohort Selection",
        ha="center", va="center", fontsize=16, fontweight="bold")

# --- Row 0: Starting cohort ---
top_y = 0.90 - box_h
draw_box(x_main_start, top_y, box_w, box_h,
         f"All adult hospitalizations\n(2018\u20132024)\nn = {n_total_hosp:,}")

# Build rows data
rows = [
    {"remaining_label": "Remaining hospitalizations\nCRRT hospitalizations",
     "remaining_n": n_crrt,
     "excluded_label": f"Excluded: No CRRT\nn = {n_total_hosp - n_crrt:,}",
     "excluded_n": n_total_hosp - n_crrt},
    {"remaining_label": "Remaining hospitalizations\nAfter ESRD exclusion",
     "remaining_n": n_no_esrd,
     "excluded_label": f"Excluded: ESRD diagnosis\nn = {n_crrt - n_no_esrd:,}",
     "excluded_n": n_crrt - n_no_esrd},
    {"remaining_label": "Remaining hospitalizations\nWith documented weight",
     "remaining_n": n_with_weight,
     "excluded_label": f"Excluded: Missing weight\nn = {n_no_esrd - n_with_weight:,}",
     "excluded_n": n_no_esrd - n_with_weight},
    {"remaining_label": "Remaining hospitalizations\nWith required baseline labs",
     "remaining_n": n_with_labs,
     "excluded_label": f"Excluded: Missing required labs\nn = {n_with_settings - n_with_labs:,}",
     "excluded_n": n_with_settings - n_with_labs},
]

# Add causal exclusion row
excl_short_label = f"Excluded: Died or off CRRT within 24h\nn = {n_excluded_short:,}"
if n_excluded_scuf > 0:
    excl_short_label += f"\n+ SCUF-only: n = {n_excluded_scuf:,}"
rows.append({
    "remaining_label": "Causal analysis cohort",
    "remaining_n": n_causal,
    "excluded_label": excl_short_label,
    "excluded_n": n_excluded_short + n_excluded_scuf,
})

for i, row in enumerate(rows):
    y_parent = top_y if i == 0 else top_y - (i * v_spacing)
    current_y = top_y - ((i + 1) * v_spacing)

    # Vertical center between the two main-column box centers
    box_center_y_top = y_parent + box_h / 2
    box_center_y_bottom = current_y + box_h / 2
    arrow_vertical_center = (box_center_y_top + box_center_y_bottom) / 2

    # Main-column remaining box
    draw_box(x_main_start, current_y, box_w, box_h,
             f"{row['remaining_label']}\nn = {row['remaining_n']:,}")

    # Vertical arrow (top box bottom -> bottom box top)
    ax.annotate("", xy=(x_main_center, current_y + box_h),
                xytext=(x_main_center, y_parent),
                arrowprops=arrow_props)

    # Exclusion box + horizontal arrow (from midpoint of vertical segment)
    if row["excluded_n"] > 0:
        draw_box(x_excl_start, arrow_vertical_center - box_h / 2, box_w, box_h,
                 row["excluded_label"])
        ax.annotate("",
                    xy=(x_excl_start - excl_arrow_gap, arrow_vertical_center),
                    xytext=(x_main_center, arrow_vertical_center),
                    arrowprops=arrow_props, annotation_clip=False)


# ---------------------------------------------------------------------------
# 4. Save
# ---------------------------------------------------------------------------
out_png = OUTPUT_DIR / f"{site}_causal_consort_diagram.png"
out_pdf = OUTPUT_DIR / f"{site}_causal_consort_diagram.pdf"
plt.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
plt.savefig(out_pdf, bbox_inches="tight", facecolor="white")
plt.close(fig)

print(f"\nSaved: {out_png}")
print(f"Saved: {out_pdf}")
