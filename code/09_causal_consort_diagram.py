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
# 2. Build the CONSORT figure
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(13, 12))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

# Geometry
box_h = 0.05
box_w = 0.32
x_main = 0.08
x_main_center = x_main + box_w / 2
x_excl = 0.56
excl_box_w = 0.30
v_spacing = 0.08
arrow_pad = 0.006  # padding so arrowheads don't touch box edges

arrow_main = dict(arrowstyle="-|>", lw=1.8, color="#333333",
                  mutation_scale=12)
arrow_excl = dict(arrowstyle="-|>", lw=1.2, color="#999999",
                  mutation_scale=10)


def draw_box(x, y, w, h, text, fontsize=10, weight="normal", facecolor="white",
             edgecolor="black", linewidth=1.5, zorder=2):
    rect = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.012",
        linewidth=linewidth,
        edgecolor=edgecolor,
        facecolor=facecolor,
        zorder=zorder,
    )
    ax.add_patch(rect)
    ax.text(x + w / 2, y + h / 2, text,
            ha="center", va="center", fontsize=fontsize, fontweight=weight,
            linespacing=1.25, zorder=zorder + 1)
    return x + w / 2, y


def draw_flow_arrow(x, y_from_box, y_to_box):
    """Vertical arrow from bottom of upper box to top of lower box, with padding."""
    ax.annotate(
        "", xy=(x, y_to_box + box_h + arrow_pad),
        xytext=(x, y_from_box - arrow_pad),
        arrowprops=arrow_main, zorder=1,
    )


def draw_excl_arrow_and_box(y_upper_box, y_lower_box, excl_text):
    """Horizontal arrow from right edge of vertical segment to exclusion box."""
    mid_y = (y_upper_box + y_lower_box + box_h) / 2
    # Arrow starts from right edge of main box column, ends at exclusion box
    ax.annotate(
        "", xy=(x_excl, mid_y),
        xytext=(x_main + box_w + arrow_pad, mid_y),
        arrowprops=arrow_excl, annotation_clip=False, zorder=1,
    )
    draw_box(x_excl, mid_y - box_h / 2, excl_box_w, box_h, excl_text,
             fontsize=9, facecolor="#fef2f2", edgecolor="#dc2626", linewidth=1)


# ---------------------------------------------------------------------------
# 3. Draw the diagram
# ---------------------------------------------------------------------------

# Title
fig.text(0.5, 0.965, f"CONSORT Flow Diagram \u2014 {site}",
         ha="center", va="center", fontsize=16, fontweight="bold")
fig.text(0.5, 0.94, "Descriptive Cohort through Causal Analysis",
         ha="center", va="center", fontsize=10.5, color="#555555")

# Row 0: Starting cohort
y = 0.90
draw_box(x_main, y, box_w, box_h,
         f"Adult hospitalizations (2018\u20132024)\nn = {n_total_hosp:,}",
         fontsize=10, weight="bold")

# Row 1: CRRT
y1 = y - v_spacing
draw_flow_arrow(x_main_center, y, y1)
draw_box(x_main, y1, box_w, box_h,
         f"CRRT hospitalizations\nn = {n_crrt:,}")
draw_excl_arrow_and_box(y, y1,
                        f"No CRRT\nn = {n_total_hosp - n_crrt:,}")

# Row 2: No ESRD
y2 = y1 - v_spacing
draw_flow_arrow(x_main_center, y1, y2)
draw_box(x_main, y2, box_w, box_h,
         f"After ESRD exclusion\nn = {n_no_esrd:,}")
draw_excl_arrow_and_box(y1, y2,
                        f"ESRD diagnosis\nn = {n_crrt - n_no_esrd:,}")

# Row 3: Weight
y3 = y2 - v_spacing
draw_flow_arrow(x_main_center, y2, y3)
draw_box(x_main, y3, box_w, box_h,
         f"With documented weight\nn = {n_with_weight:,}")
draw_excl_arrow_and_box(y2, y3,
                        f"Missing weight\nn = {n_no_esrd - n_with_weight:,}")

# Row 4: Labs
y4 = y3 - v_spacing
draw_flow_arrow(x_main_center, y3, y4)
draw_box(x_main, y4, box_w, box_h,
         f"With required baseline labs\nn = {n_with_labs:,}",
         fontsize=10, weight="bold", facecolor="#dcfce7", edgecolor="#16a34a")
if n_with_settings - n_with_labs > 0:
    draw_excl_arrow_and_box(y3, y4,
                            f"Missing required labs\nn = {n_with_settings - n_with_labs:,}")

# Descriptive cohort label
ax.text(x_main + box_w + 0.015, y4 + box_h / 2,
        "Descriptive\nCohort",
        fontsize=9, fontweight="bold", color="#16a34a", va="center")

# Row 5: Causal cohort
y5 = y4 - v_spacing - 0.01
draw_flow_arrow(x_main_center, y4, y5)
draw_box(x_main, y5, box_w, box_h,
         f"Causal analysis cohort\nn = {n_causal:,}",
         fontsize=10, weight="bold", facecolor="#dbeafe", edgecolor="#2563eb")
excl_short_text = f"Died or off CRRT within 24h\nn = {n_excluded_short:,}"
if n_excluded_scuf > 0:
    excl_short_text += f"\n+ SCUF-only: n = {n_excluded_scuf:,}"
draw_excl_arrow_and_box(y4, y5, excl_short_text)

# Causal cohort label
ax.text(x_main + box_w + 0.015, y5 + box_h / 2,
        "Causal\nCohort",
        fontsize=9, fontweight="bold", color="#2563eb", va="center")


# ---------------------------------------------------------------------------
# 4. Save
# ---------------------------------------------------------------------------
out_png = OUTPUT_DIR / f"{site}_causal_consort_diagram.png"
out_pdf = OUTPUT_DIR / f"{site}_causal_consort_diagram.pdf"
plt.subplots_adjust(top=0.92)
plt.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
plt.savefig(out_pdf, bbox_inches="tight", facecolor="white")
plt.close(fig)

print(f"\nSaved: {out_png}")
print(f"Saved: {out_pdf}")
