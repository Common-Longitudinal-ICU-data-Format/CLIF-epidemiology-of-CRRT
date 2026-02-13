"""
CRRT Visualizations — Clinical variables over the first 7 days post-initiation.

1. CRRT dose (ml/kg/hr) hourly.
2. Lab distributions (lactate, pH, bicarbonate, potassium, phosphate) in 12h bins.
3. Mean arterial pressure (MAP) in 12h bins.
4. Respiratory support: FiO2 and IMV proportion in 12h bins.

Saves aggregated (non-patient-level) data and figures to output/final/graphs/.

Usage: uv run python code/crrt_visualizations.py
"""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
project_root = Path(__file__).resolve().parent.parent
with open(project_root / "config" / "config.json") as f:
    config = json.load(f)

INTERMEDIATE_DIR = project_root / "output" / "intermediate"
GRAPHS_DIR = project_root / "output" / "final" / "graphs"
GRAPHS_DIR.mkdir(parents=True, exist_ok=True)

SITE_NAME = config["site_name"]
HAS_CRRT_SETTINGS = config.get("has_crrt_settings", True)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
print("Loading data …")
crrt_initiation = pd.read_parquet(INTERMEDIATE_DIR / "crrt_initiation.parquet")
index_crrt_df = pd.read_parquet(INTERMEDIATE_DIR / "index_crrt_df.parquet")
eb_map = index_crrt_df[["hospitalization_id", "encounter_block"]].copy()

# Load only needed columns from wide_df
if HAS_CRRT_SETTINGS:
    crrt_cols = [
        "crrt_blood_flow_rate", "crrt_dialysate_flow_rate",
        "crrt_pre_filter_replacement_fluid_rate",
        "crrt_post_filter_replacement_fluid_rate",
        "crrt_ultrafiltration_out", "crrt_mode_category",
    ]
else:
    crrt_cols = []
lab_cols = ["lab_lactate", "lab_ph_arterial", "lab_bicarbonate", "lab_potassium", "lab_phosphate"]
other_cols = ["vital_map", "med_cont_nee", "resp_fio2_set", "resp_device_category"]
available = {f.name for f in pq.read_schema(INTERMEDIATE_DIR / "wide_df.parquet")}
needed = ["hospitalization_id", "event_dttm"] + \
    [c for c in crrt_cols + lab_cols + other_cols if c in available]

wide_df = pd.read_parquet(INTERMEDIATE_DIR / "wide_df.parquet", columns=needed)
wide_df = wide_df.merge(eb_map, on="hospitalization_id", how="inner")
wide_df = wide_df.merge(
    crrt_initiation[["encounter_block", "crrt_initiation_time"]],
    on="encounter_block", how="inner",
)
wide_df["hours_from_crrt"] = (
    (wide_df["event_dttm"] - wide_df["crrt_initiation_time"]).dt.total_seconds() / 3600
)
print(f"  wide_df: {wide_df.shape}")

# Weight per encounter (for dose calculation)
weight = index_crrt_df[["encounter_block", "weight_kg"]].copy()
wide_df = wide_df.merge(weight, on="encounter_block", how="left")

# ===================================================================
# 1. CRRT Dose over time (hourly, first 7 days = 168 hours)
# ===================================================================
MAX_HOURS = 168
if not HAS_CRRT_SETTINGS:
    print("CRRT settings not available — skipping dose visualization")

if HAS_CRRT_SETTINGS:
    print("Computing hourly CRRT dose …")

    # Filter to first 7 days post-initiation, only rows with CRRT data
    crrt_rows = wide_df[
        (wide_df["hours_from_crrt"] >= 0)
        & (wide_df["hours_from_crrt"] <= MAX_HOURS)
        & (wide_df["crrt_mode_category"].notna())
    ].copy()

    # Compute dose same way as 00_cohort: mode-specific total flow rate / weight
    conditions = [
        crrt_rows["crrt_mode_category"] == "cvvhd",
        crrt_rows["crrt_mode_category"] == "cvvh",
        crrt_rows["crrt_mode_category"] == "cvvhdf",
    ]
    choices = [
        crrt_rows["crrt_dialysate_flow_rate"],
        crrt_rows["crrt_pre_filter_replacement_fluid_rate"].fillna(0)
        + crrt_rows["crrt_post_filter_replacement_fluid_rate"].fillna(0),
        crrt_rows["crrt_dialysate_flow_rate"].fillna(0)
        + crrt_rows["crrt_pre_filter_replacement_fluid_rate"].fillna(0)
        + crrt_rows["crrt_post_filter_replacement_fluid_rate"].fillna(0),
    ]
    crrt_rows["total_flow_rate"] = np.select(conditions, choices, default=np.nan)
    crrt_rows["crrt_dose_ml_kg_hr"] = np.where(
        (crrt_rows["weight_kg"] > 0) & (crrt_rows["total_flow_rate"] > 0),
        crrt_rows["total_flow_rate"] / crrt_rows["weight_kg"],
        np.nan,
    )

    # Bin into integer hours and compute hourly aggregates
    crrt_rows["hour_bin"] = crrt_rows["hours_from_crrt"].astype(int)
    hourly_dose = (
        crrt_rows.groupby("hour_bin")["crrt_dose_ml_kg_hr"]
        .agg(["mean", "median", "std", "count"])
        .reset_index()
    )
    hourly_dose.columns = ["hour", "mean_dose", "median_dose", "std_dose", "n_patients"]

    # Save aggregated data
    hourly_dose.to_csv(GRAPHS_DIR / "crrt_dose_hourly.csv", index=False)
    print(f"  Saved: {GRAPHS_DIR / 'crrt_dose_hourly.csv'} ({len(hourly_dose)} hours)")

    # Also compute IQR for plotting
    hourly_dose_iqr = (
        crrt_rows.groupby("hour_bin")["crrt_dose_ml_kg_hr"]
        .agg([lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)])
        .reset_index()
    )
    hourly_dose_iqr.columns = ["hour", "q25", "q75"]
    hourly_dose = hourly_dose.merge(hourly_dose_iqr, on="hour", how="left")

    # Re-save with IQR
    hourly_dose["site"] = SITE_NAME
    hourly_dose.to_csv(GRAPHS_DIR / "crrt_dose_hourly.csv", index=False)

    # Plot (median + IQR — robust to outliers)
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(hourly_dose["hour"], hourly_dose["median_dose"], color="#2196F3", linewidth=1.5,
            label="Median dose")
    ax.fill_between(
        hourly_dose["hour"],
        hourly_dose["q25"],
        hourly_dose["q75"],
        alpha=0.2, color="#2196F3", label="IQR",
    )
    ax.set_xlabel("Hours from CRRT initiation")
    ax.set_ylabel("CRRT Dose (ml/kg/hr)")
    ax.set_title(f"CRRT Dose Over First 7 Days — {SITE_NAME}")
    ax.set_xlim(0, MAX_HOURS)
    ax.axhline(y=25, color="grey", linestyle="--", linewidth=0.8, alpha=0.5, label="25 ml/kg/hr (KDIGO)")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(GRAPHS_DIR / "crrt_dose_over_time.png", dpi=150)
    plt.close(fig)
    print(f"  Saved: {GRAPHS_DIR / 'crrt_dose_over_time.png'}")


# ===================================================================
# 2. Lab distributions over CRRT course
# ===================================================================
print("Computing lab distributions over CRRT course …")

lab_info = {
    "lab_lactate":      ("Lactate (mmol/L)", "lactate"),
    "lab_ph_arterial":  ("Arterial pH", "ph_arterial"),
    "lab_bicarbonate":  ("Bicarbonate (mEq/L)", "bicarbonate"),
    "lab_potassium":    ("Potassium (mEq/L)", "potassium"),
    "lab_phosphate":    ("Phosphate (mg/dL)", "phosphate"),
}

# Time bins: 12-hour intervals, using midpoint as continuous x-axis
BIN_WIDTH_H = 12
time_bins = list(range(0, MAX_HOURS + BIN_WIDTH_H, BIN_WIDTH_H))

post_crrt = wide_df[
    (wide_df["hours_from_crrt"] >= 0) & (wide_df["hours_from_crrt"] <= MAX_HOURS)
].copy()
post_crrt["hour_bin"] = (post_crrt["hours_from_crrt"] // BIN_WIDTH_H).astype(int) * BIN_WIDTH_H
post_crrt["hour_mid"] = post_crrt["hour_bin"] + BIN_WIDTH_H / 2  # midpoint

# Compute per-bin aggregates for each lab and save
all_lab_agg = []
for col, (display_name, short_name) in lab_info.items():
    if col not in post_crrt.columns:
        print(f"  Skipping {col} (not in wide_df)")
        continue
    agg = (
        post_crrt.groupby("hour_mid", observed=True)[col]
        .agg(["median", lambda x: x.quantile(0.25), lambda x: x.quantile(0.75), "count"])
        .reset_index()
    )
    agg.columns = ["hour", "median", "q25", "q75", "n"]
    agg["lab"] = short_name
    all_lab_agg.append(agg)

lab_agg_df = pd.concat(all_lab_agg, ignore_index=True)
lab_agg_df["site"] = SITE_NAME
lab_agg_df.to_csv(GRAPHS_DIR / "lab_distributions_over_crrt.csv", index=False)
print(f"  Saved: {GRAPHS_DIR / 'lab_distributions_over_crrt.csv'}")

# Plot: 3x2 grid (5 labs, hide 6th panel)
fig, axes = plt.subplots(3, 2, figsize=(14, 13))
axes = axes.flatten()

for idx, (col, (display_name, short_name)) in enumerate(lab_info.items()):
    if col not in post_crrt.columns:
        continue
    ax = axes[idx]
    agg = lab_agg_df[lab_agg_df["lab"] == short_name]
    ax.plot(agg["hour"], agg["median"], color="#4CAF50", linewidth=1.5, marker="o", markersize=3)
    ax.fill_between(agg["hour"], agg["q25"], agg["q75"], alpha=0.2, color="#4CAF50", label="IQR")
    ax.set_xlim(0, MAX_HOURS)
    ax.set_xlabel("Hours from CRRT initiation")
    ax.set_ylabel(display_name)
    ax.set_title(display_name)
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=8)

# Hide unused 6th panel
axes[5].set_visible(False)

fig.suptitle(f"Lab Distributions Over CRRT Course (First 7 Days) — {SITE_NAME}",
             fontsize=13, y=1.01)
fig.tight_layout()
fig.savefig(GRAPHS_DIR / "lab_distributions_over_crrt.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {GRAPHS_DIR / 'lab_distributions_over_crrt.png'}")

# ===================================================================
# 3. MAP over CRRT course
# ===================================================================
if "vital_map" in post_crrt.columns:
    print("Computing MAP over CRRT course …")
    map_agg = (
        post_crrt.groupby("hour_mid", observed=True)["vital_map"]
        .agg(["median", lambda x: x.quantile(0.25), lambda x: x.quantile(0.75), "count"])
        .reset_index()
    )
    map_agg.columns = ["hour", "median", "q25", "q75", "n"]
    map_agg["site"] = SITE_NAME
    map_agg.to_csv(GRAPHS_DIR / "map_over_crrt.csv", index=False)
    print(f"  Saved: {GRAPHS_DIR / 'map_over_crrt.csv'}")

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(map_agg["hour"], map_agg["median"], color="#E91E63", linewidth=1.5,
            marker="o", markersize=3, label="Median MAP")
    ax.fill_between(map_agg["hour"], map_agg["q25"], map_agg["q75"],
                    alpha=0.2, color="#E91E63", label="IQR")
    ax.set_xlabel("Hours from CRRT initiation")
    ax.set_ylabel("MAP (mmHg)")
    ax.set_title(f"Mean Arterial Pressure Over CRRT Course (First 7 Days) — {SITE_NAME}")
    ax.set_xlim(0, MAX_HOURS)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(GRAPHS_DIR / "map_over_crrt.png", dpi=150)
    plt.close(fig)
    print(f"  Saved: {GRAPHS_DIR / 'map_over_crrt.png'}")
else:
    print("  Skipping MAP figure (vital_map not in wide_df)")

# ===================================================================
# 3b. Norepinephrine Equivalence (NEE) over CRRT course
# ===================================================================
if "med_cont_nee" in post_crrt.columns:
    print("Computing NEE over CRRT course …")
    nee_agg = (
        post_crrt.groupby("hour_mid", observed=True)["med_cont_nee"]
        .agg(["median", lambda x: x.quantile(0.25), lambda x: x.quantile(0.75), "count"])
        .reset_index()
    )
    nee_agg.columns = ["hour", "median", "q25", "q75", "n"]
    nee_agg["site"] = SITE_NAME
    nee_agg.to_csv(GRAPHS_DIR / "nee_over_crrt.csv", index=False)
    print(f"  Saved: {GRAPHS_DIR / 'nee_over_crrt.csv'}")

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(nee_agg["hour"], nee_agg["median"], color="#FF5722", linewidth=1.5,
            marker="o", markersize=3, label="Median NEE")
    ax.fill_between(nee_agg["hour"], nee_agg["q25"], nee_agg["q75"],
                    alpha=0.2, color="#FF5722", label="IQR")
    ax.set_xlabel("Hours from CRRT initiation")
    ax.set_ylabel("NEE (mcg/kg/min)")
    ax.set_title(f"Norepinephrine Equivalence Over CRRT Course (First 7 Days) — {SITE_NAME}")
    ax.set_xlim(0, MAX_HOURS)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(GRAPHS_DIR / "nee_over_crrt.png", dpi=150)
    plt.close(fig)
    print(f"  Saved: {GRAPHS_DIR / 'nee_over_crrt.png'}")
else:
    print("  Skipping NEE figure (med_cont_nee not in wide_df)")

# ===================================================================
# 4. Respiratory support (FiO2 + IMV proportion)
# ===================================================================
print("Computing respiratory support over CRRT course …")

resp_dfs = []

# Panel 1: FiO2 time-series
if "resp_fio2_set" in post_crrt.columns:
    fio2_agg = (
        post_crrt.groupby("hour_mid", observed=True)["resp_fio2_set"]
        .agg(["median", lambda x: x.quantile(0.25), lambda x: x.quantile(0.75), "count"])
        .reset_index()
    )
    fio2_agg.columns = ["hour", "median", "q25", "q75", "n"]
    fio2_agg["variable"] = "fio2"
    resp_dfs.append(fio2_agg)

# Panel 2: IMV proportion
if "resp_device_category" in post_crrt.columns:
    resp_any = post_crrt[post_crrt["resp_device_category"].notna()].copy()
    resp_any["is_imv"] = (resp_any["resp_device_category"] == "imv").astype(int)
    imv_agg = (
        resp_any.groupby("hour_mid", observed=True)["is_imv"]
        .agg(["mean", "count"])
        .reset_index()
    )
    imv_agg.columns = ["hour", "imv_proportion", "n"]
    imv_agg["imv_proportion"] *= 100  # convert to percentage
    imv_agg["variable"] = "imv_proportion"
    resp_dfs.append(imv_agg)

if resp_dfs:
    resp_agg_df = pd.concat(resp_dfs, ignore_index=True)
    resp_agg_df["site"] = SITE_NAME
    resp_agg_df.to_csv(GRAPHS_DIR / "respiratory_over_crrt.csv", index=False)
    print(f"  Saved: {GRAPHS_DIR / 'respiratory_over_crrt.csv'}")

    fig, axes = plt.subplots(2, 1, figsize=(12, 9))

    # Top panel: FiO2
    if "resp_fio2_set" in post_crrt.columns:
        ax = axes[0]
        fio2 = resp_agg_df[resp_agg_df["variable"] == "fio2"]
        ax.plot(fio2["hour"], fio2["median"], color="#FF9800", linewidth=1.5,
                marker="o", markersize=3, label="Median FiO2")
        ax.fill_between(fio2["hour"], fio2["q25"], fio2["q75"],
                        alpha=0.2, color="#FF9800", label="IQR")
        ax.set_xlim(0, MAX_HOURS)
        ax.set_xlabel("Hours from CRRT initiation")
        ax.set_ylabel("FiO2 (fraction)")
        ax.set_title("FiO2 Over CRRT Course")
        ax.legend(fontsize=8)
        ax.grid(axis="y", alpha=0.3)
    else:
        axes[0].set_visible(False)

    # Bottom panel: IMV proportion
    if "resp_device_category" in post_crrt.columns:
        ax = axes[1]
        imv = resp_agg_df[resp_agg_df["variable"] == "imv_proportion"]
        ax.plot(imv["hour"], imv["imv_proportion"], color="#9C27B0", linewidth=1.5,
                marker="o", markersize=3)
        ax.set_xlim(0, MAX_HOURS)
        ax.set_ylim(0, 100)
        ax.set_xlabel("Hours from CRRT initiation")
        ax.set_ylabel("Patients on IMV (%)")
        ax.set_title("IMV Proportion Over CRRT Course")
        ax.grid(axis="y", alpha=0.3)
    else:
        axes[1].set_visible(False)

    fig.suptitle(f"Respiratory Support Over CRRT Course (First 7 Days) — {SITE_NAME}",
                 fontsize=13, y=1.01)
    fig.tight_layout()
    fig.savefig(GRAPHS_DIR / "respiratory_over_crrt.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {GRAPHS_DIR / 'respiratory_over_crrt.png'}")
else:
    print("  Skipping respiratory figure (no resp columns in wide_df)")

# ===================================================================
# 5. Patient State Stacked Area (On IMV, Off IMV, Dead, Discharged)
# ===================================================================
print("Computing patient state proportions over CRRT course …")

import matplotlib.colors as mcolors

outcomes_df = pd.read_parquet(INTERMEDIATE_DIR / "outcomes_df.parquet")

# Merge outcome timing with crrt_initiation_time
# Use `died` (discharge_category in expired/hospice) rather than the stricter
# `in_hosp_death` so that all expired/hospice patients count as dead.
outcome_cols = ["encounter_block", "in_hosp_death", "died", "death_dttm", "last_vital_dttm"]
if "discharge_dttm" in outcomes_df.columns:
    outcome_cols.append("discharge_dttm")
outcome_timing = outcomes_df[outcome_cols].merge(
    crrt_initiation[["encounter_block", "crrt_initiation_time"]],
    on="encounter_block", how="inner",
)

# Hours from CRRT to death (for patients who died — expired or hospice)
outcome_timing["hours_to_death"] = np.where(
    outcome_timing["died"] == 1,
    (outcome_timing["death_dttm"].fillna(outcome_timing["last_vital_dttm"])
     - outcome_timing["crrt_initiation_time"]).dt.total_seconds() / 3600,
    np.nan,
)

# Hours from CRRT to discharge alive (for true survivors only)
# Use actual discharge_dttm from hospitalization table when available,
# fall back to last_vital_dttm
discharge_col = "discharge_dttm" if "discharge_dttm" in outcome_timing.columns else "last_vital_dttm"
outcome_timing["hours_to_discharge"] = np.where(
    outcome_timing["died"] == 0,
    (outcome_timing[discharge_col]
     - outcome_timing["crrt_initiation_time"]).dt.total_seconds() / 3600,
    np.nan,
)

# --- Get hourly IMV status from wide_df respiratory rows ---
if "resp_device_category" in wide_df.columns:
    resp_rows = wide_df[
        (wide_df["resp_device_category"].notna())
        & (wide_df["hours_from_crrt"] >= 0)
        & (wide_df["hours_from_crrt"] <= MAX_HOURS)
    ].copy()

    resp_rows["hour"] = resp_rows["hours_from_crrt"].astype(int)
    resp_rows["is_imv"] = (resp_rows["resp_device_category"] == "imv").astype(int)

    # Last IMV status per patient per hour
    hourly_imv = (
        resp_rows.sort_values(["encounter_block", "hours_from_crrt"])
        .groupby(["encounter_block", "hour"])["is_imv"]
        .last()
        .reset_index()
    )
else:
    print("  Warning: resp_device_category not in wide_df — IMV status unknown")
    hourly_imv = pd.DataFrame(columns=["encounter_block", "hour", "is_imv"])

# --- Build per-hour state for all patients ---
all_patients = crrt_initiation["encounter_block"].unique()
n_patients = len(all_patients)
hours = list(range(0, MAX_HOURS + 1))

# Pre-compute outcome sets for each hour to avoid repeated filtering
dead_ebs = outcome_timing.loc[outcome_timing["hours_to_death"].notna(), ["encounter_block", "hours_to_death"]]
discharged_ebs = outcome_timing.loc[outcome_timing["hours_to_discharge"].notna(), ["encounter_block", "hours_to_discharge"]]

# Disposition lookup for all patients
discharge_disposition = outcomes_df.set_index("encounter_block")["discharge_category"].to_dict()

# Disposition buckets — separate maps for dead vs discharged-alive
DEAD_DISPOSITION_MAP = {
    "expired": "expired",
    "hospice": "hospice",
}
ALIVE_DISPOSITION_MAP = {
    "home": "home",
    "skilled nursing facility (snf)": "snf",
    "long term care hospital (ltach)": "ltach",
    "acute inpatient rehab facility": "rehab",
    "acute care hospital": "other_facility",
    "against medical advice (ama)": "other_facility",
    "psychiatric hospital": "other_facility",
    "jail": "other_facility",
}

state_counts = []
for h in hours:
    # Outcome states (cumulative, absorbing)
    dead_by_h = set(dead_ebs.loc[dead_ebs["hours_to_death"] <= h, "encounter_block"])
    discharged_by_h = set(discharged_ebs.loc[discharged_ebs["hours_to_discharge"] <= h, "encounter_block"])

    # IMV status at this hour
    imv_at_h = hourly_imv[hourly_imv["hour"] == h].set_index("encounter_block")["is_imv"]

    n_dead = 0
    n_discharged = 0
    n_imv = 0
    n_off_imv = 0
    # Dead disposition counters
    dead_disp = {"expired": 0, "hospice": 0}
    # Alive discharge disposition counters
    alive_disp = {"home": 0, "snf": 0, "ltach": 0, "rehab": 0, "other_facility": 0}

    for eb in all_patients:
        if eb in dead_by_h:
            n_dead += 1
            raw_cat = discharge_disposition.get(eb, "")
            bucket = DEAD_DISPOSITION_MAP.get(raw_cat, "expired")
            dead_disp[bucket] += 1
        elif eb in discharged_by_h:
            n_discharged += 1
            raw_cat = discharge_disposition.get(eb, "")
            bucket = ALIVE_DISPOSITION_MAP.get(raw_cat, "other_facility")
            alive_disp[bucket] += 1
        elif eb in imv_at_h.index:
            val = imv_at_h.at[eb] if not isinstance(imv_at_h.loc[eb], pd.Series) else imv_at_h.loc[eb].iloc[-1]
            if val == 1:
                n_imv += 1
            else:
                n_off_imv += 1
        else:
            n_off_imv += 1

    state_counts.append({
        "hour": h,
        "n_total": n_patients,
        "n_dead": n_dead,
        "n_discharged": n_discharged,
        "n_imv": n_imv,
        "n_off_imv": n_off_imv,
        "dead_expired": dead_disp["expired"],
        "dead_hospice": dead_disp["hospice"],
        "discharged_home": alive_disp["home"],
        "discharged_snf": alive_disp["snf"],
        "discharged_ltach": alive_disp["ltach"],
        "discharged_rehab": alive_disp["rehab"],
        "discharged_other_facility": alive_disp["other_facility"],
    })

state_df = pd.DataFrame(state_counts)
state_df["prop_dead"] = state_df["n_dead"] / state_df["n_total"] * 100
state_df["prop_discharged"] = state_df["n_discharged"] / state_df["n_total"] * 100
state_df["prop_imv"] = state_df["n_imv"] / state_df["n_total"] * 100
state_df["prop_off_imv"] = state_df["n_off_imv"] / state_df["n_total"] * 100

state_df["site"] = SITE_NAME

# Save aggregate data
state_df.to_csv(GRAPHS_DIR / "patient_state_over_crrt.csv", index=False)
print(f"  Saved: {GRAPHS_DIR / 'patient_state_over_crrt.csv'}")

# --- Stacked area plot ---
fig, ax = plt.subplots(figsize=(14, 7))

colors = {
    "On IMV": mcolors.to_rgba("#ff8a80", alpha=0.7),
    "Off IMV": mcolors.to_rgba("#90caf9", alpha=0.7),
    "Discharged Alive": mcolors.to_rgba("#ce93d8", alpha=0.7),
    "Dead": mcolors.to_rgba("#bdbdbd", alpha=0.7),
}

ax.stackplot(
    state_df["hour"],
    state_df["prop_imv"],
    state_df["prop_off_imv"],
    state_df["prop_discharged"],
    state_df["prop_dead"],
    labels=list(colors.keys()),
    colors=list(colors.values()),
)

ax.set_xlabel("Hours from CRRT Initiation", fontsize=12)
ax.set_ylabel("Proportion of Patients (%)", fontsize=12)
ax.set_title(f"Patient State Over CRRT Course (First 7 Days) — {SITE_NAME}", fontsize=13)
ax.set_xlim(0, MAX_HOURS)
ax.set_ylim(0, 100)
ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=10)
ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(GRAPHS_DIR / "patient_state_over_crrt.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {GRAPHS_DIR / 'patient_state_over_crrt.png'}")

# ===================================================================
# 6. CRRT Course Summary (milestone table + time-to-event stats)
# ===================================================================
print("Generating CRRT course summary …")

# --- 6a: Hourly milestones (every hour, with discharge disposition) ---
milestone_rows = []
for _, row in state_df.iterrows():
    h = int(row["hour"])
    milestone_rows.append({
        "hour": h,
        "n_total": int(row["n_total"]),
        "n_still_hospitalized": int(row["n_imv"] + row["n_off_imv"]),
        "n_on_imv": int(row["n_imv"]),
        "n_off_imv": int(row["n_off_imv"]),
        "n_dead": int(row["n_dead"]),
        "dead_expired": int(row["dead_expired"]),
        "dead_hospice": int(row["dead_hospice"]),
        "n_discharged": int(row["n_discharged"]),
        "discharged_home": int(row["discharged_home"]),
        "discharged_snf": int(row["discharged_snf"]),
        "discharged_ltach": int(row["discharged_ltach"]),
        "discharged_rehab": int(row["discharged_rehab"]),
        "discharged_other_facility": int(row["discharged_other_facility"]),
        "pct_on_imv": round(row["prop_imv"], 1),
        "pct_off_imv": round(row["prop_off_imv"], 1),
        "pct_dead": round(row["prop_dead"], 1),
        "pct_discharged": round(row["prop_discharged"], 1),
    })
milestone_df = pd.DataFrame(milestone_rows)

# --- 6b: Time-to-event statistics ---
time_to_death = outcome_timing.loc[
    outcome_timing["hours_to_death"].notna() & (outcome_timing["hours_to_death"] > 0),
    "hours_to_death",
]
time_to_discharge = outcome_timing.loc[
    outcome_timing["hours_to_discharge"].notna() & (outcome_timing["hours_to_discharge"] > 0),
    "hours_to_discharge",
]

# Time to first extubation (first hour where patient transitions from IMV to off-IMV)
imv_sorted = hourly_imv.sort_values(["encounter_block", "hour"])
imv_sorted["was_imv"] = imv_sorted.groupby("encounter_block")["is_imv"].shift(1)
first_extubation = (
    imv_sorted[(imv_sorted["was_imv"] == 1) & (imv_sorted["is_imv"] == 0)]
    .groupby("encounter_block")["hour"]
    .first()
)

def _tte_row(name, series, unit="hours"):
    """Build a summary row for a time-to-event variable."""
    return {
        "metric": name,
        "n": int(series.count()),
        "median": round(series.median(), 1),
        "q25": round(series.quantile(0.25), 1),
        "q75": round(series.quantile(0.75), 1),
        "unit": unit,
    }

tte_rows = [
    _tte_row("Time to death (from CRRT)", time_to_death),
    _tte_row("Time to discharge alive (from CRRT)", time_to_discharge),
    _tte_row("Time to first extubation (from CRRT)", first_extubation),
    _tte_row("Time to death (days)", time_to_death / 24, unit="days"),
    _tte_row("Time to discharge alive (days)", time_to_discharge / 24, unit="days"),
]
tte_df = pd.DataFrame(tte_rows)

# Save both tables
milestone_df["site"] = SITE_NAME
tte_df["site"] = SITE_NAME
milestone_df.to_csv(GRAPHS_DIR / "crrt_course_milestones.csv", index=False)
tte_df.to_csv(GRAPHS_DIR / "crrt_course_time_to_event.csv", index=False)
print(f"  Saved: {GRAPHS_DIR / 'crrt_course_milestones.csv'}")
print(f"  Saved: {GRAPHS_DIR / 'crrt_course_time_to_event.csv'}")

# Print summary to console
print("\n" + "=" * 70)
print("CRRT COURSE SUMMARY")
print("=" * 70)
print(f"\nTotal patients: {n_patients}")
print(f"  Died in-hospital: {int(outcome_timing['in_hosp_death'].sum())} "
      f"({outcome_timing['in_hosp_death'].mean() * 100:.1f}%)")
key_hours = [24, 48, 72, 120, 168]
key_milestones = milestone_df[milestone_df["hour"].isin(key_hours)]
print(f"\nPatient state at key milestones:")
print(key_milestones[["hour", "n_still_hospitalized", "n_on_imv",
                       "n_dead", "n_discharged"]].to_string(index=False))
print(f"\nDisposition at day 7 (hour 168):")
day7 = milestone_df[milestone_df["hour"] == 168].iloc[0]
print(f"  Dead — expired: {int(day7['dead_expired'])}, hospice: {int(day7['dead_hospice'])}")
print(f"  Discharged alive — home: {int(day7['discharged_home'])}, "
      f"SNF: {int(day7['discharged_snf'])}, LTACH: {int(day7['discharged_ltach'])}, "
      f"rehab: {int(day7['discharged_rehab'])}, other: {int(day7['discharged_other_facility'])}")
print(f"\nTime-to-event (median [IQR]):")
for _, r in tte_df.iterrows():
    print(f"  {r['metric']}: {r['median']} [{r['q25']}, {r['q75']}] "
          f"{r['unit']} (n={r['n']})")
print("=" * 70)

print("\nDone!")
