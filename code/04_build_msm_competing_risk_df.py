"""
Build the MSM / competing-risk analysis dataframe (59 columns).

Extends the archived 04_build_competing_risk_df.py with:
  - Oxygenation index at t=0 and t=12 (P/F or S/F)
  - Norepinephrine equivalent (NEE) at t=0 and t=12
  - IMV status at t=0 and t=12
  - 19 CCI binary components (from clifpy.calculate_cci + ICD-10 cancer split)
  - 30-day censoring (columns named time_to_event_90d / censored_at_90d
    to match R script expectations)
  - Sensitivity analysis columns for 24h/48h MSM: CRRT dose (0-24h, 24-48h),
    labs/SOFA/oxygenation/NEE/IMV at t=24, and 48h eligibility flag

Usage: uv run python code/04_build_msm_competing_risk_df.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import pyarrow.parquet as pq

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root / "code"))

with open(project_root / "config" / "config.json") as f:
    config = json.load(f)

from pipeline_helpers import validate_config, load_intermediate  # noqa: E402
config = validate_config(config)

from sofa_calculator import compute_sofa_polars  # noqa: E402

INTERMEDIATE_DIR = project_root / "output" / "intermediate"
OUTPUT_DIR = project_root / "output" / "intermediate"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

TABLES_PATH = config["tables_path"]
FILE_TYPE = config["file_type"]
TIMEZONE = config["timezone"]

# ===================================================================
# STEP 1: Load intermediate files
# ===================================================================
print("Step 1: Loading intermediate files …")
outcomes_df = load_intermediate(INTERMEDIATE_DIR / "outcomes_df.parquet")
index_crrt_df = load_intermediate(INTERMEDIATE_DIR / "index_crrt_df.parquet")
crrt_initiation = load_intermediate(INTERMEDIATE_DIR / "crrt_initiation.parquet")
tableone_df = load_intermediate(INTERMEDIATE_DIR / "tableone_analysis_df.parquet")

eb_map = index_crrt_df[["hospitalization_id", "encounter_block"]].copy()
print(f"  {len(outcomes_df)} encounters")

# ===================================================================
# STEP 2: Demographics + static columns
# ===================================================================
print("Step 2: Demographics and static columns …")
result = outcomes_df[["encounter_block", "age_at_admission", "sex_category",
                       "race_category", "ethnicity_category"]].copy()

result = result.merge(
    index_crrt_df[["encounter_block", "weight_kg", "crrt_mode_category",
                    "duration_days", "imv_duration_days"]],
    on="encounter_block", how="left",
)
result = result.rename(columns={
    "duration_days": "crrt_duration_days",
})
# Never on IMV → 0 days (NaN means no IMV records, not missing data)
result["imv_duration_days"] = result["imv_duration_days"].fillna(0)
print(f"  {len(result)} rows after merge")

# ===================================================================
# STEP 3: Baseline labs (t=0) from tableone_analysis_df
# ===================================================================
print("Step 3: Baseline labs (t=0) from tableone …")
baseline_lab_map = {
    "creatinine_baseline": "creatinine_0",
    "lactate_baseline": "lactate_0",
    "bicarbonate_baseline": "bicarbonate_0",
    "potassium_baseline": "potassium_0",
}
baseline_cols = ["encounter_block"] + list(baseline_lab_map.keys())
available_baseline = [c for c in baseline_cols if c in tableone_df.columns]
baseline_labs = tableone_df[available_baseline].rename(columns=baseline_lab_map)
result = result.merge(baseline_labs, on="encounter_block", how="left")

# SOFA at t=0
result = result.merge(
    tableone_df[["encounter_block", "sofa_total"]].rename(
        columns={"sofa_total": "sofa_total_0"}
    ),
    on="encounter_block", how="left",
)
print(f"  Baseline labs + SOFA attached")

# ===================================================================
# STEP 4: CRRT dose — at initiation, 0-12h mean, 12-24h mean
# ===================================================================
print("Step 4: Computing CRRT doses from wide_df …")

crrt_cols = [
    "crrt_blood_flow_rate", "crrt_dialysate_flow_rate",
    "crrt_pre_filter_replacement_fluid_rate",
    "crrt_post_filter_replacement_fluid_rate",
    "crrt_mode_category",
]
available = {f.name for f in pq.read_schema(INTERMEDIATE_DIR / "wide_df.parquet")}
needed = ["hospitalization_id", "event_dttm"] + [c for c in crrt_cols if c in available]

wide_df = load_intermediate(INTERMEDIATE_DIR / "wide_df.parquet", columns=needed)
wide_df = wide_df.merge(eb_map, on="hospitalization_id", how="inner")
wide_df = wide_df.merge(
    crrt_initiation[["encounter_block", "crrt_initiation_time"]],
    on="encounter_block", how="inner",
)
wide_df["hours_from_crrt"] = (
    (wide_df["event_dttm"] - wide_df["crrt_initiation_time"]).dt.total_seconds() / 3600
)

# Merge weight for dose calc
wide_df = wide_df.merge(
    index_crrt_df[["encounter_block", "weight_kg"]], on="encounter_block", how="left"
)

# Forward-fill CRRT settings within each patient to cover charting gaps
wide_df = wide_df.sort_values(["encounter_block", "event_dttm"])
crrt_ffill_cols = [c for c in crrt_cols if c in wide_df.columns]
wide_df[crrt_ffill_cols] = wide_df.groupby("encounter_block")[crrt_ffill_cols].ffill(limit=8)
print(f"  Forward-filled CRRT settings (limit=8 rows)")

# Mode-specific dose formula
crrt_rows = wide_df[wide_df["crrt_mode_category"].notna()].copy()
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

# 4a: Dose at initiation
init_dose = (
    crrt_rows[crrt_rows["hours_from_crrt"] >= 0]
    .sort_values("event_dttm")
    .groupby("encounter_block")["crrt_dose_ml_kg_hr"]
    .first()
    .reset_index()
    .rename(columns={"crrt_dose_ml_kg_hr": "crrt_dose_ml_kg_hr_0"})
)

# 4b: Mean dose 0-12h
dose_0_12 = (
    crrt_rows[
        (crrt_rows["hours_from_crrt"] >= 0) & (crrt_rows["hours_from_crrt"] < 12)
    ]
    .groupby("encounter_block")["crrt_dose_ml_kg_hr"]
    .mean()
    .reset_index()
    .rename(columns={"crrt_dose_ml_kg_hr": "crrt_dose_0_12"})
)

# 4c: Mean dose 12-24h
dose_12_24 = (
    crrt_rows[
        (crrt_rows["hours_from_crrt"] >= 12) & (crrt_rows["hours_from_crrt"] < 24)
    ]
    .groupby("encounter_block")["crrt_dose_ml_kg_hr"]
    .mean()
    .reset_index()
    .rename(columns={"crrt_dose_ml_kg_hr": "crrt_dose_12_24"})
)

result = result.merge(init_dose, on="encounter_block", how="left")
result = result.merge(dose_0_12, on="encounter_block", how="left")
result = result.merge(dose_12_24, on="encounter_block", how="left")
# Fallback: if dose 12-24h still missing, carry forward from 0-12h
n_miss_before = result["crrt_dose_12_24"].isna().sum()
result["crrt_dose_12_24"] = result["crrt_dose_12_24"].fillna(result["crrt_dose_0_12"])
n_miss_after = result["crrt_dose_12_24"].isna().sum()
print(f"  Dose at initiation: {init_dose['crrt_dose_ml_kg_hr_0'].notna().sum()} encounters")
print(f"  Dose 0-12h: {dose_0_12['crrt_dose_0_12'].notna().sum()} encounters")
print(f"  Dose 12-24h: {dose_12_24['crrt_dose_12_24'].notna().sum()} encounters"
      f" (+{n_miss_before - n_miss_after} filled from 0-12h)")

# 4d: Mean dose 0-24h (sensitivity)
dose_0_24 = (
    crrt_rows[
        (crrt_rows["hours_from_crrt"] >= 0) & (crrt_rows["hours_from_crrt"] < 24)
    ]
    .groupby("encounter_block")["crrt_dose_ml_kg_hr"]
    .mean()
    .reset_index()
    .rename(columns={"crrt_dose_ml_kg_hr": "crrt_dose_0_24"})
)

# 4e: Mean dose 24-48h (sensitivity)
dose_24_48 = (
    crrt_rows[
        (crrt_rows["hours_from_crrt"] >= 24) & (crrt_rows["hours_from_crrt"] < 48)
    ]
    .groupby("encounter_block")["crrt_dose_ml_kg_hr"]
    .mean()
    .reset_index()
    .rename(columns={"crrt_dose_ml_kg_hr": "crrt_dose_24_48"})
)

result = result.merge(dose_0_24, on="encounter_block", how="left")
result = result.merge(dose_24_48, on="encounter_block", how="left")
# Fallback: if dose 24-48h missing, carry forward from 0-24h
n_miss_24_48_before = result["crrt_dose_24_48"].isna().sum()
result["crrt_dose_24_48"] = result["crrt_dose_24_48"].fillna(result["crrt_dose_0_24"])
n_miss_24_48_after = result["crrt_dose_24_48"].isna().sum()
print(f"  Dose 0-24h: {dose_0_24['crrt_dose_0_24'].notna().sum()} encounters")
print(f"  Dose 24-48h: {dose_24_48['crrt_dose_24_48'].notna().sum()} encounters"
      f" (+{n_miss_24_48_before - n_miss_24_48_after} filled from 0-24h)")

# Flag patients with any non-zero dialysate flow in 0-24h (for SCUF exclusion)
has_dialysate_24h = (
    wide_df[
        (wide_df["hours_from_crrt"] >= 0) & (wide_df["hours_from_crrt"] < 24)
        & (wide_df["crrt_dialysate_flow_rate"].notna())
        & (wide_df["crrt_dialysate_flow_rate"] > 0)
    ]
    .groupby("encounter_block")["crrt_dialysate_flow_rate"]
    .count()
    .reset_index()
    .rename(columns={"crrt_dialysate_flow_rate": "has_dialysate_24h"})
)
has_dialysate_24h["has_dialysate_24h"] = 1
result = result.merge(has_dialysate_24h, on="encounter_block", how="left")
result["has_dialysate_24h"] = result["has_dialysate_24h"].fillna(0).astype(int)
n_no_dialysate = (result["has_dialysate_24h"] == 0).sum()
print(f"  Patients with no dialysate in 0-24h (SCUF-only candidates): {n_no_dialysate}")

# 4f: Dialysate check extended to 48h (sensitivity)
has_dialysate_48h = (
    wide_df[
        (wide_df["hours_from_crrt"] >= 0) & (wide_df["hours_from_crrt"] < 48)
        & (wide_df["crrt_dialysate_flow_rate"].notna())
        & (wide_df["crrt_dialysate_flow_rate"] > 0)
    ]
    .groupby("encounter_block")["crrt_dialysate_flow_rate"]
    .count()
    .reset_index()
    .rename(columns={"crrt_dialysate_flow_rate": "has_dialysate_48h"})
)
has_dialysate_48h["has_dialysate_48h"] = 1
result = result.merge(has_dialysate_48h, on="encounter_block", how="left")
result["has_dialysate_48h"] = result["has_dialysate_48h"].fillna(0).astype(int)
n_no_dialysate_48 = (result["has_dialysate_48h"] == 0).sum()
print(f"  Patients with no dialysate in 0-48h: {n_no_dialysate_48}")

del crrt_rows, wide_df
import gc; gc.collect()

# ===================================================================
# STEP 5: Labs at 12h — last value in +3h to +12h window
# ===================================================================
print("Step 5: Labs at t=12 (last obs in +3h to +12h, baseline fallback) …")

lab_cols_wide = ["lab_lactate", "lab_bicarbonate", "lab_potassium"]
needed_labs = ["hospitalization_id", "event_dttm"] + [c for c in lab_cols_wide if c in available]

labs_wide = load_intermediate(INTERMEDIATE_DIR / "wide_df.parquet", columns=needed_labs)
labs_wide = labs_wide.merge(eb_map, on="hospitalization_id", how="inner")
labs_wide = labs_wide.merge(
    crrt_initiation[["encounter_block", "crrt_initiation_time"]],
    on="encounter_block", how="inner",
)
labs_wide["hours_from_crrt"] = (
    (labs_wide["event_dttm"] - labs_wide["crrt_initiation_time"]).dt.total_seconds() / 3600
)

labs_12h = labs_wide[
    (labs_wide["hours_from_crrt"] >= 3) & (labs_wide["hours_from_crrt"] <= 12)
].copy()

lab_12h_rename = {
    "lab_lactate": "lactate_12",
    "lab_bicarbonate": "bicarbonate_12",
    "lab_potassium": "potassium_12",
}
baseline_fallback = {
    "lactate_12": "lactate_0",
    "bicarbonate_12": "bicarbonate_0",
    "potassium_12": "potassium_0",
}

for col, new_name in lab_12h_rename.items():
    if col not in labs_12h.columns:
        result[new_name] = np.nan
        print(f"  {col} not in wide_df, filling {new_name} with NaN")
        continue
    lab_vals = (
        labs_12h[labs_12h[col].notna()]
        .sort_values("event_dttm")
        .groupby("encounter_block")[col]
        .last()
        .reset_index()
        .rename(columns={col: new_name})
    )
    result = result.merge(lab_vals, on="encounter_block", how="left")
    fallback_col = baseline_fallback[new_name]
    if fallback_col in result.columns:
        n_before = result[new_name].notna().sum()
        result[new_name] = result[new_name].fillna(result[fallback_col])
        n_after = result[new_name].notna().sum()
        print(f"  {new_name}: {n_before}/{len(result)} from +3-12h window, "
              f"{n_after - n_before} filled from baseline, {n_after}/{len(result)} total")

# 5b: Labs at 24h (sensitivity) — last obs in +12h to +24h, baseline fallback
print("Step 5b: Labs at t=24 (last obs in +12h to +24h, baseline fallback) …")

labs_24h = labs_wide[
    (labs_wide["hours_from_crrt"] >= 12) & (labs_wide["hours_from_crrt"] <= 24)
].copy()

lab_24h_rename = {
    "lab_lactate": "lactate_24",
    "lab_bicarbonate": "bicarbonate_24",
    "lab_potassium": "potassium_24",
}
baseline_fallback_24 = {
    "lactate_24": "lactate_0",
    "bicarbonate_24": "bicarbonate_0",
    "potassium_24": "potassium_0",
}

for col, new_name in lab_24h_rename.items():
    if col not in labs_24h.columns:
        result[new_name] = np.nan
        print(f"  {col} not in wide_df, filling {new_name} with NaN")
        continue
    lab_vals = (
        labs_24h[labs_24h[col].notna()]
        .sort_values("event_dttm")
        .groupby("encounter_block")[col]
        .last()
        .reset_index()
        .rename(columns={col: new_name})
    )
    result = result.merge(lab_vals, on="encounter_block", how="left")
    fallback_col = baseline_fallback_24[new_name]
    if fallback_col in result.columns:
        n_before = result[new_name].notna().sum()
        result[new_name] = result[new_name].fillna(result[fallback_col])
        n_after = result[new_name].notna().sum()
        print(f"  {new_name}: {n_before}/{len(result)} from +12-24h window, "
              f"{n_after - n_before} filled from baseline, {n_after}/{len(result)} total")

del labs_wide, labs_12h, labs_24h
gc.collect()

# ===================================================================
# STEP 6: SOFA at 12h
# ===================================================================
print("Step 6: Computing SOFA at t=12 (0h to +12h window) …")
sofa_cohort_12h = crrt_initiation.merge(eb_map, on="encounter_block")
sofa_cohort_12h["start_dttm"] = sofa_cohort_12h["crrt_initiation_time"]
sofa_cohort_12h["end_dttm"] = sofa_cohort_12h["crrt_initiation_time"] + pd.Timedelta(hours=12)
sofa_cohort_12h = sofa_cohort_12h[["hospitalization_id", "encounter_block", "start_dttm", "end_dttm"]]

sofa_12h_scores = compute_sofa_polars(
    data_directory=TABLES_PATH,
    cohort_df=pl.from_pandas(sofa_cohort_12h),
    filetype=FILE_TYPE,
    id_name="encounter_block",
    extremal_type="worst",
    fill_na_scores_with_zero=False,
    remove_outliers=True,
    timezone=TIMEZONE,
).to_pandas()

sofa_12h = sofa_12h_scores[["encounter_block", "sofa_total"]].rename(
    columns={"sofa_total": "sofa_total_12"}
)
result = result.merge(sofa_12h, on="encounter_block", how="left")
n_sofa = result["sofa_total_12"].notna().sum()
print(f"  SOFA at 12h: {n_sofa}/{len(result)} encounters")

gc.collect()

# ===================================================================
# STEP 6b: SOFA at 24h (sensitivity)
# ===================================================================
print("Step 6b: Computing SOFA at t=24 (0h to +24h window) …")
sofa_cohort_24h = crrt_initiation.merge(eb_map, on="encounter_block")
sofa_cohort_24h["start_dttm"] = sofa_cohort_24h["crrt_initiation_time"]
sofa_cohort_24h["end_dttm"] = sofa_cohort_24h["crrt_initiation_time"] + pd.Timedelta(hours=24)
sofa_cohort_24h = sofa_cohort_24h[["hospitalization_id", "encounter_block", "start_dttm", "end_dttm"]]

sofa_24h_scores = compute_sofa_polars(
    data_directory=TABLES_PATH,
    cohort_df=pl.from_pandas(sofa_cohort_24h),
    filetype=FILE_TYPE,
    id_name="encounter_block",
    extremal_type="worst",
    fill_na_scores_with_zero=False,
    remove_outliers=True,
    timezone=TIMEZONE,
).to_pandas()

sofa_24h = sofa_24h_scores[["encounter_block", "sofa_total"]].rename(
    columns={"sofa_total": "sofa_total_24"}
)
result = result.merge(sofa_24h, on="encounter_block", how="left")
n_sofa_24 = result["sofa_total_24"].notna().sum()
print(f"  SOFA at 24h: {n_sofa_24}/{len(result)} encounters")

gc.collect()

# ===================================================================
# STEP 7: Oxygenation index, NEE, IMV status at t=0 and t=12
# ===================================================================
print("Step 7: Oxygenation index, NEE, IMV status …")

# Columns needed from wide_df
oxy_nee_cols = [
    "lab_po2_arterial", "vital_spo2", "resp_fio2_set",
    "med_cont_nee", "resp_device_category", "resp_lpm_set",
]
needed_extra = ["hospitalization_id", "event_dttm"] + [c for c in oxy_nee_cols if c in available]
extra_df = load_intermediate(INTERMEDIATE_DIR / "wide_df.parquet", columns=needed_extra)
extra_df = extra_df.merge(eb_map, on="hospitalization_id", how="inner")
extra_df = extra_df.merge(
    crrt_initiation[["encounter_block", "crrt_initiation_time"]],
    on="encounter_block", how="inner",
)
extra_df["hours_from_crrt"] = (
    (extra_df["event_dttm"] - extra_df["crrt_initiation_time"]).dt.total_seconds() / 3600
)

# --- FiO2 imputation: LPM→FiO2 for nasal cannula, room air default ---
if "resp_fio2_set" in extra_df.columns:
    # Nasal cannula LPM→FiO2 conversion where FiO2 is missing
    if "resp_lpm_set" in extra_df.columns and "resp_device_category" in extra_df.columns:
        lpm_mask = (
            extra_df["resp_fio2_set"].isna()
            & extra_df["resp_lpm_set"].notna()
            & (extra_df["resp_lpm_set"] >= 1)
            & (extra_df["resp_lpm_set"] <= 10)
            & extra_df["resp_device_category"].str.lower().str.contains("nasal", na=False)
        )
        lpm_rounded = extra_df.loc[lpm_mask, "resp_lpm_set"].round().astype(int)
        lpm_to_fio2 = {1: 0.24, 2: 0.28, 3: 0.32, 4: 0.36, 5: 0.40,
                       6: 0.44, 7: 0.48, 8: 0.52, 9: 0.56, 10: 0.60}
        extra_df.loc[lpm_mask, "resp_fio2_set"] = lpm_rounded.map(lpm_to_fio2)
        n_lpm = lpm_mask.sum()
        print(f"  FiO2 imputed from LPM (nasal cannula): {n_lpm} rows")

    # Room air default — if device is room air and FiO2 still missing, set 0.21
    if "resp_device_category" in extra_df.columns:
        room_air_mask = (
            extra_df["resp_fio2_set"].isna()
            & extra_df["resp_device_category"].str.lower().str.contains("room air", na=False)
        )
        extra_df.loc[room_air_mask, "resp_fio2_set"] = 0.21
        n_ra = room_air_mask.sum()
        print(f"  FiO2 set to 0.21 for room air: {n_ra} rows")

    # Forward-fill FiO2 within each patient (max 8 rows ≈ ~4-8h in dense ICU data)
    n_fio2_before = extra_df["resp_fio2_set"].notna().sum()
    extra_df = extra_df.sort_values(["encounter_block", "event_dttm"])
    extra_df["resp_fio2_set"] = extra_df.groupby("encounter_block")["resp_fio2_set"].ffill(limit=8)
    n_fio2_after = extra_df["resp_fio2_set"].notna().sum()
    print(f"  FiO2 forward-filled: {n_fio2_after - n_fio2_before} additional rows")

# --- Oxygenation index: P/F if PaO2 available, else S/F = SpO2 / (FiO2/100) ---
if "resp_fio2_set" in extra_df.columns:
    fio2_frac = extra_df["resp_fio2_set"].copy()
    # If FiO2 > 1, assume it's in %, convert to fraction
    fio2_frac = np.where(fio2_frac > 1, fio2_frac / 100.0, fio2_frac)
    fio2_frac = np.where((fio2_frac > 0) & (fio2_frac <= 1), fio2_frac, np.nan)
else:
    fio2_frac = np.full(len(extra_df), np.nan)

has_pao2 = "lab_po2_arterial" in extra_df.columns
has_spo2 = "vital_spo2" in extra_df.columns

if has_pao2:
    pf_ratio = extra_df["lab_po2_arterial"] / fio2_frac
else:
    pf_ratio = pd.Series(np.nan, index=extra_df.index)

if has_spo2:
    sf_ratio = extra_df["vital_spo2"] / fio2_frac
else:
    sf_ratio = pd.Series(np.nan, index=extra_df.index)

# Prefer P/F, fallback to S/F
extra_df["oxygenation_index"] = np.where(pf_ratio.notna(), pf_ratio, sf_ratio)

# --- IMV status: 1 if device_category contains 'imv' (case-insensitive) ---
if "resp_device_category" in extra_df.columns:
    extra_df["imv_status"] = (
        extra_df["resp_device_category"]
        .str.lower()
        .str.strip()
        .eq("imv")
        .astype("Int32")
    )
else:
    extra_df["imv_status"] = pd.array([pd.NA] * len(extra_df), dtype="Int32")

# --- Extract t=0 (-12h to +3h) and t=12 (+3h to +12h) for each variable ---
def extract_last_in_window(df, col, h_low, h_high):
    """Last non-null value per encounter_block in the time window."""
    mask = (
        (df["hours_from_crrt"] >= h_low)
        & (df["hours_from_crrt"] <= h_high)
        & df[col].notna()
    )
    return (
        df.loc[mask]
        .sort_values("event_dttm")
        .groupby("encounter_block")[col]
        .last()
        .reset_index()
    )

for var, col in [
    ("oxygenation_index", "oxygenation_index"),
    ("norepinephrine_equivalent", "med_cont_nee"),
    ("imv_status", "imv_status"),
]:
    if col not in extra_df.columns:
        result[f"{var}_0"] = np.nan
        result[f"{var}_12"] = np.nan
        print(f"  {col} not available, {var}_0 and {var}_12 set to NaN")
        continue

    # t=0: last value in -12h to +3h
    t0 = extract_last_in_window(extra_df, col, -12, 3)
    t0 = t0.rename(columns={col: f"{var}_0"})
    result = result.merge(t0, on="encounter_block", how="left")

    # t=12: last value in +3h to +12h, fallback to t=0
    t12 = extract_last_in_window(extra_df, col, 3, 12)
    t12 = t12.rename(columns={col: f"{var}_12"})
    result = result.merge(t12, on="encounter_block", how="left")

    # Fallback t=12 → t=0
    result[f"{var}_12"] = result[f"{var}_12"].fillna(result[f"{var}_0"])

    n0 = result[f"{var}_0"].notna().sum()
    n12 = result[f"{var}_12"].notna().sum()
    print(f"  {var}_0: {n0}/{len(result)}, {var}_12: {n12}/{len(result)}")

# --- 7b. Extract t=24 values (sensitivity) ---
print("\n  Extracting t=24 values (sensitivity) …")
for var, col in [
    ("oxygenation_index", "oxygenation_index"),
    ("norepinephrine_equivalent", "med_cont_nee"),
    ("imv_status", "imv_status"),
]:
    if col not in extra_df.columns:
        result[f"{var}_24"] = np.nan
        print(f"  {col} not available, {var}_24 set to NaN")
        continue

    # t=24: last value in +12h to +24h, fallback to t=0
    t24 = extract_last_in_window(extra_df, col, 12, 24)
    t24 = t24.rename(columns={col: f"{var}_24"})
    result = result.merge(t24, on="encounter_block", how="left")

    # Fallback t=24 → t=0
    result[f"{var}_24"] = result[f"{var}_24"].fillna(result[f"{var}_0"])

    n24 = result[f"{var}_24"].notna().sum()
    print(f"  {var}_24: {n24}/{len(result)}")

# --- 7a. NEE recovery: impute 0 for patients not on vasopressors, widen window ---
print("\n  NEE imputation …")
if "med_cont_nee" in extra_df.columns:
    n_miss_before = result["norepinephrine_equivalent_0"].isna().sum()

    # Patients with ANY vasopressor row in wide_df
    ebs_with_nee = extra_df.loc[
        extra_df["med_cont_nee"].notna(), "encounter_block"
    ].unique()

    # Step 1: Patients with NO vasopressor rows → NEE = 0
    no_vasopressor_mask = (
        ~result["encounter_block"].isin(ebs_with_nee)
        & result["norepinephrine_equivalent_0"].isna()
    )
    result.loc[no_vasopressor_mask, "norepinephrine_equivalent_0"] = 0.0
    result.loc[no_vasopressor_mask, "norepinephrine_equivalent_12"] = 0.0
    n_no_pressor = no_vasopressor_mask.sum()
    print(f"    NEE=0 for {n_no_pressor} patients with no vasopressor records")

    # Step 2: Widen window to [-24h, +3h] for remaining missing
    still_missing = result["norepinephrine_equivalent_0"].isna()
    if still_missing.any():
        t0_wide = extract_last_in_window(extra_df, "med_cont_nee", -24, 3)
        t0_wide = t0_wide.rename(columns={"med_cont_nee": "nee_0_wide"})
        result = result.merge(t0_wide, on="encounter_block", how="left")
        result["norepinephrine_equivalent_0"] = result["norepinephrine_equivalent_0"].fillna(
            result["nee_0_wide"]
        )
        result.drop(columns=["nee_0_wide"], inplace=True)
        n_recovered_wide = n_miss_before - n_no_pressor - result["norepinephrine_equivalent_0"].isna().sum()
        print(f"    NEE recovered {n_recovered_wide} via widened [-24h, +3h] window")

    # Step 3: Any remaining → 0 (pressors likely stopped well before CRRT)
    still_missing = result["norepinephrine_equivalent_0"].isna()
    if still_missing.any():
        n_remaining = still_missing.sum()
        result.loc[still_missing, "norepinephrine_equivalent_0"] = 0.0
        print(f"    NEE=0 for {n_remaining} remaining (pressors likely stopped)")

    # Cascade to _12
    result["norepinephrine_equivalent_12"] = result["norepinephrine_equivalent_12"].fillna(
        result["norepinephrine_equivalent_0"]
    )

    n_miss_after = result["norepinephrine_equivalent_0"].isna().sum()
    print(f"    NEE_0 missing: {n_miss_before} → {n_miss_after}")

    # Cascade NEE imputation to _24 (sensitivity)
    result.loc[no_vasopressor_mask, "norepinephrine_equivalent_24"] = 0.0
    result["norepinephrine_equivalent_24"] = result["norepinephrine_equivalent_24"].fillna(
        result["norepinephrine_equivalent_0"]
    )
    n_nee_24_miss = result["norepinephrine_equivalent_24"].isna().sum()
    print(f"    NEE_24 missing after cascade: {n_nee_24_miss}")

# --- 7b. Oxygenation index recovery: widen window to [-24h, +3h] ---
print("\n  Oxygenation index imputation …")
if "oxygenation_index" in extra_df.columns:
    n_miss_before = result["oxygenation_index_0"].isna().sum()

    if n_miss_before > 0:
        t0_wide = extract_last_in_window(extra_df, "oxygenation_index", -24, 3)
        t0_wide = t0_wide.rename(columns={"oxygenation_index": "oxy_0_wide"})
        result = result.merge(t0_wide, on="encounter_block", how="left")
        result["oxygenation_index_0"] = result["oxygenation_index_0"].fillna(
            result["oxy_0_wide"]
        )
        result.drop(columns=["oxy_0_wide"], inplace=True)

        # Cascade to _12
        result["oxygenation_index_12"] = result["oxygenation_index_12"].fillna(
            result["oxygenation_index_0"]
        )

    n_miss_after = result["oxygenation_index_0"].isna().sum()
    print(f"    oxygenation_index_0 missing: {n_miss_before} → {n_miss_after}")

    # Cascade oxygenation imputation to _24 (sensitivity)
    result["oxygenation_index_24"] = result["oxygenation_index_24"].fillna(
        result["oxygenation_index_0"]
    )
    n_oxy_24_miss = result["oxygenation_index_24"].isna().sum()
    print(f"    oxygenation_index_24 missing after cascade: {n_oxy_24_miss}")

# --- 7c. Lactate: impute normal value (1.0 mmol/L) for missing ---
print("\n  Lactate imputation …")
n_miss_before = result["lactate_0"].isna().sum()
if n_miss_before > 0:
    result["lactate_0"] = result["lactate_0"].fillna(1.0)
    # Also cascade to lactate_12 if it was using baseline fallback
    result["lactate_12"] = result["lactate_12"].fillna(result["lactate_0"])
    # Also cascade to lactate_24 (sensitivity)
    result["lactate_24"] = result["lactate_24"].fillna(result["lactate_0"])
    print(f"    lactate_0 missing: {n_miss_before} → {result['lactate_0'].isna().sum()} (imputed normal=1.0)")
else:
    print(f"    lactate_0: no missing values")

del extra_df
gc.collect()

# ===================================================================
# STEP 8: CCI components (17 binary columns from clifpy)
# ===================================================================
print("Step 8: Computing CCI components from hospital_diagnosis (POA only) …")

import clifpy

# Load hospital_diagnosis table
diag_df = pd.read_parquet(
    Path(TABLES_PATH) / f"clif_hospital_diagnosis.{FILE_TYPE}"
) if FILE_TYPE == "parquet" else pd.read_csv(
    Path(TABLES_PATH) / f"clif_hospital_diagnosis.csv"
)

# Filter to our cohort and present-on-admission only
hosp_ids = set(eb_map["hospitalization_id"].unique())
diag_df = diag_df[
    (diag_df["hospitalization_id"].isin(hosp_ids))
    & (diag_df["poa_present"] == 1)
].copy()
print(f"  {len(diag_df)} POA diagnosis rows for cohort")

# Pre-clean diagnosis codes: remove dots so "E11.65" → "E1165"
# (clifpy internally splits on "." and takes only the part before it,
#  which loses sub-code precision needed for CCI prefix matching)
diag_df["diagnosis_code"] = (
    diag_df["diagnosis_code"]
    .astype(str)
    .str.replace(".", "", regex=False)
)

# Calculate CCI using clifpy (returns 17 condition columns + cci_score)
cci_result = clifpy.calculate_cci(diag_df, hierarchy=True)

if isinstance(cci_result, pl.DataFrame):
    cci_result = cci_result.to_pandas()

# Map hospitalization_id → encounter_block
cci_result["hospitalization_id"] = cci_result["hospitalization_id"].astype(str)
eb_map_str = eb_map.copy()
eb_map_str["hospitalization_id"] = eb_map_str["hospitalization_id"].astype(str)
cci_result = cci_result.merge(eb_map_str, on="hospitalization_id", how="inner")

# Add cci_ prefix to all 17 condition columns (keep names as clifpy returns them)
clifpy_condition_cols = [
    "myocardial_infarction", "congestive_heart_failure",
    "peripheral_vascular_disease", "cerebrovascular_disease",
    "dementia", "chronic_pulmonary_disease",
    "connective_tissue_disease", "peptic_ulcer_disease",
    "mild_liver_disease", "diabetes_uncomplicated",
    "diabetes_with_complications", "hemiplegia",
    "renal_disease", "cancer",
    "moderate_severe_liver_disease", "metastatic_solid_tumor",
    "aids",
]
cci_rename = {c: f"cci_{c}" for c in clifpy_condition_cols if c in cci_result.columns}
cci_result = cci_result.rename(columns=cci_rename)

# Final CCI columns = encounter_block + all cci_ prefixed condition cols
cci_cols = ["encounter_block"] + [f"cci_{c}" for c in clifpy_condition_cols
                                   if f"cci_{c}" in cci_result.columns]

cci_final = cci_result[cci_cols].copy()

# Ensure binary int
for c in cci_cols[1:]:
    cci_final[c] = cci_final[c].fillna(0).astype(int)

result = result.merge(cci_final, on="encounter_block", how="left")

# Fill CCI = 0 for encounters with no POA diagnosis records
for c in cci_cols[1:]:
    result[c] = result[c].fillna(0).astype(int)

n_any_cci = (result[cci_cols[1:]].sum(axis=1) > 0).sum()
print(f"  CCI: {n_any_cci}/{len(result)} encounters have at least one CCI condition")
print(f"  CCI columns: {cci_cols[1:]}")

del diag_df, cci_result
gc.collect()

# ===================================================================
# STEP 9: Outcome encoding (competing risk, 30-day censoring)
# ===================================================================
print("Step 9: Outcome encoding …")

outcomes = outcomes_df[["encounter_block", "discharge_category", "death_dttm",
                         "last_vital_dttm", "final_outcome_dttm"]].copy()
outcomes = outcomes.merge(
    crrt_initiation[["encounter_block", "crrt_initiation_time"]],
    on="encounter_block", how="left",
)

# If expired but no death_dttm, use last_vital_dttm as proxy
mask_expired = outcomes["discharge_category"].isin(["expired", "hospice"])
mask_no_death_dttm = mask_expired & outcomes["death_dttm"].isna()
outcomes.loc[mask_no_death_dttm, "death_dttm"] = outcomes.loc[mask_no_death_dttm, "last_vital_dttm"]

# Time from CRRT initiation to event (days)
outcomes["days_to_death"] = (
    (outcomes["death_dttm"] - outcomes["crrt_initiation_time"]).dt.total_seconds() / 86400
)
outcomes["discharge_dttm"] = outcomes["final_outcome_dttm"].fillna(outcomes["last_vital_dttm"])
outcomes["days_to_discharge"] = (
    (outcomes["discharge_dttm"] - outcomes["crrt_initiation_time"]).dt.total_seconds() / 86400
)

# Competing risk: 2=died, 1=discharged alive, 0=censored at 30d
# Column names use "90d" to match R script expectations
CENSOR_DAYS = 30

outcomes["outcome"] = 0
outcomes["time_to_event_90d"] = float(CENSOR_DAYS)

# Died within window
died_mask = outcomes["days_to_death"].notna() & (outcomes["days_to_death"] <= CENSOR_DAYS)
outcomes.loc[died_mask, "outcome"] = 2
outcomes.loc[died_mask, "time_to_event_90d"] = outcomes.loc[died_mask, "days_to_death"]

# Discharged alive within window
discharged_alive_mask = (
    ~died_mask
    & outcomes["days_to_discharge"].notna()
    & (outcomes["days_to_discharge"] <= CENSOR_DAYS)
    & ~mask_expired
)
outcomes.loc[discharged_alive_mask, "outcome"] = 1
outcomes.loc[discharged_alive_mask, "time_to_event_90d"] = outcomes.loc[
    discharged_alive_mask, "days_to_discharge"
]

outcomes["censored_at_90d"] = (outcomes["outcome"] == 0).astype(int)
outcomes["time_to_event_90d"] = outcomes["time_to_event_90d"].clip(lower=0)

# Also create a clearly-named 30-day version
outcomes["time_to_event_30d"] = outcomes["time_to_event_90d"].copy()

result = result.merge(
    outcomes[["encounter_block", "time_to_event_90d", "time_to_event_30d", "outcome", "censored_at_90d"]],
    on="encounter_block", how="left",
)

print(f"  Outcome distribution:")
print(f"    0 (censored):          {(result['outcome'] == 0).sum()}")
print(f"    1 (discharged alive):  {(result['outcome'] == 1).sum()}")
print(f"    2 (died):              {(result['outcome'] == 2).sum()}")

# ===================================================================
# STEP 9b: Exclude patients who died or were off CRRT within 24h
# ===================================================================
EXCLUDE_SHORT_CRRT = True  # Set False to keep all patients

if EXCLUDE_SHORT_CRRT:
    short_mask = (
        ((result["outcome"] == 2) & (result["time_to_event_90d"] <= 1.0))  # died ≤24h
        | (result["crrt_duration_days"] < 1.0)  # CRRT lasted <24h
    )
    n_excluded_short = short_mask.sum()
    result = result[~short_mask].reset_index(drop=True)
    print(f"\nStep 9b: Excluded {n_excluded_short} patients (died or off CRRT within 24h)")
    print(f"  Remaining: {len(result)}")

    # Exclude SCUF-only patients (no dialysate clearance in first 24h)
    scuf_no_clearance = (result["has_dialysate_24h"] == 0)
    n_excluded_scuf = scuf_no_clearance.sum()
    result = result[~scuf_no_clearance].reset_index(drop=True)
    print(f"  Excluded {n_excluded_scuf} SCUF-only patients (no dialysate in 0-24h)")
    print(f"  Remaining: {len(result)}")

# ===================================================================
# STEP 9c: Eligibility flag for 48h sensitivity analysis
# ===================================================================
print("\nStep 9c: Computing eligible_48h_sensitivity flag …")
eligible_48h_mask = (
    # Survived at least 48h (not dead within 48h)
    ~((result["outcome"] == 2) & (result["time_to_event_90d"] <= 2.0))
    # CRRT duration >= 48h (2 days)
    & (result["crrt_duration_days"] >= 2.0)
    # Has dialysate flow in 0-48h window
    & (result["has_dialysate_48h"] == 1)
)
result["eligible_48h_sensitivity"] = eligible_48h_mask.astype(int)
n_eligible = result["eligible_48h_sensitivity"].sum()
n_total = len(result)
print(f"  Eligible for 48h sensitivity: {n_eligible}/{n_total} "
      f"({n_eligible / n_total * 100:.1f}%)")
if n_eligible < 200:
    print(f"  WARNING: N={n_eligible} may be underpowered for sensitivity MSM")

# ===================================================================
# STEP 10: Final column ordering per schema (59 columns)
# ===================================================================
# Drop helper columns used for exclusion
for helper_col in ["has_dialysate_24h", "has_dialysate_48h"]:
    if helper_col in result.columns:
        result.drop(columns=[helper_col], inplace=True)

print("Step 10: Final column ordering …")
final_cols = [
    "encounter_block",
    "age_at_admission",
    "sex_category",
    "race_category",
    "ethnicity_category",
    "weight_kg",
    "crrt_mode_category",
    "crrt_dose_ml_kg_hr_0",
    "crrt_dose_0_12",
    "crrt_dose_12_24",
    "crrt_duration_days",
    "imv_duration_days",
    "creatinine_0",
    "lactate_0",
    "bicarbonate_0",
    "potassium_0",
    "sofa_total_0",
    "lactate_12",
    "bicarbonate_12",
    "potassium_12",
    "sofa_total_12",
    "time_to_event_90d",
    "time_to_event_30d",
    "outcome",
    "censored_at_90d",
    "oxygenation_index_0",
    "norepinephrine_equivalent_0",
    "imv_status_0",
    "oxygenation_index_12",
    "norepinephrine_equivalent_12",
    "imv_status_12",
    "cci_myocardial_infarction",
    "cci_congestive_heart_failure",
    "cci_peripheral_vascular_disease",
    "cci_cerebrovascular_disease",
    "cci_dementia",
    "cci_chronic_pulmonary_disease",
    "cci_connective_tissue_disease",
    "cci_peptic_ulcer_disease",
    "cci_mild_liver_disease",
    "cci_diabetes_uncomplicated",
    "cci_diabetes_with_complications",
    "cci_hemiplegia",
    "cci_renal_disease",
    "cci_cancer",
    "cci_moderate_severe_liver_disease",
    "cci_metastatic_solid_tumor",
    "cci_aids",
    # --- Sensitivity analysis columns (24h/48h windows) ---
    "crrt_dose_0_24",
    "crrt_dose_24_48",
    "lactate_24",
    "bicarbonate_24",
    "potassium_24",
    "sofa_total_24",
    "oxygenation_index_24",
    "norepinephrine_equivalent_24",
    "imv_status_24",
    "eligible_48h_sensitivity",
]

# Verify all columns present
missing = [c for c in final_cols if c not in result.columns]
if missing:
    print(f"  WARNING: Missing columns: {missing}")
    for c in missing:
        result[c] = np.nan

result = result[final_cols]

# --- Patient-level missingness diagnostic ---
covariate_cols = [c for c in final_cols if c != "encounter_block"]
diag_missing = result[covariate_cols].isna().astype(int)
diag_missing.insert(0, "encounter_block", result["encounter_block"].values)
diag_missing["total_missing"] = diag_missing.drop(columns=["encounter_block"]).sum(axis=1)
diag_path = OUTPUT_DIR / "missingness_diagnostic.csv"
diag_missing.to_csv(diag_path, index=False)
print(f"\nSaved patient-level missingness: {diag_path}")
print(f"  Patients with any missing covariate: {(diag_missing['total_missing'] > 0).sum()}/{len(diag_missing)}")

# --- Aggregate missingness summary (safe to share) ---
FINAL_DIR = project_root / "output" / "final" / "crrt_epi"
FINAL_DIR.mkdir(parents=True, exist_ok=True)
GRAPHS_DIR = FINAL_DIR / "graphs"
GRAPHS_DIR.mkdir(parents=True, exist_ok=True)
FINAL_DIR.mkdir(parents=True, exist_ok=True)

n_total = len(result)
agg_rows = []
for col in covariate_cols:
    n_miss = result[col].isna().sum()
    agg_rows.append({
        "variable": col,
        "n_total": n_total,
        "n_missing": n_miss,
        "n_present": n_total - n_miss,
        "pct_missing": round(n_miss / n_total * 100, 2),
    })
agg_df = pd.DataFrame(agg_rows).sort_values("pct_missing", ascending=False)
agg_csv_path = FINAL_DIR / "missingness_summary.csv"
agg_df.to_csv(agg_csv_path, index=False)
print(f"Saved aggregate missingness: {agg_csv_path}")

# --- Missingness heatmap ---
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Only show variables with >0% missing for a cleaner heatmap
agg_any = agg_df[agg_df["pct_missing"] > 0].copy()
if len(agg_any) > 0:
    fig, ax = plt.subplots(figsize=(8, max(3, len(agg_any) * 0.45)))
    bars = ax.barh(agg_any["variable"], agg_any["pct_missing"], color="#e74c3c", edgecolor="white")
    ax.set_xlabel("% Missing")
    ax.set_title(f"Covariate Missingness (n={n_total})")
    ax.invert_yaxis()
    for bar, pct, n in zip(bars, agg_any["pct_missing"], agg_any["n_missing"]):
        ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height() / 2,
                f"{pct}% (n={n})", va="center", fontsize=9)
    ax.set_xlim(0, agg_any["pct_missing"].max() * 1.4)
    plt.tight_layout()
    heatmap_path = GRAPHS_DIR / "missingness_heatmap.png"
    fig.savefig(heatmap_path, dpi=150)
    plt.close(fig)
    print(f"Saved missingness heatmap: {heatmap_path}")
else:
    print("  No missing covariates — heatmap skipped")

# Save
out_path = OUTPUT_DIR / "msm_competing_risk_df.parquet"
result.to_parquet(out_path, index=False)
print(f"\nSaved: {out_path} ({len(result)} rows x {len(result.columns)} cols)")

csv_path = OUTPUT_DIR / "msm_competing_risk_df.csv"
result.to_csv(csv_path, index=False)
print(f"Saved: {csv_path}")

# Summary
print("\n" + "=" * 60)
print("Column summary:")
print("=" * 60)
for col in final_cols:
    dtype = result[col].dtype
    n_miss = result[col].isna().sum()
    pct_miss = n_miss / len(result) * 100
    if result[col].dtype in ["float64", "int64", "int32", "Int32", "Int64"]:
        desc = f"median={result[col].median():.2f}" if n_miss < len(result) else "all NaN"
    else:
        desc = f"{result[col].nunique()} unique"
    print(f"  {col:40s}  {str(dtype):12s}  miss={n_miss:4d} ({pct_miss:4.1f}%)  {desc}")

print("\nDone!")
