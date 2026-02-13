"""
Build the competing-risk analysis dataframe per competing_risk_final_schema.xlsx.

Produces one row per encounter_block with:
  - Demographics and baseline characteristics
  - CRRT dose at initiation and mean dose for 0-12h / 12-24h intervals
  - Baseline labs (t=0): last value in -12h to +3h window (from tableone)
  - 12h labs (t=12): last value in +3h to +12h window, with baseline fallback
  - SOFA at t=0 and t=12
  - Competing-risk outcome (0=censored, 1=discharged, 2=died) capped at 30 days

Usage: uv run python code/build_competing_risk_df.py
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
outcomes_df = pd.read_parquet(INTERMEDIATE_DIR / "outcomes_df.parquet")
index_crrt_df = pd.read_parquet(INTERMEDIATE_DIR / "index_crrt_df.parquet")
crrt_initiation = pd.read_parquet(INTERMEDIATE_DIR / "crrt_initiation.parquet")
tableone_df = pd.read_parquet(INTERMEDIATE_DIR / "tableone_analysis_df.parquet")

eb_map = index_crrt_df[["hospitalization_id", "encounter_block"]].copy()
print(f"  {len(outcomes_df)} encounters")

# ===================================================================
# STEP 2: Demographics + static columns (from outcomes_df / index_crrt_df)
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
    "imv_duration_days": "imv_duration_days",
})
print(f"  {len(result)} rows after merge")

# ===================================================================
# STEP 3: Baseline labs (t=0) from tableone_analysis_df
#         Already computed with -12h to +3h window, last recorded value
# ===================================================================
print("Step 3: Baseline labs (t=0) from tableone …")
baseline_lab_map = {
    "creatinine_baseline": "creatinine_0",
    "lactate_baseline": "lactate_0",
    "bicarbonate_baseline": "bicarbonate_0",
    "potassium_baseline": "potassium_0",
}
baseline_cols = ["encounter_block"] + list(baseline_lab_map.keys())
baseline_labs = tableone_df[baseline_cols].rename(columns=baseline_lab_map)
result = result.merge(baseline_labs, on="encounter_block", how="left")

# SOFA at t=0 (already in tableone)
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

# Load only needed columns
crrt_cols = [
    "crrt_blood_flow_rate", "crrt_dialysate_flow_rate",
    "crrt_pre_filter_replacement_fluid_rate",
    "crrt_post_filter_replacement_fluid_rate",
    "crrt_mode_category",
]
available = {f.name for f in pq.read_schema(INTERMEDIATE_DIR / "wide_df.parquet")}
needed = ["hospitalization_id", "event_dttm"] + [c for c in crrt_cols if c in available]

wide_df = pd.read_parquet(INTERMEDIATE_DIR / "wide_df.parquet", columns=needed)
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

# Mode-specific dose formula (same as 00_cohort)
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

# 4a: Dose at initiation (first valid CRRT obs after t=0)
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
print(f"  Dose at initiation: {init_dose['crrt_dose_ml_kg_hr_0'].notna().sum()} encounters")
print(f"  Dose 0-12h: {dose_0_12['crrt_dose_0_12'].notna().sum()} encounters")
print(f"  Dose 12-24h: {dose_12_24['crrt_dose_12_24'].notna().sum()} encounters")

# Free CRRT-specific memory
del crrt_rows, wide_df

# ===================================================================
# STEP 5: Labs at 12h — last value in +3h to +12h window, baseline fallback
# ===================================================================
print("Step 5: Labs at t=12 (last obs in +3h to +12h, baseline fallback) …")

lab_cols_wide = ["lab_lactate", "lab_bicarbonate", "lab_potassium"]
needed_labs = ["hospitalization_id", "event_dttm"] + [c for c in lab_cols_wide if c in available]

labs_wide = pd.read_parquet(INTERMEDIATE_DIR / "wide_df.parquet", columns=needed_labs)
labs_wide = labs_wide.merge(eb_map, on="hospitalization_id", how="inner")
labs_wide = labs_wide.merge(
    crrt_initiation[["encounter_block", "crrt_initiation_time"]],
    on="encounter_block", how="inner",
)
labs_wide["hours_from_crrt"] = (
    (labs_wide["event_dttm"] - labs_wide["crrt_initiation_time"]).dt.total_seconds() / 3600
)

# Filter to +3h to +12h window (avoids overlap with t=0 window which is -12h to +3h)
labs_12h = labs_wide[
    (labs_wide["hours_from_crrt"] >= 3) & (labs_wide["hours_from_crrt"] <= 12)
].copy()

# For each lab, take the last non-null value per encounter in the window
# Fall back to baseline t=0 value from tableone if still missing
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
    # Fallback: fill missing t=12 values with baseline t=0
    fallback_col = baseline_fallback[new_name]
    if fallback_col in result.columns:
        n_before = result[new_name].notna().sum()
        result[new_name] = result[new_name].fillna(result[fallback_col])
        n_after = result[new_name].notna().sum()
        print(f"  {new_name}: {n_before}/{len(result)} from +3-12h window, "
              f"{n_after - n_before} filled from baseline, {n_after}/{len(result)} total")
    else:
        n_avail = result[new_name].notna().sum()
        print(f"  {new_name}: {n_avail}/{len(result)} encounters have data")

del labs_wide, labs_12h

# ===================================================================
# STEP 6: SOFA at 12h — compute from raw data using 0-12h window
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

# ===================================================================
# STEP 7: Outcome encoding (competing risk)
# ===================================================================
print("Step 7: Outcome encoding …")

# Merge crrt_initiation_time + death/discharge info
outcomes = outcomes_df[["encounter_block", "discharge_category", "death_dttm",
                         "last_vital_dttm", "final_outcome_dttm"]].copy()
outcomes = outcomes.merge(
    crrt_initiation[["encounter_block", "crrt_initiation_time"]],
    on="encounter_block", how="left",
)

# If expired but no death_dttm, use last_vital_dttm as proxy (same as tableone logic)
mask_expired = outcomes["discharge_category"].isin(["expired", "hospice"])
mask_no_death_dttm = mask_expired & outcomes["death_dttm"].isna()
outcomes.loc[mask_no_death_dttm, "death_dttm"] = outcomes.loc[mask_no_death_dttm, "last_vital_dttm"]

# Time from CRRT initiation to event (days)
outcomes["days_to_death"] = (
    (outcomes["death_dttm"] - outcomes["crrt_initiation_time"]).dt.total_seconds() / 86400
)
# Use last_vital_dttm as discharge proxy when final_outcome_dttm is missing
outcomes["discharge_dttm"] = outcomes["final_outcome_dttm"].fillna(outcomes["last_vital_dttm"])
outcomes["days_to_discharge"] = (
    (outcomes["discharge_dttm"] - outcomes["crrt_initiation_time"]).dt.total_seconds() / 86400
)

# Competing risk outcome: 2=died, 1=discharged alive, 0=censored at 30d
# Priority: death within 30d > discharge within 30d > censored
outcomes["outcome"] = 0  # default: censored
outcomes["time_to_event_30d"] = 30.0  # default: administrative censoring

# Died within 30 days
died_mask = outcomes["days_to_death"].notna() & (outcomes["days_to_death"] <= 30)
outcomes.loc[died_mask, "outcome"] = 2
outcomes.loc[died_mask, "time_to_event_30d"] = outcomes.loc[died_mask, "days_to_death"]

# Discharged alive within 30 days (not already marked as died)
discharged_alive_mask = (
    ~died_mask
    & outcomes["days_to_discharge"].notna()
    & (outcomes["days_to_discharge"] <= 30)
    & ~mask_expired
)
outcomes.loc[discharged_alive_mask, "outcome"] = 1
outcomes.loc[discharged_alive_mask, "time_to_event_30d"] = outcomes.loc[
    discharged_alive_mask, "days_to_discharge"
]

# Censored at 30 days flag
outcomes["censored_at_30d"] = (outcomes["outcome"] == 0).astype(int)

# Ensure non-negative time
outcomes["time_to_event_30d"] = outcomes["time_to_event_30d"].clip(lower=0)

result = result.merge(
    outcomes[["encounter_block", "time_to_event_30d", "outcome", "censored_at_30d"]],
    on="encounter_block", how="left",
)

print(f"  Outcome distribution:")
print(f"    0 (censored):          {(result['outcome'] == 0).sum()}")
print(f"    1 (discharged alive):  {(result['outcome'] == 1).sum()}")
print(f"    2 (died):              {(result['outcome'] == 2).sum()}")

# ===================================================================
# STEP 8: Final column ordering per schema
# ===================================================================
print("Step 8: Final column ordering …")
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
    "time_to_event_30d",
    "outcome",
    "censored_at_30d",
]

# Verify all columns present
missing = [c for c in final_cols if c not in result.columns]
if missing:
    print(f"  WARNING: Missing columns: {missing}")

result = result[final_cols]

# Save
out_path = OUTPUT_DIR / "competing_risk_df.parquet"
result.to_parquet(out_path, index=False)
print(f"\nSaved: {out_path} ({len(result)} rows x {len(result.columns)} cols)")

# Also save CSV for quick inspection
csv_path = OUTPUT_DIR / "competing_risk_df.csv"
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
    print(f"  {col:25s}  {str(dtype):12s}  miss={n_miss:4d} ({pct_miss:4.1f}%)  {desc}")

print("\nDone!")
