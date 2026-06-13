"""
Construct CRRT Table 1 matching reference/tableone_crrt_template.csv.

Produces a Table 1 with 5 column groups:
  - Baseline: characteristics at CRRT initiation (-12h to +3h)
  - At 72h - Survivors: patients alive at 72h post-CRRT
  - At 72h - Non-survivors: patients dead within 72h of CRRT
  - At discharge - Survivors: patients surviving to hospital discharge (90-day cap)
  - At discharge - Non-survivors: patients who died in-hospital (90-day cap)

Continuous variables use last recorded value in the relevant window (not mean).
Labs use unlimited forward-fill (no 12h cap).
SOFA scores use window-specific computations (baseline, 72h, post-CRRT).

Usage: uv run python code/construct_crrt_tableone.py
"""

import gc
import sys
from pathlib import Path

import pandas as pd
import polars as pl
import pyarrow.parquet as pq

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root / "code"))

# Load config directly (avoid importing utils package which conflicts with clifpy.utils)
import json  # noqa: E402

with open(project_root / "config" / "config.json") as _f:
    config = json.load(_f)

from pipeline_helpers import (  # noqa: E402
    validate_config, load_intermediate, get_tables_path, course_average_intensity,
)
config = validate_config(config)

from sofa_calculator import compute_sofa_polars  # noqa: E402

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
INTERMEDIATE_DIR = project_root / "output" / "intermediate"
FINAL_DIR = project_root / "output" / "final" / "crrt_epi"
FINAL_DIR.mkdir(parents=True, exist_ok=True)
FINAL_DIR.mkdir(parents=True, exist_ok=True)

SITE_NAME = config["site_name"]
TABLES_PATH = config["tables_path"]
FILE_TYPE = config["file_type"]
TIMEZONE = config["timezone"]
HAS_CRRT_SETTINGS = config.get("has_crrt_settings", True)

# Column families (names WITHOUT their wide_df prefix)
LAB_COLS = [
    "bicarbonate", "bun", "calcium_total", "chloride", "creatinine",
    "magnesium", "glucose_serum", "lactate", "potassium", "sodium",
    "ph_arterial", "po2_arterial", "phosphate",
]
VASO_COLS = [
    "angiotensin", "dobutamine", "dopamine", "epinephrine",
    "norepinephrine", "phenylephrine", "vasopressin",
]
NEE_COLS = ["nee"]
RESP_SETTING_COLS = [
    "fio2_set", "peep_set", "resp_rate_set", "tidal_volume_set",
    "pressure_control_set", "pressure_support_set",
    "peak_inspiratory_pressure_set",
]
CRRT_SETTING_COLS = [
    "blood_flow_rate", "pre_filter_replacement_fluid_rate",
    "post_filter_replacement_fluid_rate", "dialysate_flow_rate",
    "ultrafiltration_out",
]
SOFA_COLS = [
    "sofa_cv_97", "sofa_coag", "sofa_renal", "sofa_liver",
    "sofa_resp", "sofa_cns", "sofa_total",
]

# Columns whose values are patient-level (not encounter-level)
PATIENT_LEVEL_CATEGORICALS = {"sex_category", "race_category", "ethnicity_category"}
PATIENT_LEVEL_CONTINUOUS = {"age_at_admission"}


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------
def fmt_median_iqr(series: pd.Series) -> str:
    """Format as 'median [Q1, Q3] (n/N)'."""
    valid = series.dropna()
    n, N = len(valid), len(series)
    if n == 0:
        return "—"
    med = valid.median()
    q1 = valid.quantile(0.25)
    q3 = valid.quantile(0.75)
    return f"{med:.2f} [{q1:.2f}, {q3:.2f}] ({n}/{N})"


def fmt_n_pct(count: int, total: int) -> str:
    """Format as 'n (pct%)'."""
    if total == 0:
        return "0 (0.0%)"
    return f"{count} ({count / total * 100:.1f}%)"


# ===================================================================
# STEP 1: Load intermediate files
# ===================================================================
print("Step 1: Loading intermediate files …")
outcomes_df = load_intermediate(INTERMEDIATE_DIR / "outcomes_df.parquet")
cohort_df = load_intermediate(INTERMEDIATE_DIR / "cohort_df.parquet")
index_crrt_df = load_intermediate(INTERMEDIATE_DIR / "index_crrt_df.parquet")
crrt_initiation = load_intermediate(INTERMEDIATE_DIR / "crrt_initiation.parquet")

# N:1 mapping — all hospitalization_ids per encounter_block (needed for stitched sites)
eb_map = cohort_df[["hospitalization_id", "encounter_block"]].copy()
print(f"  {len(outcomes_df)} encounters loaded ({len(eb_map)} hospitalization_ids in eb_map)")

# Pre-filter CLIF tables to cohort (prevents OOM during SOFA/ASE)
cohort_hosp_ids = set(eb_map["hospitalization_id"].unique())
cohort_patient_ids = set(outcomes_df["patient_id"].unique()) if "patient_id" in outcomes_df.columns else None
TABLES_PATH = get_tables_path(config, cohort_hosp_ids, INTERMEDIATE_DIR, patient_ids=cohort_patient_ids)

# Merge crrt_initiation_time early (needed for 90-day mortality cap)
outcomes_df = outcomes_df.merge(
    crrt_initiation[["encounter_block", "crrt_initiation_time"]],
    on="encounter_block",
    how="left",
)

# Redefine in-hospital mortality within this script:
#   died = discharge_category == 'expired' or 'hospice'
#   If expired but no death_dttm, use last_vital_dttm as proxy for timing
#   Cap: deaths > 90 days after CRRT initiation are excluded
old_deaths = int(outcomes_df["in_hosp_death"].sum())
mask_expired = outcomes_df["discharge_category"].isin(["expired", "hospice"])
mask_no_death_dttm = mask_expired & outcomes_df["death_dttm"].isna()
outcomes_df.loc[mask_no_death_dttm, "death_dttm"] = outcomes_df.loc[mask_no_death_dttm, "last_vital_dttm"]
outcomes_df["in_hosp_death"] = mask_expired.astype(int)

# 90-day cap: if death occurred > 90 days after CRRT initiation, exclude
death_gap = (outcomes_df["death_dttm"] - outcomes_df["crrt_initiation_time"]).dt.total_seconds() / 86400
mask_beyond_90d = death_gap > 90
outcomes_df.loc[mask_beyond_90d, "in_hosp_death"] = 0
new_deaths = int(outcomes_df["in_hosp_death"].sum())
capped = int(mask_beyond_90d.sum())
print(f"  Mortality redefined: {old_deaths} → {new_deaths} in-hospital deaths ({capped} capped at 90d)")


# ===================================================================
# STEP 2: Compute SOFA scores
# ===================================================================
# 2a: SOFA at CRRT initiation (-12h to +3h) — used for Overall and Alive subgroups
print("Step 2a: Computing SOFA at initiation (-12h to +3h) …")
sofa_cohort = crrt_initiation.merge(eb_map, on="encounter_block")
sofa_cohort["start_dttm"] = sofa_cohort["crrt_initiation_time"] - pd.Timedelta(hours=12)
sofa_cohort["end_dttm"] = sofa_cohort["crrt_initiation_time"] + pd.Timedelta(hours=3)
sofa_cohort = sofa_cohort[["hospitalization_id", "encounter_block", "start_dttm", "end_dttm"]]

sofa_scores = compute_sofa_polars(
    data_directory=TABLES_PATH,
    cohort_df=pl.from_pandas(sofa_cohort),
    filetype=FILE_TYPE,
    id_name="encounter_block",
    extremal_type="worst",
    fill_na_scores_with_zero=False,
    remove_outliers=True,
    timezone=TIMEZONE,
).to_pandas()

sofa_keep = ["encounter_block"] + [c for c in SOFA_COLS if c in sofa_scores.columns]
sofa_df = sofa_scores[sofa_keep]
del sofa_scores, sofa_cohort
gc.collect()
print(f"  SOFA at initiation: {len(sofa_df)} encounters")

# 2b: SOFA 72h window (0h to 72h) — used for 72h survivor/non-survivor columns
print("Step 2b: Computing SOFA over 72h window (0h to +72h) …")
sofa_cohort_72h = crrt_initiation.merge(eb_map, on="encounter_block")
sofa_cohort_72h["start_dttm"] = sofa_cohort_72h["crrt_initiation_time"]
sofa_cohort_72h["end_dttm"] = sofa_cohort_72h["crrt_initiation_time"] + pd.Timedelta(hours=72)
sofa_cohort_72h = sofa_cohort_72h[["hospitalization_id", "encounter_block", "start_dttm", "end_dttm"]]

sofa_72h_scores = compute_sofa_polars(
    data_directory=TABLES_PATH,
    cohort_df=pl.from_pandas(sofa_cohort_72h),
    filetype=FILE_TYPE,
    id_name="encounter_block",
    extremal_type="worst",
    fill_na_scores_with_zero=False,
    remove_outliers=True,
    timezone=TIMEZONE,
).to_pandas()

sofa_72h_keep = ["encounter_block"] + [c for c in SOFA_COLS if c in sofa_72h_scores.columns]
sofa_72h_df = sofa_72h_scores[sofa_72h_keep].rename(
    columns={c: f"{c}_72h" for c in SOFA_COLS if c in sofa_72h_scores.columns}
)
del sofa_72h_scores, sofa_cohort_72h
gc.collect()
print(f"  SOFA 72h: {len(sofa_72h_df)} encounters")

# 2c: SOFA post-CRRT (0h to last_vital_dttm) — used for discharge survivor/non-survivor columns
print("Step 2c: Computing SOFA post-CRRT (0h to last_vital_dttm) …")
sofa_cohort_post = crrt_initiation.merge(eb_map, on="encounter_block")
sofa_cohort_post = sofa_cohort_post.merge(
    outcomes_df[["encounter_block", "last_vital_dttm"]],
    on="encounter_block",
    how="left",
)
sofa_cohort_post["start_dttm"] = sofa_cohort_post["crrt_initiation_time"]
sofa_cohort_post["end_dttm"] = sofa_cohort_post["last_vital_dttm"]
sofa_cohort_post = sofa_cohort_post[["hospitalization_id", "encounter_block", "start_dttm", "end_dttm"]]

sofa_post_scores = compute_sofa_polars(
    data_directory=TABLES_PATH,
    cohort_df=pl.from_pandas(sofa_cohort_post),
    filetype=FILE_TYPE,
    id_name="encounter_block",
    extremal_type="worst",
    fill_na_scores_with_zero=False,
    remove_outliers=True,
    timezone=TIMEZONE,
).to_pandas()

sofa_post_keep = ["encounter_block"] + [c for c in SOFA_COLS if c in sofa_post_scores.columns]
sofa_post_df = sofa_post_scores[sofa_post_keep].rename(
    columns={c: f"{c}_post" for c in SOFA_COLS if c in sofa_post_scores.columns}
)
del sofa_post_scores, sofa_cohort_post
gc.collect()
print(f"  SOFA post-CRRT: {len(sofa_post_df)} encounters")


# ===================================================================
# STEP 3: Compute sepsis flags (ASE ±72h of CRRT initiation)
# ===================================================================
print("Step 3: Computing sepsis flags …")
from clifpy.utils.ase import compute_ase  # noqa: E402

hosp_ids = outcomes_df["hospitalization_id"].unique().tolist()
# Use filtered config if pre-filtering was applied, otherwise original
_ase_config_path = Path(TABLES_PATH) / "config.json"
if not _ase_config_path.exists():
    _ase_config_path = project_root / "config" / "config.json"
sepsis_raw = compute_ase(
    hospitalization_ids=hosp_ids,
    config_path=str(_ase_config_path),
    apply_rit=True,
    include_lactate=True,
    verbose=True,
)

# Identify encounters with sepsis onset within ±72h of CRRT initiation
sepsis_pos = sepsis_raw[sepsis_raw["sepsis"] == 1].copy()
sepsis_pos = sepsis_pos.merge(eb_map, on="hospitalization_id", how="inner")
sepsis_pos = sepsis_pos.merge(
    crrt_initiation[["encounter_block", "crrt_initiation_time"]],
    on="encounter_block",
    how="inner",
)

ebs_with_sepsis: set = set()
onset_col = "ase_onset_w_lactate_dttm"
if onset_col in sepsis_pos.columns and len(sepsis_pos) > 0:
    sepsis_pos["onset_dttm"] = pd.to_datetime(sepsis_pos[onset_col], utc=True)
    crrt_utc = sepsis_pos["crrt_initiation_time"].dt.tz_convert("UTC")
    hours_diff = (sepsis_pos["onset_dttm"] - crrt_utc).dt.total_seconds() / 3600
    ebs_with_sepsis = set(sepsis_pos.loc[hours_diff.abs() <= 72, "encounter_block"])

sepsis_flag_df = pd.DataFrame({
    "encounter_block": crrt_initiation["encounter_block"],
    "sepsis_within_window": crrt_initiation["encounter_block"].isin(ebs_with_sepsis),
})
del sepsis_raw, sepsis_pos
gc.collect()
print(f"  Sepsis within ±72h: {sepsis_flag_df['sepsis_within_window'].sum()}/{len(sepsis_flag_df)}")


# ===================================================================
# STEP 4: Extract variables from wide_df
# ===================================================================
print("Step 4: Loading wide_df (selective columns) …")
# Determine which columns to load
lab_wide = [f"lab_{c}" for c in LAB_COLS]
vaso_wide = [f"med_cont_{c}" for c in VASO_COLS]
nee_wide = [f"med_cont_{c}" for c in NEE_COLS]
resp_wide = [f"resp_{c}" for c in RESP_SETTING_COLS] + ["resp_device_category"]
adt_wide = ["adt_location_category", "adt_location_type"]

available_cols = {f.name for f in pq.read_schema(INTERMEDIATE_DIR / "wide_df.parquet")}
needed_cols = (["hospitalization_id", "event_dttm", "crrt_ultrafiltration_out"]
               + lab_wide + vaso_wide + nee_wide + resp_wide + adt_wide)
needed_cols = [c for c in needed_cols if c in available_cols]

wide_df = load_intermediate(INTERMEDIATE_DIR / "wide_df.parquet", columns=needed_cols)
wide_df = wide_df.merge(eb_map, on="hospitalization_id", how="inner")
wide_df = wide_df.merge(
    crrt_initiation[["encounter_block", "crrt_initiation_time"]],
    on="encounter_block",
    how="inner",
)
wide_df["hours_from_crrt"] = (
    (wide_df["event_dttm"] - wide_df["crrt_initiation_time"]).dt.total_seconds() / 3600
)
print(f"  wide_df: {wide_df.shape}")

# Net UF intensity (mL/kg/hr) over the FIRST 72 h of CRRT — the Murugan / 2025
# literature definition (total UFnet / weight / duration), replacing the first-3h
# initiation snapshot. The fixed early window matches the published UFnet
# intensity exposure and avoids the duration/survivorship confound of a
# whole-course average. Shares ONE definition with 03's figures via
# pipeline_helpers (same 72 h window).
UF_INTENSITY_WINDOW_H = 72.0
_uf_intensity_by_eb = None
if HAS_CRRT_SETTINGS and "crrt_ultrafiltration_out" in wide_df.columns:
    _wt_eb = index_crrt_df[["encounter_block", "weight_kg"]].drop_duplicates("encounter_block")
    _uf_recs = wide_df[["encounter_block", "hours_from_crrt", "crrt_ultrafiltration_out"]].merge(
        _wt_eb, on="encounter_block", how="left")
    _uf_recs["uf_rate"] = (pd.to_numeric(_uf_recs["crrt_ultrafiltration_out"], errors="coerce")
                           / pd.to_numeric(_uf_recs["weight_kg"], errors="coerce"))
    _uf_int = course_average_intensity(_uf_recs, "encounter_block", "hours_from_crrt",
                                       "uf_rate", max_h=UF_INTENSITY_WINDOW_H)
    _uf_intensity_by_eb = _uf_int.set_index("encounter_block")["intensity"]
    _med = _uf_int["intensity"].median() if len(_uf_int) else float("nan")
    print(f"    net UF intensity (first 72h): {len(_uf_int)} encounters, median {_med:.2f} mL/kg/hr")

# -------------------------------------------------------------------
# 4a: 12h forward-fill for labs
# -------------------------------------------------------------------
print("  4a: Applying unlimited forward-fill for labs …")
lab_cols_present = [c for c in lab_wide if c in wide_df.columns]
vaso_present = [c for c in vaso_wide if c in wide_df.columns]
nee_present = [c for c in nee_wide if c in wide_df.columns]
resp_present = [c for c in resp_wide if c in wide_df.columns]

wide_df = wide_df.sort_values(["encounter_block", "event_dttm"])
for col in lab_cols_present:
    wide_df[col] = wide_df.groupby("encounter_block")[col].ffill()

# -------------------------------------------------------------------
# 4b: Compute per-encounter means for each time window
# -------------------------------------------------------------------
print("  4b: Computing last values per window for all continuous variables …")
continuous_wide = (
    lab_cols_present
    + [c for c in vaso_present]
    + [c for c in nee_present]
    + [c for c in resp_present if c != "resp_device_category"]
)


def _strip_prefix(col_name: str) -> str:
    """Remove wide_df prefix (lab_, resp_, med_cont_) from column name."""
    for prefix in ("lab_", "resp_", "med_cont_"):
        if col_name.startswith(prefix):
            return col_name[len(prefix):]
    return col_name


def window_last_values(df, min_h, max_h, suffix, cols):
    """Per-encounter last non-null value of continuous columns within [min_h, max_h].
    Pass None for min_h/max_h to leave that bound open (full stay)."""
    mask = pd.Series(True, index=df.index)
    if min_h is not None:
        mask &= df["hours_from_crrt"] >= min_h
    if max_h is not None:
        mask &= df["hours_from_crrt"] <= max_h
    last_vals = df[mask].groupby("encounter_block")[cols].last().reset_index()
    rename = {c: f"{_strip_prefix(c)}_{suffix}" for c in cols}
    return last_vals.rename(columns=rename)


baseline_last = window_last_values(wide_df, -12, 3, "baseline", continuous_wide)
post72h_last = window_last_values(wide_df, 0, 72, "post72h", continuous_wide)
post_crrt_last = window_last_values(wide_df, 0, None, "post_crrt", continuous_wide)

print(f"    baseline: {len(baseline_last)}, "
      f"post72h: {len(post72h_last)}, post_crrt: {len(post_crrt_last)} encounters")

# -------------------------------------------------------------------
# 4c: Extract categorical baselines (ADT location, resp device)
# -------------------------------------------------------------------
print("  4c: Categorical baselines (ADT, device) …")
pre_crrt = wide_df[wide_df["hours_from_crrt"] < 0].copy()
pre_crrt = pre_crrt.sort_values(["encounter_block", "event_dttm"])

# ADT
adt_present = [c for c in adt_wide if c in pre_crrt.columns]
if adt_present:
    pre_crrt[adt_present] = pre_crrt.groupby("encounter_block")[adt_present].ffill()
    adt_baseline = pre_crrt.groupby("encounter_block").tail(1)[["encounter_block"] + adt_present].copy()
    adt_baseline = adt_baseline.rename(columns={
        "adt_location_category": "location_category",
        "adt_location_type": "location_type",
    })
else:
    adt_baseline = pd.DataFrame({"encounter_block": crrt_initiation["encounter_block"]})

# Respiratory device (categorical only)
if "resp_device_category" in pre_crrt.columns:
    pre_crrt["resp_device_category"] = pre_crrt.groupby("encounter_block")["resp_device_category"].ffill()
    device_baseline = pre_crrt.groupby("encounter_block").tail(1)[
        ["encounter_block", "resp_device_category"]
    ].copy()
    device_baseline = device_baseline.rename(columns={"resp_device_category": "device_category"})
else:
    device_baseline = pd.DataFrame({"encounter_block": crrt_initiation["encounter_block"]})

# -------------------------------------------------------------------
# 4d: Compute IMV duration (first IMV to 3h-gap stop) from wide_df
# -------------------------------------------------------------------
print("  4d: Computing IMV duration (3h gap) …")
if "resp_device_category" in wide_df.columns:
    imv_rows = wide_df[wide_df["resp_device_category"].str.lower() == "imv"].copy()
    imv_rows = imv_rows.sort_values(["encounter_block", "event_dttm"])

    def _imv_duration(group):
        times = group["event_dttm"].dropna().sort_values()
        if len(times) == 0:
            return pd.Series({"imv_duration_hours": 0.0})
        start = times.iloc[0]
        gaps = times.diff() > pd.Timedelta(hours=3)
        if gaps.any():
            gap_pos = times.index.get_loc(gaps.idxmax())
            end = times.iloc[gap_pos - 1]
        else:
            end = times.iloc[-1]
        return pd.Series({"imv_duration_hours": (end - start).total_seconds() / 3600})

    imv_duration_df = imv_rows.groupby("encounter_block").apply(_imv_duration).reset_index()
    del imv_rows
else:
    imv_duration_df = pd.DataFrame({
        "encounter_block": crrt_initiation["encounter_block"],
        "imv_duration_hours": pd.NA,
    })
print(f"    IMV duration: {len(imv_duration_df)} encounters")

# -------------------------------------------------------------------
# 4e: Compute vasopressor duration (first vaso to 3h-gap stop)
# -------------------------------------------------------------------
print("  4e: Computing vasopressor duration (3h gap) …")
vaso_cols_present = [c for c in vaso_present if c in wide_df.columns]
if vaso_cols_present:
    vaso_rows = wide_df[wide_df[vaso_cols_present].notna().any(axis=1)].copy()
    vaso_rows = vaso_rows.sort_values(["encounter_block", "event_dttm"])

    def _vaso_duration(group):
        times = group["event_dttm"].dropna().sort_values()
        if len(times) == 0:
            return pd.Series({"vaso_duration_hours": 0.0})
        start = times.iloc[0]
        gaps = times.diff() > pd.Timedelta(hours=3)
        if gaps.any():
            gap_pos = times.index.get_loc(gaps.idxmax())
            end = times.iloc[gap_pos - 1]
        else:
            end = times.iloc[-1]
        return pd.Series({
            "vaso_duration_hours": (end - start).total_seconds() / 3600,
        })

    vaso_duration_df = vaso_rows.groupby("encounter_block").apply(_vaso_duration).reset_index()
    del vaso_rows
else:
    vaso_duration_df = pd.DataFrame({
        "encounter_block": crrt_initiation["encounter_block"],
        "vaso_duration_hours": pd.NA,
    })
print(f"    Vasopressor duration: {len(vaso_duration_df)} encounters")

# Free memory
del wide_df, pre_crrt
gc.collect()


# ===================================================================
# STEP 5: Merge into analysis-ready DataFrame
# ===================================================================
print("Step 5: Merging all data …")
analysis_df = outcomes_df.copy()

# CRRT settings, duration & dose
crrt_merge_cols = ["encounter_block"]
if HAS_CRRT_SETTINGS:
    crrt_merge_cols.append("crrt_mode_category")
    crrt_merge_cols += [c for c in CRRT_SETTING_COLS + ["duration_hours", "duration_days", "crrt_dose_ml_kg_hr"] if c in index_crrt_df.columns]
else:
    crrt_merge_cols += [c for c in ["duration_hours", "duration_days"] if c in index_crrt_df.columns]
analysis_df = analysis_df.merge(index_crrt_df[crrt_merge_cols], on="encounter_block", how="left")

# NOTE: crrt_initiation_time already merged into outcomes_df in Step 1 (90-day cap)

# SOFA at initiation + SOFA 72h + SOFA post-CRRT
analysis_df = analysis_df.merge(sofa_df, on="encounter_block", how="left")
analysis_df = analysis_df.merge(sofa_72h_df, on="encounter_block", how="left")
analysis_df = analysis_df.merge(sofa_post_df, on="encounter_block", how="left")

# Sepsis
analysis_df = analysis_df.merge(sepsis_flag_df, on="encounter_block", how="left")

# Window last values (labs, vasopressors, respiratory settings)
analysis_df = analysis_df.merge(baseline_last, on="encounter_block", how="left")
analysis_df = analysis_df.merge(post72h_last, on="encounter_block", how="left")
analysis_df = analysis_df.merge(post_crrt_last, on="encounter_block", how="left")

# Categorical baselines (ADT, device)
analysis_df = analysis_df.merge(adt_baseline, on="encounter_block", how="left")
analysis_df = analysis_df.merge(device_baseline, on="encounter_block", how="left")

# IMV & vasopressor duration
analysis_df = analysis_df.merge(imv_duration_df, on="encounter_block", how="left")
analysis_df = analysis_df.merge(vaso_duration_df, on="encounter_block", how="left")

# Derived columns
analysis_df["alive_72h"] = (
    analysis_df["death_dttm"].isna()
    | ((analysis_df["death_dttm"] - analysis_df["crrt_initiation_time"]).dt.total_seconds() / 3600 > 72)
)
analysis_df["dead_within_72h"] = ~analysis_df["alive_72h"]
if HAS_CRRT_SETTINGS:
    analysis_df["last_crrt_mode"] = analysis_df["crrt_mode_category"]

# Normalize all _category / _type columns to lowercase
for col in analysis_df.columns:
    if (col.endswith("_category") or col == "location_type") and analysis_df[col].dtype == "object":
        analysis_df[col] = analysis_df[col].str.lower()

# Save
analysis_df.to_parquet(INTERMEDIATE_DIR / "tableone_analysis_df.parquet", index=False)
print(f"  Analysis DataFrame: {analysis_df.shape}")


# ===================================================================
# STEP 6: Table 1 - baseline characteristics by CRRT dose band
# ===================================================================
# Single project Table 1: full analytic CRRT cohort, stratified by initial
# (median-first-3h) CRRT dose band - <20 / 20-30 / >30 mL/kg/hr (the 20-30
# band brackets the KDIGO-recommended delivered target). gtsummary shape
# consumed by the combine renderer (09) and manuscript builder (10). p-value
# is across the three bands (Kruskal-Wallis for continuous, chi-square for
# categorical).
print("Step 6: Generating Table 1 (baseline by CRRT dose band) ...")
from scipy import stats as _stats

_wt = index_crrt_df[["encounter_block", "weight_kg"]].drop_duplicates("encounter_block")
t1 = analysis_df.merge(_wt, on="encounter_block", how="left")
if _uf_intensity_by_eb is not None:
    # Course-average net UF intensity (mL/kg/hr) — NOT the first-3h snapshot.
    t1["net_uf_intensity"] = t1["encounter_block"].map(_uf_intensity_by_eb)

# Accurate baseline labs for Table 1: the value MEASURED nearest CRRT
# initiation within [-12, +3] h. This avoids the unlimited-forward-fill
# staleness in the analysis_df *_baseline columns (a point lab does not
# persist the way an infusion rate does), which otherwise carried stale
# pre-window values into the baseline (e.g. lactate median 3.3 vs ~4.6
# measured nearest t=0). NOTE: the *_baseline columns consumed by 04/causal
# still use the forward-filled values - addressing that is a separate,
# causal-affecting decision.
_LABS = ["creatinine", "bun", "lactate", "bicarbonate", "potassium", "phosphate"]
_lab_cols = [f"lab_{x}" for x in _LABS]
_avail = {f.name for f in pq.read_schema(INTERMEDIATE_DIR / "wide_df.parquet")}
_need = ["hospitalization_id", "event_dttm"] + [c for c in _lab_cols if c in _avail]
_wd = load_intermediate(INTERMEDIATE_DIR / "wide_df.parquet", columns=_need)
_wd = _wd.merge(eb_map, on="hospitalization_id", how="inner").merge(
    crrt_initiation[["encounter_block", "crrt_initiation_time"]], on="encounter_block", how="inner")
_wd["_h"] = (_wd["event_dttm"] - _wd["crrt_initiation_time"]).dt.total_seconds() / 3600.0
_wd = _wd[(_wd["_h"] >= -12) & (_wd["_h"] <= 3)].copy()
_wd["_abs"] = _wd["_h"].abs()
for x in _LABS:
    c = f"lab_{x}"
    if c in _wd.columns:
        near = (_wd.dropna(subset=[c]).sort_values("_abs")
                .groupby("encounter_block")[c].first())
        t1[f"{x}_t1"] = t1["encounter_block"].map(near)
del _wd

BANDS = ["<20 mL/kg/hr", "20-30 mL/kg/hr", ">30 mL/kg/hr"]
_dose = pd.to_numeric(t1.get("crrt_dose_ml_kg_hr"), errors="coerce")
# Exact band membership: <20, 20-30 (inclusive of both ends), >30. Mask-based
# (not pd.cut) so dose==30 lands in the 20-30 guideline band, not >30; missing
# dose stays NA.
t1["dose_band"] = pd.Series(pd.NA, index=t1.index, dtype="object")
t1.loc[_dose < 20, "dose_band"] = BANDS[0]
t1.loc[(_dose >= 20) & (_dose <= 30), "dose_band"] = BANDS[1]
t1.loc[_dose > 30, "dose_band"] = BANDS[2]
t1["dose_band"] = pd.Categorical(t1["dose_band"], categories=BANDS, ordered=True)

# Derived display columns (computed before building strata frames)
t1["_female"] = t1["sex_category"].astype("string").str.lower().eq("female")
t1["_imv"] = t1["device_category"].astype("string").str.lower().eq("imv") if "device_category" in t1.columns else False
t1["_death30"] = (pd.to_numeric(t1["death_30d"], errors="coerce") == 1) if "death_30d" in t1.columns else False
_r = t1["race_category"].astype("string").str.lower().fillna("unknown")
t1["_race_grp"] = _r.map(lambda s: "Black" if "black" in s else ("White" if s == "white" else "Other"))

GROUPS = ["Overall"] + BANDS
gframe = {"Overall": t1}
for _b in BANDS:
    gframe[_b] = t1[t1["dose_band"] == _b]
gn = {g: len(gframe[g]) for g in GROUPS}
GHDR = {g: f"**{g}**  N = {gn[g]:,}" for g in GROUPS}
_dosed = t1[t1["dose_band"].isin(BANDS)]
COL = "**Characteristic**"
PCOL = "**p-value**"
columns = [COL] + [GHDR[g] for g in GROUPS] + [PCOL]


def _iqr(s, dec):
    s = pd.to_numeric(s, errors="coerce").dropna()
    if s.empty:
        return "NA"
    return f"{s.median():.{dec}f} ({s.quantile(0.25):.{dec}f}, {s.quantile(0.75):.{dec}f})"


def _npct(count, total):
    return f"{int(count):,} ({100 * count / total:.0f}%)" if total else "NA"


def _fmt_p(p):
    if p is None or pd.isna(p):
        return "NA"
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}" if p < 0.01 else f"{p:.2f}"


def _p_cont(col):
    arrs = []
    for _b in BANDS:
        a = pd.to_numeric(gframe[_b][col], errors="coerce").dropna()
        if len(a) >= 2:
            arrs.append(a)
    if len(arrs) < 2:
        return None
    try:
        return _stats.kruskal(*arrs).pvalue
    except Exception:
        return None


def _p_cat(series):
    ct = pd.crosstab(series, _dosed["dose_band"])
    if ct.shape[0] < 2 or ct.shape[1] < 2:
        return None
    try:
        return _stats.chi2_contingency(ct)[1]
    except Exception:
        return None


_rows = []


def _row_cont(label, col, dec):
    if col not in t1.columns:
        return
    r = {COL: f"__{label}__", PCOL: _fmt_p(_p_cont(col))}
    for g in GROUPS:
        r[GHDR[g]] = _iqr(gframe[g][col], dec)
    _rows.append(r)


def _row_binary(label, boolcol):
    r = {COL: f"__{label}__", PCOL: _fmt_p(_p_cat(_dosed[boolcol]))}
    for g in GROUPS:
        f = gframe[g]
        r[GHDR[g]] = _npct(int(f[boolcol].sum()), len(f))
    _rows.append(r)


def _row_multi(label, col, level_order):
    _rows.append({COL: f"__{label}__", PCOL: _fmt_p(_p_cat(_dosed[col])),
                  **{GHDR[g]: "" for g in GROUPS}})
    for lv, disp in level_order:
        r = {COL: disp, PCOL: ""}
        for g in GROUPS:
            f = gframe[g]
            r[GHDR[g]] = _npct(int((f[col] == lv).sum()), len(f))
        _rows.append(r)


# Demographics
_row_cont("Age at Admission (years)", "age_at_admission", 0)
_row_binary("Female (%)", "_female")
_row_multi("Race", "_race_grp", [("Black", "Black"), ("White", "White"), ("Other", "Other")])
# Severity and labs at CRRT initiation
_row_cont("SOFA Score", "sofa_total", 0)
_row_cont("Creatinine (mg/dL)", "creatinine_t1", 2)
_row_cont("BUN (mg/dL)", "bun_t1", 0)
_row_cont("Lactate (mmol/L)", "lactate_t1", 1)
_row_cont("Bicarbonate (mEq/L)", "bicarbonate_t1", 1)
_row_cont("Potassium (mEq/L)", "potassium_t1", 2)
_row_cont("Phosphate (mg/dL)", "phosphate_t1", 1)
_row_cont("NE Equivalent (mcg/kg/min)", "nee_baseline", 2)
_row_binary("On IMV (%)", "_imv")
# CRRT practice descriptors
_row_cont("CRRT Dose (mL/kg/hr)", "crrt_dose_ml_kg_hr", 1)
if HAS_CRRT_SETTINGS:
    _modes = list(_dosed["crrt_mode_category"].dropna().value_counts().index)
    if _modes:
        _row_multi("CRRT Modality", "crrt_mode_category", [(m, str(m).upper()) for m in _modes])
    _row_cont("Net UF Intensity (mL/kg/hr, first 72h)", "net_uf_intensity", 2)
# Outcome
_row_binary("30-Day Mortality (%)", "_death30")


# ===================================================================
# STEP 7: Write Table 1 (gtsummary-shaped CSV + standalone HTML)
# ===================================================================
import re as _re
print("Step 7: Writing Table 1 ...")
table1_df = pd.DataFrame(_rows, columns=columns)
output_csv = FINAL_DIR / f"{SITE_NAME}_table1_crrt.csv"
table1_df.to_csv(output_csv, index=False)
print(f"  CSV: {output_csv}  ({len(table1_df)} rows; N {gn})")

_disp = table1_df.copy()
_disp.columns = [_re.sub(r"\*+", "", c).strip() for c in _disp.columns]
_disp[_disp.columns[0]] = _disp[_disp.columns[0]].apply(
    lambda v: _re.sub(r"^__(.+?)__$", r"<b>\1</b>",
                      ("&nbsp;&nbsp;&nbsp;" + v if not str(v).startswith("__") else v)))
html_path = FINAL_DIR / f"{SITE_NAME}_table1_crrt.html"
with open(html_path, "w", encoding="utf-8") as _f:
    _f.write(
        "<!DOCTYPE html><html><head><meta charset='UTF-8'>"
        "<title>Table 1 - CRRT cohort by dose band</title>"
        "<style>body{font-family:Arial,sans-serif;margin:20px}"
        "table{border-collapse:collapse}th,td{border:1px solid #ddd;padding:6px 10px;text-align:left}"
        "th{background:#1e417c;color:#fff}</style></head><body>"
        f"<h2>Table 1: Baseline Characteristics by CRRT Dose Band - {SITE_NAME}</h2>"
        "<p><em>Continuous: median (Q1, Q3); categorical: n (%). Labs and severity "
        "measured at or nearest CRRT initiation (-12 to +3 h); age at hospital admission. "
        "p-value across the three dose bands (Kruskal-Wallis / chi-square).</em></p>"
        + _disp.to_html(index=False, escape=False) + "</body></html>")
print(f"  HTML: {html_path}")
print("Done!")
