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

from sofa_calculator import compute_sofa_polars  # noqa: E402

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
INTERMEDIATE_DIR = project_root / "output" / "intermediate"
FINAL_DIR = project_root / "output" / "final"
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
outcomes_df = pd.read_parquet(INTERMEDIATE_DIR / "outcomes_df.parquet")
index_crrt_df = pd.read_parquet(INTERMEDIATE_DIR / "index_crrt_df.parquet")
crrt_initiation = pd.read_parquet(INTERMEDIATE_DIR / "crrt_initiation.parquet")

# 1:1 hospitalization_id <-> encounter_block
eb_map = index_crrt_df[["hospitalization_id", "encounter_block"]].copy()
print(f"  {len(outcomes_df)} encounters loaded")

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
print(f"  SOFA post-CRRT: {len(sofa_post_df)} encounters")


# ===================================================================
# STEP 3: Compute sepsis flags (ASE ±72h of CRRT initiation)
# ===================================================================
print("Step 3: Computing sepsis flags …")
from clifpy.utils.ase import compute_ase  # noqa: E402

hosp_ids = outcomes_df["hospitalization_id"].unique().tolist()
sepsis_raw = compute_ase(
    hospitalization_ids=hosp_ids,
    config_path=str(project_root / "config" / "config.json"),
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
needed_cols = ["hospitalization_id", "event_dttm"] + lab_wide + vaso_wide + nee_wide + resp_wide + adt_wide
needed_cols = [c for c in needed_cols if c in available_cols]

wide_df = pd.read_parquet(INTERMEDIATE_DIR / "wide_df.parquet", columns=needed_cols)
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
# STEP 6: Generate Table 1
# ===================================================================
print("Step 6: Generating Table 1 …")

# --- Encounter-level subgroups (for clinical measurements) ---
enc_subgroups = {
    "Baseline": analysis_df,
    "At 72h - Survivors": analysis_df[analysis_df["alive_72h"]],
    "At 72h - Non-survivors": analysis_df[analysis_df["dead_within_72h"]],
    "At discharge - Survivors": analysis_df[analysis_df["in_hosp_death"] == 0],
    "At discharge - Non-survivors": analysis_df[analysis_df["in_hosp_death"] == 1],
}

# --- Patient-level subgroups (for static demographics & death) ---
# Deduplicate: one row per patient, keeping first CRRT encounter
patient_df = analysis_df.sort_values("crrt_initiation_time").drop_duplicates(subset="patient_id", keep="first")
pat_subgroups = {
    "Baseline": patient_df,
    "At 72h - Survivors": patient_df[patient_df["alive_72h"]],
    "At 72h - Non-survivors": patient_df[patient_df["dead_within_72h"]],
    "At discharge - Survivors": patient_df[patient_df["in_hosp_death"] == 0],
    "At discharge - Non-survivors": patient_df[patient_df["in_hosp_death"] == 1],
}

SG_NAMES = list(enc_subgroups.keys())

# Build rows as 7-element lists matching template structure
rows: list[list[str]] = []

# Machine-readable long-format rows for multi-site aggregation
long_rows: list[dict] = []

# ── Header rows ──
rows.append(["CRRT Population (excluding ESRD patients)"] + [""] * (len(SG_NAMES) + 1))
rows.append(["", "", SITE_NAME] + [""] * (len(SG_NAMES) - 1))
rows.append(["", ""] + SG_NAMES)


# ── Helper: add categorical variable rows ──
def add_categorical(var_name: str, levels: list, col_name: str | None = None,
                    patient_level: bool = False):
    """Add rows for a categorical variable. All string levels are lowercase."""
    col = col_name or var_name
    sgs = pat_subgroups if patient_level else enc_subgroups
    for level in levels:
        if level is None:
            level_str = "None"
        elif isinstance(level, bool):
            level_str = str(level).upper()
        else:
            level_str = str(level)
        row = [f"{var_name}, n (%)", level_str]
        for sg_name, sg_df in sgs.items():
            if col not in sg_df.columns:
                row.append("—")
                continue
            count = int(sg_df[col].isna().sum()) if level is None else int((sg_df[col] == level).sum())
            total = len(sg_df)
            row.append(fmt_n_pct(count, total))
            long_rows.append({
                "variable": var_name, "level": level_str, "subgroup": sg_name,
                "stat_type": "categorical", "n": count, "total": total,
                "median": None, "q25": None, "q75": None,
            })
        rows.append(row)


# ── Helper: add continuous variable row ──
def add_continuous(label: str, col_name: str, col_map: dict | None = None,
                   patient_level: bool = False):
    """Add a row for a continuous variable.

    col_map: optional dict {subgroup_name: column_name} to use a different
    column for specific subgroups (e.g., window-specific means).
    Falls back to col_name for any subgroup not in col_map.
    """
    sgs = pat_subgroups if patient_level else enc_subgroups
    row = [label, ""]
    for sg_name, sg_df in sgs.items():
        col = col_map.get(sg_name, col_name) if col_map else col_name
        if col in sg_df.columns:
            valid = sg_df[col].dropna()
            row.append(fmt_median_iqr(sg_df[col]))
            long_rows.append({
                "variable": label, "level": "", "subgroup": sg_name,
                "stat_type": "continuous", "n": len(valid), "total": len(sg_df),
                "median": float(valid.median()) if len(valid) > 0 else None,
                "q25": float(valid.quantile(0.25)) if len(valid) > 0 else None,
                "q75": float(valid.quantile(0.75)) if len(valid) > 0 else None,
            })
        else:
            row.append("—")
    rows.append(row)


# ── Per-subgroup column mapping for windowed continuous variables ──
WINDOW_SUFFIXES = {
    "Baseline": "baseline",
    "At 72h - Survivors": "post72h",
    "At 72h - Non-survivors": "post72h",
    "At discharge - Survivors": "post_crrt",
    "At discharge - Non-survivors": "post_crrt",
}


def _windowed_col_map(base: str) -> dict:
    """Build {subgroup: column} mapping for a windowed variable."""
    return {sg: f"{base}_{suffix}" for sg, suffix in WINDOW_SUFFIXES.items()}


# ── Counts ──
n_pat_row = ["n patients", ""]
for sg_name, sg_df in pat_subgroups.items():
    n_pat_row.append(str(len(sg_df)))
    long_rows.append({
        "variable": "n_patients", "level": "", "subgroup": sg_name,
        "stat_type": "count", "n": len(sg_df), "total": len(sg_df),
        "median": None, "q25": None, "q75": None,
    })
rows.append(n_pat_row)

n_hosp_row = ["n hospitalizations ", ""]
for sg_name, sg_df in enc_subgroups.items():
    n_hosp_row.append(str(len(sg_df)))
    long_rows.append({
        "variable": "n_hospitalizations", "level": "", "subgroup": sg_name,
        "stat_type": "count", "n": len(sg_df), "total": len(sg_df),
        "median": None, "q25": None, "q75": None,
    })
rows.append(n_hosp_row)

# ── Categorical variables (all levels lowercase to match data) ──
# Patient-level: demographics & death
add_categorical("sex_category", ["female", "male", None], patient_level=True)
add_categorical("race_category", [
    "american indian or alaska native", "asian", "black or african american",
    "native hawaiian or other pacific islander", None, "other", "unknown", "white",
], patient_level=True)
add_categorical("ethnicity_category", ["hispanic", "non-hispanic", None, "unknown"], patient_level=True)

# Encounter-level: clinical context at CRRT initiation
add_categorical("location_category", [None, "icu", "procedural", "ward", "ed"])
add_categorical("location_type", [
    None, "general_icu", "medical_icu", "mixed_cardiothoracic_icu",
    "mixed_neuro_icu", "surgical_icu",
])
add_categorical("device_category", [
    None, "cpap", "face mask", "high flow nc", "imv", "nasal cannula",
    "nippv", "other", "room air", "trach collar",
])
# Mortality within window — contextual per subgroup:
#   Baseline / Discharge columns use in-hospital mortality
#   72h columns use 72h mortality
_wm_row_false = ["Mortality within window, n (%)", "FALSE"]
_wm_row_true = ["Mortality within window, n (%)", "TRUE"]
_wm_col_map = {
    "Baseline": "in_hosp_death",
    "At 72h - Survivors": "dead_within_72h",
    "At 72h - Non-survivors": "dead_within_72h",
    "At discharge - Survivors": "in_hosp_death",
    "At discharge - Non-survivors": "in_hosp_death",
}
for sg_name, sg_df in pat_subgroups.items():
    col = _wm_col_map[sg_name]
    total = len(sg_df)
    true_count = int(sg_df[col].sum())
    false_count = total - true_count
    _wm_row_false.append(fmt_n_pct(false_count, total))
    _wm_row_true.append(fmt_n_pct(true_count, total))
    long_rows.append({
        "variable": "Mortality within window", "level": "FALSE", "subgroup": sg_name,
        "stat_type": "categorical", "n": false_count, "total": total,
        "median": None, "q25": None, "q75": None,
    })
    long_rows.append({
        "variable": "Mortality within window", "level": "TRUE", "subgroup": sg_name,
        "stat_type": "categorical", "n": true_count, "total": total,
        "median": None, "q25": None, "q75": None,
    })
rows.append(_wm_row_false)
rows.append(_wm_row_true)

# Mortality overall — in-hospital mortality across all subgroups
add_categorical("Mortality overall", [0, 1], col_name="in_hosp_death", patient_level=True)
add_categorical("sepsis_within_window", [False, True])

# ── LOS ──
add_continuous("hospital LOS (days), median [Q1,Q3]", "hosp_los_days")
add_continuous("ICU LOS (days), median [Q1,Q3]", "icu_los_days")

# ── Age (patient-level) ──
add_continuous("age_at_admission, median [Q1,Q3]", "age_at_admission", patient_level=True)

# ── Vasopressors (windowed per subgroup) ──
for v in VASO_COLS:
    add_continuous(f"{v}, median [Q1,Q3]", f"{v}_baseline", col_map=_windowed_col_map(v))

# ── NEE (windowed per subgroup) ──
for n in NEE_COLS:
    add_continuous(f"NEE (mcg/kg/min), median [Q1,Q3]", f"{n}_baseline", col_map=_windowed_col_map(n))

# ── Labs (windowed per subgroup) ──
for lab in LAB_COLS:
    add_continuous(f"{lab}, median [Q1,Q3]", f"{lab}_baseline", col_map=_windowed_col_map(lab))

# ── Respiratory settings (windowed per subgroup) ──
for r in RESP_SETTING_COLS:
    add_continuous(f"{r}, median [Q1,Q3]", f"{r}_baseline", col_map=_windowed_col_map(r))

# ── SOFA scores (initiation for Baseline, 72h for 72h cols, post for discharge cols) ──
for s in SOFA_COLS:
    add_continuous(f"{s}, median [Q1,Q3]", s, col_map={
        "At 72h - Survivors": f"{s}_72h",
        "At 72h - Non-survivors": f"{s}_72h",
        "At discharge - Survivors": f"{s}_post",
        "At discharge - Non-survivors": f"{s}_post",
    })

# ── Duration of IMV and vasopressor support ──
add_continuous("Duration of IMV (hours), median [Q1,Q3]", "imv_duration_hours")
add_continuous("Duration of vasopressor support (hours), median [Q1,Q3]", "vaso_duration_hours")

# ── CRRT details ──
rows.append(["CRRT details", ""] + [""] * len(SG_NAMES))

# Duration of CRRT
add_continuous("Duration of CRRT (hours), median [Q1,Q3]", "duration_hours")

if HAS_CRRT_SETTINGS:
    # CRRT mode distribution
    modes = sorted(analysis_df["last_crrt_mode"].dropna().unique())
    add_categorical("last_crrt_mode", modes)


    # ── Helper: add continuous row filtered by CRRT mode ──
    def add_continuous_by_mode(label: str, col_name: str, mode_value: str):
        """Add a row for a continuous variable, filtered to a specific CRRT mode."""
        row = [label, ""]
        for sg_name, sg_df in enc_subgroups.items():
            mode_df = sg_df[sg_df["last_crrt_mode"] == mode_value]
            if col_name in mode_df.columns and len(mode_df) > 0:
                valid = mode_df[col_name].dropna()
                row.append(fmt_median_iqr(mode_df[col_name]))
                long_rows.append({
                    "variable": label, "level": mode_value, "subgroup": sg_name,
                    "stat_type": "continuous", "n": len(valid), "total": len(mode_df),
                    "median": float(valid.median()) if len(valid) > 0 else None,
                    "q25": float(valid.quantile(0.25)) if len(valid) > 0 else None,
                    "q75": float(valid.quantile(0.75)) if len(valid) > 0 else None,
                })
            else:
                row.append("—")
        rows.append(row)


    # Per-mode CRRT settings and dose
    crrt_detail_cols = [
        ("crrt_dose_ml_kg_hr", "crrt_dose_ml_kg_hr, median [Q1,Q3]"),
        ("blood_flow_rate", "blood_flow_rate, median [Q1,Q3]"),
        ("pre_filter_replacement_fluid_rate", "pre_filter_replacement_fluid_rate, median [Q1,Q3]"),
        ("post_filter_replacement_fluid_rate", "post_filter_replacement_fluid_rate, median [Q1,Q3]"),
        ("dialysate_flow_rate", "dialysate_flow_rate, median [Q1,Q3]"),
        ("ultrafiltration_out", "ultrafiltration_out, median [Q1,Q3]"),
    ]

    for mode in modes:
        rows.append([mode, ""] + [""] * len(SG_NAMES))
        for col, label in crrt_detail_cols:
            add_continuous_by_mode(label, col, mode)
else:
    rows.append(["CRRT settings not available at this site", ""] + [""] * len(SG_NAMES))


# ===================================================================
# STEP 7: Write output
# ===================================================================
print("Step 7: Writing output …")

table1_df = pd.DataFrame(rows)

# CSV (formatted, for display)
output_csv = FINAL_DIR / "table1_crrt.csv"
table1_df.to_csv(output_csv, index=False, header=False)
print(f"  CSV: {output_csv}")

# Long-format CSV (machine-readable, for multi-site aggregation)
long_df = pd.DataFrame(long_rows)
long_df["site"] = SITE_NAME
output_long_csv = FINAL_DIR / "table1_crrt_long.csv"
long_df.to_csv(output_long_csv, index=False)
print(f"  Long CSV: {output_long_csv}")

# HTML
html_path = FINAL_DIR / "table1_crrt.html"
# Use rows 3+ as the actual data table (skip the 3 header meta-rows)
data_df = pd.DataFrame(rows[3:], columns=["Variable", "Level"] + SG_NAMES)

html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Table 1 - CRRT Population</title>
    <style>
        body {{ font-family: 'Arial', sans-serif; margin: 20px; background: #f5f5f5; }}
        .container {{ max-width: 1400px; margin: 0 auto; background: white; padding: 30px;
                      box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        h1 {{ color: #333; border-bottom: 3px solid #4CAF50; padding-bottom: 10px; }}
        table {{ width: 100%; border-collapse: collapse; margin-top: 20px; font-size: 13px; }}
        th {{ background: #4CAF50; color: white; padding: 12px; text-align: left;
              border: 1px solid #ddd; }}
        td {{ padding: 10px; border: 1px solid #ddd; vertical-align: top; }}
        tr:nth-child(even) {{ background: #f9f9f9; }}
        tr:hover {{ background: #f0f0f0; }}
        .footer {{ margin-top: 20px; padding-top: 10px; border-top: 1px solid #ddd;
                   font-size: 11px; color: #666; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Table 1: CRRT Population (excluding ESRD patients) — {SITE_NAME}</h1>
        <p><em>Continuous variables: median [Q1, Q3] (n/N with data)</em></p>
        <p><em>Categorical variables: n (%)</em></p>
        {data_df.to_html(index=False, escape=False, classes='table')}
        <div class="footer">
            <p>Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
    </div>
</body>
</html>"""

with open(html_path, "w", encoding="utf-8") as f:
    f.write(html_content)
print(f"  HTML: {html_path}")

print("\nDone!")
