"""
Shared helpers for the CRRT epidemiology pipeline.

Provides:
  - validate_config()      : validate config.json fields
  - load_intermediate()    : read parquet with clear missing-file errors
  - safe_load_clif_table() : wrap clifpy table loading with clear errors
  - prefilter_clif_tables(): filter CLIF tables to cohort hospitalization_ids
  - get_tables_path()      : auto-detect whether to use original or filtered tables
"""

import json
import os
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

# ---------------------------------------------------------------------------
# Map intermediate files to the script that produces them
# ---------------------------------------------------------------------------
_PRODUCERS = {
    "outcomes_df.parquet": "00_cohort.py",
    "cohort_df.parquet": "00_cohort.py",
    "index_crrt_df.parquet": "00_cohort.py",
    "crrt_initiation.parquet": "00_cohort.py",
    "crrt_at_initiation.parquet": "00_cohort.py",
    "crrt_cohort.parquet": "00_cohort.py",
    "weight_df.parquet": "00_cohort.py",
    "dropped_missing_labs_blocks.parquet": "00_cohort.py",
    "wide_df.parquet": "01_create_wide_df.py",
    "tableone_analysis_df.parquet": "02_construct_crrt_tableone.py",
    "hospital_diagnosis_df.parquet": "00_cohort.py",
    "msm_competing_risk_df.parquet": "04_build_msm_competing_risk_df.py",
}

_PIPELINE_ORDER = (
    "00_cohort.py -> 01_create_wide_df.py -> 02_construct_crrt_tableone.py "
    "-> 03_crrt_visualizations.py -> 04_build_msm_competing_risk_df.py"
)

_VALID_FILE_TYPES = {"parquet", "csv", "fst"}

_REQUIRED_CONFIG_FIELDS = {
    "site_name": str,
    "tables_path": str,
    "file_type": str,
    "timezone": str,
    "project_root": str,
    "has_crrt_settings": bool,
}


# ---------------------------------------------------------------------------
# Config validation
# ---------------------------------------------------------------------------
def validate_config(config: dict) -> dict:
    """Validate config.json and set defaults for optional fields.

    Returns the config dict (with defaults filled in for optional fields).
    Raises KeyError, TypeError, ValueError, or FileNotFoundError on problems.
    """
    # --- Required fields ---
    for field, expected_type in _REQUIRED_CONFIG_FIELDS.items():
        if field not in config:
            raise KeyError(
                f"Missing required config field: '{field}'.\n"
                f"  See config/config_template.json for the expected format."
            )
        if not isinstance(config[field], expected_type):
            raise TypeError(
                f"Config field '{field}' must be {expected_type.__name__}, "
                f"got {type(config[field]).__name__}.\n"
                f"  See config/config_template.json for the expected format."
            )

    # --- tables_path must exist ---
    if not os.path.isdir(config["tables_path"]):
        raise FileNotFoundError(
            f"tables_path directory does not exist: {config['tables_path']}\n"
            f"  Check the 'tables_path' value in config/config.json."
        )

    # --- file_type must be valid ---
    if config["file_type"] not in _VALID_FILE_TYPES:
        raise ValueError(
            f"Config field 'file_type' must be one of {sorted(_VALID_FILE_TYPES)}, "
            f"got '{config['file_type']}'."
        )

    # --- Optional fields with defaults ---
    config.setdefault("admission_year_start", 2018)
    config.setdefault("admission_year_end", None)

    for field in ("admission_year_start", "admission_year_end"):
        val = config[field]
        if val is not None and not isinstance(val, int):
            raise TypeError(
                f"Config field '{field}' must be an integer or null, "
                f"got {type(val).__name__}."
            )

    return config


# ---------------------------------------------------------------------------
# Intermediate file loading
# ---------------------------------------------------------------------------
def load_intermediate(filepath, **kwargs) -> pd.DataFrame:
    """Read a parquet file with a clear error if it does not exist."""
    filepath = Path(filepath)
    if not filepath.exists():
        producer = _PRODUCERS.get(filepath.name, "an earlier pipeline step")
        raise FileNotFoundError(
            f"Required intermediate file not found:\n"
            f"  {filepath}\n"
            f"This file is produced by: {producer}\n"
            f"Please run that step first. Pipeline order:\n"
            f"  {_PIPELINE_ORDER}"
        )
    return pd.read_parquet(filepath, **kwargs)


# ---------------------------------------------------------------------------
# CLIF table loading
# ---------------------------------------------------------------------------
def safe_load_clif_table(clif, table_name, tables_path=None, **kwargs):
    """Wrap clif.load_table() with clear error messaging.

    Returns the loaded table object (e.g., clif.vitals).
    """
    try:
        clif.load_table(table_name, **kwargs)
        table = getattr(clif, table_name)
        print(f"  Loaded {table_name}: {len(table.df):,} rows")
        return table
    except FileNotFoundError:
        expected = (
            f"{tables_path}/clif_{table_name}.*" if tables_path
            else f"clif_{table_name}.*"
        )
        raise FileNotFoundError(
            f"CLIF table '{table_name}' not found.\n"
            f"  Expected location: {expected}\n"
            f"  Check that the file exists and that 'tables_path' in "
            f"config/config.json is correct.\n"
            f"  See config/clif_data_requirements.yaml for required tables."
        ) from None
    except Exception as e:
        raise type(e)(
            f"Failed to load CLIF table '{table_name}': {e}\n"
            f"  Check config/clif_data_requirements.yaml for required "
            f"columns and categories."
        ) from None


# ---------------------------------------------------------------------------
# CLIF table pre-filtering for SOFA/ASE
# ---------------------------------------------------------------------------
_CLIF_TABLES_FOR_SOFA_ASE = [
    "clif_hospitalization",
    "clif_patient",
    "clif_adt",
    "clif_labs",
    "clif_vitals",
    "clif_patient_assessments",
    "clif_medication_admin_continuous",
    "clif_medication_admin_intermittent",
    "clif_respiratory_support",
    "clif_microbiology_culture",
    "clif_hospital_diagnosis",
]


def prefilter_clif_tables(
    tables_path: str,
    file_type: str,
    hospitalization_ids: set,
    output_dir: Path,
    config: dict,
    patient_ids: set | None = None,
) -> str:
    """Filter CLIF tables to only the given hospitalization_ids (and patient_ids).

    Writes filtered parquet files to output_dir and a matching config.json.
    Returns the path to the filtered tables directory (as string).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    hosp_ids_str = {str(x) for x in hospitalization_ids}
    pat_ids_str = {str(x) for x in patient_ids} if patient_ids else None

    for tname in _CLIF_TABLES_FOR_SOFA_ASE:
        src = Path(tables_path) / f"{tname}.{file_type}"
        dst = output_dir / f"{tname}.parquet"
        if not src.exists():
            continue

        try:
            table = pq.read_table(src)
            df = table.to_pandas()
            orig_rows = len(df)

            if "hospitalization_id" in df.columns:
                df["hospitalization_id"] = df["hospitalization_id"].astype(str)
                df = df[df["hospitalization_id"].isin(hosp_ids_str)]
            elif "patient_id" in df.columns and pat_ids_str:
                df["patient_id"] = df["patient_id"].astype(str)
                df = df[df["patient_id"].isin(pat_ids_str)]

            df.to_parquet(dst, index=False)
            print(f"    {tname}: {orig_rows:,} -> {len(df):,} rows")
        except Exception as e:
            print(f"    {tname}: skipped ({e})")

    # Write filtered config pointing to the new directory
    filtered_config = config.copy()
    filtered_config["tables_path"] = str(output_dir)
    filtered_config["file_type"] = "parquet"
    config_path = output_dir / "config.json"
    with open(config_path, "w") as f:
        json.dump(filtered_config, f)

    return str(output_dir)


def get_tables_path(
    config: dict,
    hospitalization_ids: set,
    intermediate_dir: Path,
    patient_ids: set | None = None,
) -> str:
    """Pre-filter CLIF tables to the cohort's hospitalization_ids and patient_ids.

    Always filters to ensure consistent behavior across sites and prevent OOM.
    Reuses filtered tables from a previous run if available.
    Returns the path to the filtered tables directory.
    """
    tables_path = config["tables_path"]
    file_type = config["file_type"]
    filtered_dir = intermediate_dir / "filtered_clif"

    # Reuse if already filtered from a previous run
    if filtered_dir.exists() and (filtered_dir / "config.json").exists():
        print(f"  Using pre-filtered CLIF tables from {filtered_dir}")
        return str(filtered_dir)

    print(f"  Pre-filtering CLIF tables to {len(hospitalization_ids):,} cohort hospitalization_ids …")
    return prefilter_clif_tables(
        tables_path, file_type, hospitalization_ids, filtered_dir, config,
        patient_ids=patient_ids,
    )
