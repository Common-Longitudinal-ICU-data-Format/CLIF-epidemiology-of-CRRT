"""
Shared helpers for the CRRT epidemiology pipeline.

Provides:
  - validate_config()      : validate config.json fields
  - load_intermediate()    : read parquet with clear missing-file errors
  - safe_load_clif_table() : wrap clifpy table loading with clear errors
"""

import os
from pathlib import Path

import pandas as pd

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
