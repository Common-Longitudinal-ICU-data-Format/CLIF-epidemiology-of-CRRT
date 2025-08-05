"""
pyCLIF Utility Module for CRRT Epidemiology Study

This module provides utility functions for working with CLIF data.
Based on the pyCLIF module from the eligibility-for-mobilization project.
"""

import pandas as pd
import numpy as np
import json
import os
import shutil
import duckdb
import pytz
from datetime import datetime
from typing import Union, List, Dict, Optional, Tuple


# Initialize DuckDB connection
conn = duckdb.connect(database=':memory:')


def load_config():
    """
    Load configuration from config.json file.
    
    Returns:
        dict: Configuration dictionary
    """
    project_root = os.path.dirname(os.path.dirname(__file__))
    json_path = os.path.join(project_root, 'config', 'config.json')
    
    with open(json_path, 'r', encoding='utf-8') as file:
        config = json.load(file)
    print("Loaded configuration from config.json")
    return config


# Load configuration on module import
helper = load_config()


def _cast_id_cols_to_string(df: pd.DataFrame) -> pd.DataFrame:
    """
    Cast all ID columns (ending with '_id') to string type.
    
    Args:
        df: Input DataFrame
        
    Returns:
        DataFrame with ID columns cast to string
    """
    id_cols = [c for c in df.columns if c.endswith("_id")]
    if id_cols:
        df[id_cols] = df[id_cols].astype("string")
    return df


def load_parquet_with_tz(file_path: str, columns: Optional[List[str]] = None, 
                        filters: Optional[Dict] = None, sample_size: Optional[int] = None) -> pd.DataFrame:
    """
    Load parquet file with timezone support using DuckDB.
    
    Args:
        file_path: Path to parquet file
        columns: List of columns to load (None for all)
        filters: Dictionary of filters to apply
        sample_size: Number of rows to sample (None for all)
        
    Returns:
        pd.DataFrame: Loaded data
    """
    con = duckdb.connect()
    # DuckDB >=0.9 understands the original zone if we ask for TIMESTAMPTZ
    con.execute("SET timezone = 'UTC';")          # read & return in UTC
    con.execute("SET pandas_analyze_sample=0;")   # avoid sampling issues

    sel = "*" if columns is None else ", ".join(columns)
    query = f"SELECT {sel} FROM parquet_scan('{file_path}')"

    if filters:                                  # optional WHERE clause
        clauses = []
        for col, val in filters.items():
            if isinstance(val, list):
                vals = ", ".join([f"'{v}'" for v in val])
                clauses.append(f"{col} IN ({vals})")
            else:
                clauses.append(f"{col} = '{val}'")
        query += " WHERE " + " AND ".join(clauses)
    if sample_size:
        query += f" LIMIT {sample_size}"

    df = con.execute(query).fetchdf()            # pandas DataFrame
    con.close()
    df = _cast_id_cols_to_string(df)         # cast id columns to string
    return df


def load_data(table: str, sample_size: Optional[int] = None, 
              columns: Optional[List[str]] = None, filters: Optional[Dict] = None) -> pd.DataFrame:
    """
    Load data from a file in the specified directory with the option to select specific columns and apply filters.

    Parameters:
        table: The name of the table to load.
        sample_size: Number of rows to load.
        columns: List of column names to load.
        filters: Dictionary of filters to apply.

    Returns:
        pd.DataFrame: DataFrame containing the requested data.
    """
    # Determine the file path based on the directory and filetype
    file_name = f"{table}.{helper['file_type']}"
    file_path = os.path.join(helper['tables_path'], file_name)
    
    # Load the data based on filetype
    if os.path.exists(file_path):
        if helper['file_type'] == 'csv':
            # For CSV, we can use DuckDB to read specific columns and apply filters efficiently
            con = duckdb.connect()
            # Build the SELECT clause
            select_clause = "*" if not columns else ", ".join(columns)
            # Start building the query
            query = f"SELECT {select_clause} FROM read_csv_auto('{file_path}')"
            # Apply filters
            if filters:
                filter_clauses = []
                for column, values in filters.items():
                    if isinstance(values, list):
                        # Escape single quotes and wrap values in quotes
                        values_list = ', '.join(["'" + str(value).replace("'", "''") + "'" for value in values])
                        filter_clauses.append(f"{column} IN ({values_list})")
                    else:
                        value = str(values).replace("'", "''")
                        filter_clauses.append(f"{column} = '{value}'")
                if filter_clauses:
                    query += " WHERE " + " AND ".join(filter_clauses)
            # Apply sample size limit
            if sample_size:
                query += f" LIMIT {sample_size}"
            # Execute the query and fetch the data
            df = con.execute(query).fetchdf()
            con.close()
        elif helper['file_type'] == 'parquet':
            df = load_parquet_with_tz(file_path, columns, filters, sample_size)
        else:
            raise ValueError("Unsupported filetype. Only 'csv' and 'parquet' are supported.")
        print(f"Data loaded successfully from {file_path}")
        df = _cast_id_cols_to_string(df) # Cast id columns to string
        return df
    else:
        raise FileNotFoundError(f"The file {file_path} does not exist in the specified directory.")


def convert_datetime_columns_to_site_tz(df: pd.DataFrame, site_tz_str: str, verbose: bool = True) -> pd.DataFrame:
    """
    Convert all datetime columns in the DataFrame to the specified site timezone.

    Parameters:
        df: Input DataFrame.
        site_tz_str: Timezone string, e.g., "America/New_York" or "US/Central"
        verbose: Whether to print detailed output (default: True).

    Returns:
        pd.DataFrame: Modified DataFrame with datetime columns converted.
    """
    site_tz = pytz.timezone(site_tz_str)

    # Identify datetime-related columns
    dttm_columns = [col for col in df.columns if 'dttm' in col]

    for col in dttm_columns:
        df[col] = pd.to_datetime(df[col], errors='coerce')
        if pd.api.types.is_datetime64tz_dtype(df[col]):
            current_tz = df[col].dt.tz
            if current_tz == site_tz:
                if verbose:
                    print(f"{col}: Already in your timezone ({current_tz}), no conversion needed.")
            elif current_tz == pytz.UTC:
                print(f"{col}: null count before conversion= {df[col].isna().sum()}")
                df[col] = df[col].dt.tz_convert(site_tz)
                if verbose:
                    print(f"{col}: Converted from UTC to your timezone ({site_tz}).")
                    print(f"{col}: null count after conversion= {df[col].isna().sum()}")
            else:
                print(f"{col}: null count before conversion= {df[col].isna().sum()}")
                df[col] = df[col].dt.tz_convert(site_tz)
                if verbose:
                    print(f"{col}: Your timezone is {current_tz}, Converting to your site timezone ({site_tz}).")
                    print(f"{col}: null count after conversion= {df[col].isna().sum()}")
        elif pd.api.types.is_datetime64_any_dtype(df[col]):
            if verbose:
                df[col] = df[col].dt.tz_localize(site_tz, ambiguous=True, nonexistent='shift_forward')
                print(f"WARNING: {col}: Naive datetime, NOT converting. Assuming it's in your LOCAL ZONE. Please check ETL!")
        else:
            if verbose:
                print(f"WARNING: {col}: Not a datetime column. Please check ETL and run again!")
    return df


def count_unique_encounters(df: pd.DataFrame, encounter_column: str = 'hospitalization_id') -> int:
    """
    Counts the unique encounters in a DataFrame.
    
    Parameters:
        df: The DataFrame to analyze.
        encounter_column: The name of the column containing encounter IDs (default is 'hospitalization_id').
    
    Returns:
        int: The number of unique encounters.
    """
    return df[encounter_column].nunique()


def remove_duplicates(df: pd.DataFrame, subset: List[str], table_name: str) -> pd.DataFrame:
    """
    Remove duplicate rows from a DataFrame based on specified columns.
    
    Args:
        df: Input DataFrame
        subset: List of columns to consider for identifying duplicates
        table_name: Name of the table (for logging)
        
    Returns:
        pd.DataFrame: DataFrame with duplicates removed
    """
    original_count = len(df)
    df_clean = df.drop_duplicates(subset=subset, keep='first')
    removed_count = original_count - len(df_clean)
    
    if removed_count > 0:
        print(f"Removed {removed_count} duplicate rows from {table_name} table")
    
    return df_clean


def create_summary_table(df: pd.DataFrame, numeric_col: str, group_by_cols: Union[str, List[str]] = None) -> pd.DataFrame:
    """
    Create a summary statistics table for a numeric column, optionally grouped by other columns.
    
    Args:
        df: Input DataFrame
        numeric_col: Name of the numeric column to summarize
        group_by_cols: Column(s) to group by (optional)
        
    Returns:
        pd.DataFrame: Summary statistics table
    """
    if group_by_cols is None:
        # Overall summary
        summary = pd.DataFrame({
            'N': [df[numeric_col].count()],
            'missing': [df[numeric_col].isna().sum()],
            'min': [df[numeric_col].min()],
            'q25': [df[numeric_col].quantile(0.25)],
            'median': [df[numeric_col].median()],
            'q75': [df[numeric_col].quantile(0.75)],
            'mean': [df[numeric_col].mean()],
            'max': [df[numeric_col].max()]
        })
    else:
        # Grouped summary
        if isinstance(group_by_cols, str):
            group_by_cols = [group_by_cols]
            
        summary = df.groupby(group_by_cols)[numeric_col].agg([
            ('N', 'count'),
            ('missing', lambda x: x.isna().sum()),
            ('min', 'min'),
            ('q25', lambda x: x.quantile(0.25)),
            ('median', 'median'),
            ('q75', lambda x: x.quantile(0.75)),
            ('mean', 'mean'),
            ('max', 'max')
        ]).reset_index()
    
    return summary.round(2)


def apply_outlier_thresholds(df: pd.DataFrame, column: str, min_val: float, max_val: float) -> pd.DataFrame:
    """
    Apply outlier thresholds to a column by setting values outside the range to NaN.
    
    Args:
        df: Input DataFrame
        column: Column name to apply thresholds to
        min_val: Minimum valid value
        max_val: Maximum valid value
        
    Returns:
        pd.DataFrame: DataFrame with outliers set to NaN
    """
    if column in df.columns:
        original_count = df[column].notna().sum()
        df.loc[df[column] < min_val, column] = np.nan
        df.loc[df[column] > max_val, column] = np.nan
        outlier_count = original_count - df[column].notna().sum()
        
        if outlier_count > 0:
            print(f"Set {outlier_count} outlier values in {column} to NaN")
    
    return df


def setup_output_folders():
    """
    Set up output folder structure, backing up existing output if present.
    
    Returns:
        str: Path to output folder
    """
    print("=== Output Folder Management ===")
    
    output_folder = '../output'
    output_old_folder = '../output_old'
    
    # Check if output folder exists
    if os.path.exists(output_folder):
        print(f"Existing output folder found: {output_folder}")
        
        # If output_old already exists, remove it first
        if os.path.exists(output_old_folder):
            print(f"Removing existing output_old folder...")
            shutil.rmtree(output_old_folder)
        
        # Rename current output to output_old
        print(f"Renaming {output_folder} → {output_old_folder}")
        os.rename(output_folder, output_old_folder)
        
        # Log what was backed up
        if os.path.exists(output_old_folder):
            backup_size = sum(
                os.path.getsize(os.path.join(dirpath, filename))
                for dirpath, dirnames, filenames in os.walk(output_old_folder)
                for filename in filenames
            ) / (1024 * 1024)  # Convert to MB
            print(f"Backup created: {backup_size:.1f} MB")
    
    # Create fresh output directory structure
    print(f"Creating fresh output directory structure...")
    os.makedirs(output_folder, exist_ok=True)
    os.makedirs(f'{output_folder}/final', exist_ok=True)
    os.makedirs(f'{output_folder}/intermediate', exist_ok=True)
    os.makedirs(f'{output_folder}/final/graphs', exist_ok=True)
    
    print(f"Output directory structure ready:")
    print(f"   {output_folder}/")
    print(f"   ├── final/")
    print(f"   │   └── graphs/")
    print(f"   └── intermediate/")
    
    return output_folder



## meds dose conversion helpers

# Define medications and their unit conversion information
meds_list = [
    "norepinephrine", "epinephrine", "phenylephrine",
    "vasopressin", "dopamine", "angiotensin", "metaraminol", "dobutamine"
]

med_unit_info = {
    'norepinephrine': {
        'required_unit': 'mcg/kg/min',
        'acceptable_units': ['mcg/kg/min', 'mcg/kg/hr', 'mg/kg/hr', 'mcg/min', 'mg/hr'],
    },
    'epinephrine': {
        'required_unit': 'mcg/kg/min',
        'acceptable_units': ['mcg/kg/min', 'mcg/kg/hr', 'mg/kg/hr', 'mcg/min', 'mg/hr'],
    },
    'phenylephrine': {
        'required_unit': 'mcg/kg/min',
        'acceptable_units': ['mcg/kg/min', 'mcg/kg/hr', 'mg/kg/hr', 'mcg/min', 'mg/hr'],
    },
    'dopamine': {
        'required_unit': 'mcg/kg/min',
        'acceptable_units': ['mcg/kg/min', 'mcg/kg/hr', 'mg/kg/hr', 'mcg/min', 'mg/hr'],
    },
    'dobutamine': {
        'required_unit': 'mcg/kg/min',
        'acceptable_units': ['mcg/kg/min', 'mcg/kg/hr', 'mg/kg/hr', 'mcg/min', 'mg/hr'],
    },
    'metaraminol': {
        'required_unit': 'mcg/kg/min',
        'acceptable_units': ['mg/hr', 'mcg/min'],
    },
    'angiotensin': {
        'required_unit': 'mcg/kg/min',
        'acceptable_units': ['ng/kg/min', 'ng/kg/hr'],
    },
    'vasopressin': {
        'required_unit': 'units/min',
        'acceptable_units': ['units/min', 'units/hr', 'units/hour','milliunits/min', 'milliunits/hr', 'milli-units/kg/hr',
                             'milli-units/min', 'milli-units/hr','milli-units/kg/h','milliunits/kg/h',
                             'milli-units/kg/min', 'milliunits/kg/min' ],
    },
}

def check_dose_unit(row):
    med_category = row['med_category']
    med_dose_unit = row['med_dose_unit']
    # Check if med_category exists in med_unit_info
    if med_category in med_unit_info:
        # Check if med_dose_unit is in the acceptable units
        if med_dose_unit in med_unit_info[med_category]['acceptable_units']:
            return "Valid"
        else:
            return "Not an acceptable unit"
    else:
        return "Not a vasoactive"
    
def has_per_hour_or_min(unit):
    if pd.isnull(unit):
        return False
    unit = unit.lower()
    return '/hr' in unit or '/min' in unit or '/hour' in unit

def is_dose_within_range(row, outlier_dict):
    '''
    Check if med_dose_converted is within the outlier-configured range for this med_category.
    Parameters:
        row (pd.Series): A row from a DataFrame, must include 'med_category' and 'med_dose_converted'.
        outlier_dict (dict): Dictionary of min/max pairs from outlier_config.json.
    Returns:
        bool: True if the dose is within range or if med_category is not found, False otherwise.
    '''
    med_category = row['med_category']
    med_dose_converted = row['med_dose_converted']
    dose_range = outlier_dict.get(med_category, None)
    if dose_range is None:
        return False
    min_dose, max_dose = dose_range
    return min_dose <= med_dose_converted <= max_dose

def get_conversion_factor(med_category: str,
                          med_dose_unit: str,
                          weight_kg: Union[float, None]) -> Union[float, None]:
    """
    Return a multiplier that converts *med_dose* from its current unit
    to the **required** unit for that medication.

    For units that need weight (mcg/min  → mcg/kg/min, mg/hr → mcg/kg/min …)
    we return *None* when weight_kg is missing so that the caller can decide
    what to do (keep the row with NaN, drop it, etc.).
    """
    med_info = med_unit_info.get(med_category)
    if med_info is None:
        return None                          # not a vaso we care about

    med_dose_unit = (med_dose_unit or "").lower()

    # ── helpers ───────────────────────────────────────────────
    has_weight = weight_kg is not None and not pd.isna(weight_kg)
    def w_needed(factor_if_known):
        return factor_if_known if has_weight else None
    # ──────────────────────────────────────────────────────────

    if med_category in ["norepinephrine", "epinephrine",
                        "phenylephrine", "dopamine",
                        "dobutamine", "metaraminol"]:
        if med_dose_unit == "mcg/kg/min": return 1.0
        elif med_dose_unit == "mcg/kg/hr": return 1/60
        elif med_dose_unit == "mg/kg/hr": return 1000/60
        elif med_dose_unit == "mg/kg/min": return 1000
        elif med_dose_unit == "mcg/min": return w_needed(1/weight_kg)
        elif med_dose_unit == "mg/hr": return w_needed(1000/60/weight_kg)
    elif med_category == "angiotensin":
        if med_dose_unit == "ng/kg/min": return 1/1_000
        elif med_dose_unit == "ng/kg/hr": return 1/1_000/60
        elif med_dose_unit == "mcg/kg/min": return 1.0
    elif med_category == "vasopressin":
        if med_dose_unit == "units/min": return 1.0
        elif med_dose_unit in ["units/hr", "units/hour"]: return 1/60
        elif med_dose_unit == "milliunits/min": return 1/1_000
        elif med_dose_unit == "milli-units/min": return 1/1_000
        elif med_dose_unit == "milliunits/hr": return 1/1_000/60
        elif med_dose_unit == "milli-units/hr": return 1/1_000/60
        elif med_dose_unit == "milliunits/kg/hr": return w_needed(1/1_000/60/weight_kg)
        elif med_dose_unit == "milli-units/kg/hr": return w_needed(1/1_000/60/weight_kg)
        elif med_dose_unit == "milliunits/kg/min": return w_needed(1/1_000/weight_kg)
        elif med_dose_unit == "milli-units/kg/min": return w_needed(1/1_000/weight_kg)
    return None                               # unit not recognised


def convert_dose(row: pd.Series) -> Union[float, None]:
    """
    Convert `row.med_dose` to the medication's **required** unit.

    Returns `np.nan` when:
    * the unit is unrecognised
    * weight is required but missing
    """
    factor = get_conversion_factor(
        row["med_category"],
        row["med_dose_unit"],
        row["weight_kg"]
    )
    return np.nan if factor is None else round(row["med_dose"] * factor, 5)


def stitch_encounters(hospitalization, adt, time_interval=6):
    """
    Stitches together related hospital encounters that occur within a specified time interval of each other.
    
    Args:
        hospitalization (pd.DataFrame): Hospitalization table with required columns
        adt (pd.DataFrame): ADT table with required columns
        time_interval (int, optional): Number of hours between encounters to consider them linked. Defaults to 6.
        
    Returns:
        pd.DataFrame: Stitched encounters with encounter blocks
    """
    hospitalization_filtered = hospitalization[["patient_id","hospitalization_id","admission_dttm",
                                                "discharge_dttm","age_at_admission", "admission_type_category", "discharge_category"]].copy()
    hospitalization_filtered['admission_dttm'] = pd.to_datetime(hospitalization_filtered['admission_dttm'])
    hospitalization_filtered['discharge_dttm'] = pd.to_datetime(hospitalization_filtered['discharge_dttm'])

    hosp_adt_join = pd.merge(hospitalization_filtered[["patient_id","hospitalization_id","age_at_admission","admission_type_category",
                                                       "admission_dttm","discharge_dttm",
                                                        "discharge_category"]], 
                      adt[["hospitalization_id","in_dttm","out_dttm","location_category","hospital_id"]],
                 on="hospitalization_id",how="left")

    hospital_cat = hosp_adt_join[["hospitalization_id","in_dttm","out_dttm","hospital_id"]]

    # Step 1: Sort by patient_id and admission_dttm
    hospital_block = hosp_adt_join[["patient_id","hospitalization_id","admission_dttm","discharge_dttm", "age_at_admission",  "discharge_category", "admission_type_category"]]
    hospital_block = hospital_block.drop_duplicates()
    hospital_block = hospital_block.sort_values(by=["patient_id", "admission_dttm"]).reset_index(drop=True)
    hospital_block = hospital_block[["patient_id","hospitalization_id","admission_dttm","discharge_dttm", "age_at_admission",  "discharge_category", "admission_type_category"]]

    # Step 2: Calculate time between discharge and next admission
    hospital_block["next_admission_dttm"] = hospital_block.groupby("patient_id")["admission_dttm"].shift(-1)
    hospital_block["discharge_to_next_admission_hrs"] = (
        (hospital_block["next_admission_dttm"] - hospital_block["discharge_dttm"]).dt.total_seconds() / 3600
    )

    # Step 3: Create linked column based on time_interval
    hospital_block["linked6hrs"] = hospital_block["discharge_to_next_admission_hrs"] < time_interval

    # Sort values to ensure correct order
    hospital_block = hospital_block.sort_values(by=["patient_id", "admission_dttm"]).reset_index(drop=True)

    # Initialize encounter_block with row indices + 1
    hospital_block['encounter_block'] = hospital_block.index + 1

    # Iteratively propagate the encounter_block values
    while True:
        shifted = hospital_block['encounter_block'].shift(-1)
        mask = hospital_block['linked6hrs'] & (hospital_block['patient_id'] == hospital_block['patient_id'].shift(-1))
        hospital_block.loc[mask, 'encounter_block'] = shifted[mask]
        if hospital_block['encounter_block'].equals(hospital_block['encounter_block'].bfill()):
            break

    hospital_block['encounter_block'] = hospital_block['encounter_block'].bfill(downcast='int')
    hospital_block = pd.merge(hospital_block,hospital_cat,how="left",on="hospitalization_id")
    hospital_block = hospital_block.sort_values(by=["patient_id", "admission_dttm","in_dttm","out_dttm"]).reset_index(drop=True)
    hospital_block = hospital_block.drop_duplicates()

    hospital_block2 = hospital_block.groupby(['patient_id','encounter_block']).agg(
        admission_dttm=pd.NamedAgg(column='admission_dttm', aggfunc='min'),
        discharge_dttm=pd.NamedAgg(column='discharge_dttm', aggfunc='max'),
        admission_type_category=pd.NamedAgg(column='admission_type_category', aggfunc='first'),
        discharge_category=pd.NamedAgg(column='discharge_category', aggfunc='last'),
        hospital_id = pd.NamedAgg(column='hospital_id', aggfunc='last'),
        age_at_admission=pd.NamedAgg(column='age_at_admission', aggfunc='last'),
        list_hospitalization_id=pd.NamedAgg(column='hospitalization_id', aggfunc=lambda x: sorted(x.unique()))
    ).reset_index()

    df = pd.merge(hospital_block[["patient_id",
                                  "hospitalization_id",
                                  "encounter_block"]].drop_duplicates(),
             hosp_adt_join[["hospitalization_id","location_category","in_dttm","out_dttm"]], on="hospitalization_id",how="left")

    df = pd.merge(df,hospital_block2[["encounter_block",
                                      "admission_dttm",
                                      "discharge_dttm",
                                      "discharge_category",
                                      "admission_type_category",
                                      "age_at_admission",
                                      "hospital_id",
                                     "list_hospitalization_id"]],on="encounter_block",how="left")
    df = df.drop_duplicates(subset=["patient_id","encounter_block","in_dttm","out_dttm","location_category"])
    
    return df


def categorize_device(row):
    if pd.notna(row['device_category']):
        return row['device_category']
    elif row['mode_category'] in ["simv", "pressure-regulated volume control", "assist control-volume control"]:
        return "vent"
    elif pd.isna(row['device_category']) and row['fio2_set'] == 0.21 and pd.isna(row['lpm_set']) and pd.isna(row['peep_set']) and pd.isna(row['tidal_volume_set']):
        return "room air"
    elif pd.isna(row['device_category']) and pd.isna(row['fio2_set']) and row['lpm_set'] == 0 and pd.isna(row['peep_set']) and pd.isna(row['tidal_volume_set']):
        return "room air"
    elif pd.isna(row['device_category']) and pd.isna(row['fio2_set']) and (0 < row['lpm_set'] <= 20) and pd.isna(row['peep_set']) and pd.isna(row['tidal_volume_set']):
        return "nasal cannula"
    elif pd.isna(row['device_category']) and pd.isna(row['fio2_set']) and row['lpm_set'] > 20 and pd.isna(row['peep_set']) and pd.isna(row['tidal_volume_set']):
        return "high flow nc"
    elif row['device_category'] == "nasal cannula" and pd.isna(row['fio2_set']) and row['lpm_set'] > 20:
        return "high flow nc"
    else:
        return row['device_category']  # Keep original value if no condition is met

# Try to fill in FiO2 based on other values    
def refill_fio2(row):
    if pd.notna(row['fio2_set']):
        return row['fio2_set']/100
    elif pd.isna(row['fio2_set']) and row['device_category'] == "room air":
        return 0.21 
    elif pd.isna(row['fio2_set']) and row['device_category'] == "nasal cannula" and pd.notna(row['lpm_set']):
        return (0.24 + (0.04 * row['lpm_set'])) 
    else:
        return np.nan