#!/usr/bin/env python
# -*- coding: utf-8 -*-


################################################################################
# CLIF CRRT Epidemiology — Build Wide Dataset
################################################################################

################################################################################
# Setup
# * Load libraries, create required directories, and load the config file
################################################################################

import sys
import os
from pathlib import Path
# Add project root to sys.path FIRST so the utils/ package is found
# before code/utils.py (which Python auto-adds to sys.path)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import polars as pl
import matplotlib.pyplot as plt
import pandas as pd
import clifpy
from utils.config import config
from clifpy.utils.outlier_handler import apply_outlier_handling
import gc
import yaml
import numpy as np
from clifpy.utils.unit_converter import convert_dose_units_by_med_category

site_name = config['site_name']
tables_path = config['tables_path']
file_type = config['file_type']
project_root = config['project_root']
has_crrt_settings = config.get('has_crrt_settings', True)
sys.path.insert(0, project_root)
print(f"Site Name: {site_name}")
print(f"Tables Path: {tables_path}")
print(f"File Type: {file_type}")
PROJECT_ROOT = Path(config['project_root'])
UTILS_DIR = PROJECT_ROOT / "utils"
OUTPUT_DIR = PROJECT_ROOT / "output"
OUTPUT_FINAL_DIR = OUTPUT_DIR / "final"
OUTPUT_INTERMEDIATE_DIR = OUTPUT_DIR / "intermediate"

# Load config - fix path to work from any directory
script_dir = Path(__file__).parent
config_path = script_dir.parent / "config" / "clif_data_requirements.yaml"

if config_path.exists():
    with open(config_path, 'r') as f:
        clif_config = yaml.safe_load(f)
else:
    # If config file doesn't exist, use default configuration
    print(f"Warning: {config_path} not found. Using default configuration.")
    clif_config = {
        'data_requirements': {
            'labs': True,
            'vitals': True,
            'medications': True,
            'respiratory': True
        },
        'patient_required_columns': [
            'patient_id', 'sex_category', 'race_category',
            'ethnicity_category', 'date_of_birth'
        ],
        'hospitalization_required_columns': [
            'hospitalization_id', 'patient_id', 'admission_dttm',
            'discharge_dttm', 'death_dttm', 'age_at_admission'
        ],
        'labs_required_columns': [
            'hospitalization_id', 'lab_result_dttm',
            'lab_category', 'lab_value_numeric'
        ],
        'vitals_required_columns': [
            'hospitalization_id', 'recorded_dttm',
            'vital_category', 'vital_value'
        ],
        'meds_continuous_required_columns': [
            'hospitalization_id', 'admin_dttm',
            'med_category', 'med_dose'
        ],
        'respiratory_support_required_columns': [
            'hospitalization_id', 'recorded_dttm',
            'ventilator_mode', 'device_category'
        ],
        'adt_required_columns': [
            'hospitalization_id', 'in_dttm', 'out_dttm',
            'location_category', 'location_name'
        ],
        'labs_of_interest': [
            'hemoglobin', 'platelet_count', 'white_blood_cell_count',
            'creatinine', 'blood_urea_nitrogen', 'glucose',
            'sodium', 'potassium', 'chloride', 'bicarbonate',
            'bilirubin_total', 'albumin', 'lactate'
        ],
        'vitals_of_interest': [
            'heart_rate', 'respiratory_rate', 'temperature',
            'sbp', 'dbp', 'map', 'spo2'
        ],
        'meds_continuous_of_interest': [
            'norepinephrine', 'epinephrine', 'dopamine',
            'dobutamine', 'vasopressin', 'milrinone'
        ]
    }

################################################################################
# Load Cohort
################################################################################

cohort_df_fp = str(OUTPUT_INTERMEDIATE_DIR / "outcomes_df.parquet")
cohort_df = pd.read_parquet(cohort_df_fp)
# Ensure correct dtypes for join keys
cohort_df['hospitalization_id'] = cohort_df['hospitalization_id'].astype(str)
cohort_df['encounter_block'] = cohort_df['encounter_block'].astype('int32')

# Get all patient_ids and hospitalization_ids for filtering
final_patient_ids = cohort_df['patient_id'].tolist()
final_hosp_ids = cohort_df['hospitalization_id'].tolist()

################################################################################
# CLIF Core
# * Patient, Hospitalization
################################################################################

print("\n" + "=" * 80)
print("Loading CLIF Tables")
print("=" * 80)

from clifpy.clif_orchestrator import ClifOrchestrator
# Initialize ClifOrchestrator
clif = ClifOrchestrator(
    data_directory=config['tables_path'],
    filetype=config['file_type'],
    timezone=config['timezone'],
    output_directory = OUTPUT_DIR
)

clif.load_table(
        'patient',
        filters={'patient_id': list(final_patient_ids)}
    )

clif.load_table(
        'hospitalization',
        filters={'hospitalization_id': list(final_hosp_ids)}
    )

patient_required_columns = clif_config['patient_required_columns']
hosp_required_columns = clif_config['hospitalization_required_columns']

patient_df = clif.patient.df.loc[:, patient_required_columns]
hosp_df = clif.hospitalization.df.loc[:, hosp_required_columns]

# Merge patient and hospitalization tables on 'patient_id'
static_df = hosp_df.merge(
    patient_df,
    on='patient_id',
    how='left'
)

# cleanup
print(clif.get_tables_obj_list())
# Delete the patient and hospitalization table objects
clif.patient = None
clif.hospitalization = None
del patient_df, hosp_df
# Verify they're gone
print(clif.get_tables_obj_list())

# Merge static_df with cohort_df to add cohort-level columns
static_df = static_df.merge(
    cohort_df,
    on='hospitalization_id',
    how='left'
)

################################################################################
# Wide Dataset
################################################################################

# Define columns for each table
columns_config = {
    'labs': clif_config['labs_required_columns'],
    'vitals': clif_config['vitals_required_columns'],
    'medication_admin_continuous': clif_config['meds_continuous_required_columns'],
}

# Define filters for each table (hospitalization_id + category filters)
filters_config = {
    'labs': {
        'hospitalization_id': list(final_hosp_ids),
        'lab_category': clif_config['labs_of_interest']
    },
    'vitals': {
        'hospitalization_id': list(final_hosp_ids),
        'vital_category': clif_config['vitals_of_interest']
    },
    'medication_admin_continuous': {
        'hospitalization_id': list(final_hosp_ids),
        'med_category': clif_config['meds_continuous_of_interest']
    },
}

# Load all tables in one call
clif.initialize(
    tables=['labs', 'vitals', 'medication_admin_continuous'],
    columns=columns_config,
    filters=filters_config
)
print("clifpy version:", getattr(clifpy, '__version__', 'unknown'))
apply_outlier_handling(clif.labs)
apply_outlier_handling(clif.vitals)

# =============================================================================
# Extract DataFrames from loaded tables
# =============================================================================
labs_df = clif.labs.df.copy()
vitals_df = clif.vitals.df.copy()
meds_cont_df = clif.medication_admin_continuous.df.copy()

# =============================================================================
# Convert vasopressor units to standard doses using clifpy
# =============================================================================
vaso_categories = ['norepinephrine', 'epinephrine', 'phenylephrine',
                   'dopamine', 'vasopressin', 'angiotensin']
vaso_preferred_units = {
    'dopamine': 'mcg/kg/min',
    'norepinephrine': 'mcg/kg/min',
    'epinephrine': 'mcg/kg/min',
    'phenylephrine': 'mcg/kg/min',
    'angiotensin': 'mcg/kg/min',
    'vasopressin': 'u/min',
}

vaso_mask = meds_cont_df['med_category'].isin(vaso_categories)
vaso_df = meds_cont_df[vaso_mask].copy()
non_vaso_df = meds_cont_df[~vaso_mask].copy()

if len(vaso_df) > 0:
    vaso_converted, _ = convert_dose_units_by_med_category(
        vaso_df, vitals_df=vitals_df,
        preferred_units=vaso_preferred_units, override=True
    )
    vaso_converted['med_dose'] = vaso_converted['med_dose_converted']
    meds_cont_df = pd.concat([vaso_converted, non_vaso_df], ignore_index=True)
    print(f"Vasopressor unit conversion complete: {len(vaso_converted)} rows converted")
else:
    print("No vasopressor rows found, skipping unit conversion")

# =============================================================================
# Pivot Vitals: narrow → wide
# =============================================================================
vitals_wide = vitals_df.pivot_table(
    index=['hospitalization_id', 'recorded_dttm'],
    columns='vital_category',
    values='vital_value',
    aggfunc='first'
).reset_index()

# Rename datetime column to event_dttm and flatten column names
vitals_wide = vitals_wide.rename(columns={'recorded_dttm': 'event_dttm'})
vitals_wide.columns = ['hospitalization_id', 'event_dttm'] + \
    [f'vital_{col}' for col in vitals_wide.columns[2:]]

# =============================================================================
# Pivot Labs: narrow → wide
# =============================================================================
labs_wide = labs_df.pivot_table(
    index=['hospitalization_id', 'lab_result_dttm'],
    columns='lab_category',
    values='lab_value_numeric',
    aggfunc='first'
).reset_index()

# Rename datetime column to event_dttm
labs_wide = labs_wide.rename(columns={'lab_result_dttm': 'event_dttm'})
labs_wide.columns = ['hospitalization_id', 'event_dttm'] + \
    [f'lab_{col}' for col in labs_wide.columns[2:]]

# =============================================================================
# Pivot Meds Continuous: narrow → wide
# =============================================================================
meds_cont_wide = meds_cont_df.pivot_table(
    index=['hospitalization_id', 'admin_dttm'],
    columns='med_category',
    values='med_dose',
    aggfunc='first'
).reset_index()

# Rename datetime column to event_dttm
meds_cont_wide = meds_cont_wide.rename(columns={'admin_dttm': 'event_dttm'})
if not pd.api.types.is_datetime64_any_dtype(meds_cont_wide['event_dttm']):
    meds_cont_wide['event_dttm'] = pd.to_datetime(meds_cont_wide['event_dttm'], utc=True)
meds_cont_wide.columns = ['hospitalization_id', 'event_dttm'] + \
    [f'med_cont_{col}' for col in meds_cont_wide.columns[2:]]

# =============================================================================
# Compute Norepinephrine Equivalent (NEE) from unit-converted vasopressors
# =============================================================================
nee_cols = {
    'med_cont_norepinephrine': 1.0,
    'med_cont_epinephrine': 1.0,
    'med_cont_phenylephrine': 1/10.0,
    'med_cont_dopamine': 1/100.0,
    'med_cont_vasopressin': 2.5,
    'med_cont_angiotensin': 10.0,
}
nee_terms = [meds_cont_wide[c].fillna(0) * w
             for c, w in nee_cols.items() if c in meds_cont_wide.columns]
if nee_terms:
    meds_cont_wide['med_cont_nee'] = sum(nee_terms, pd.Series(0.0, index=meds_cont_wide.index))
else:
    meds_cont_wide['med_cont_nee'] = np.nan
# Set NEE to NaN when no vasopressors were recorded at all
all_vaso_null = meds_cont_wide[[c for c in nee_cols if c in meds_cont_wide.columns]].isna().all(axis=1)
meds_cont_wide.loc[all_vaso_null, 'med_cont_nee'] = np.nan
print(f"NEE computed: {meds_cont_wide['med_cont_nee'].notna().sum()} non-null values")

# =============================================================================
# Merge all wide tables on hospitalization_id + event_dttm
# =============================================================================
wide_df = vitals_wide.merge(
    labs_wide,
    on=['hospitalization_id', 'event_dttm'],
    how='outer'
)

wide_df = wide_df.merge(
    meds_cont_wide,
    on=['hospitalization_id', 'event_dttm'],
    how='outer'
)

# Sort by hospitalization and time
wide_df = wide_df.sort_values(['hospitalization_id', 'event_dttm']).reset_index(drop=True)

print(f"Final wide dataframe: {wide_df.shape}")
print(f"Columns: {list(wide_df.columns)}")

del vitals_wide, labs_wide, meds_cont_wide, labs_df, vitals_df, meds_cont_df, vaso_df, non_vaso_df

################################################################################
# Respiratory Support
################################################################################

clif.load_table(
        'respiratory_support',
        columns= clif_config['respiratory_support_required_columns'], 
        filters={'hospitalization_id': list(final_hosp_ids)}
    )

clif.respiratory_support = clif.respiratory_support.waterfall()

# =============================================================================
# Respiratory Support (already wide format)
# =============================================================================
resp_df = clif.respiratory_support.df.copy()

# Rename datetime column to event_dttm
resp_df = resp_df.rename(columns={'recorded_dttm': 'event_dttm'})

# Add prefix to all columns except hospitalization_id and event_dttm
resp_cols_to_rename = [col for col in resp_df.columns if col not in ['hospitalization_id', 'event_dttm']]
resp_df = resp_df.rename(columns={col: f'resp_{col}' for col in resp_cols_to_rename})

# Merge with wide_df
wide_df = wide_df.merge(
    resp_df,
    on=['hospitalization_id', 'event_dttm'],
    how='outer'
)

# Re-sort by hospitalization and time
wide_df = wide_df.sort_values(['hospitalization_id', 'event_dttm']).reset_index(drop=True)

print(wide_df.columns)

################################################################################
# CRRT Therapy
################################################################################

if has_crrt_settings:
    crrt_cols_to_load = clif_config['crrt_required_columns']
else:
    crrt_cols_to_load = clif_config.get('crrt_minimal_columns', ['hospitalization_id', 'recorded_dttm'])

clif.load_table(
        'crrt_therapy',
        columns=crrt_cols_to_load,
        filters={'hospitalization_id': list(final_hosp_ids)}
    )

# =============================================================================
# CRRT Therapy (already wide format)
# =============================================================================
crrt_therapy_df = clif.crrt_therapy.df.copy()

# Rename datetime column to event_dttm
crrt_therapy_df = crrt_therapy_df.rename(columns={'recorded_dttm': 'event_dttm'})

# Add crrt_ prefix to setting columns that don't already have it
crrt_cols_to_rename = [col for col in crrt_therapy_df.columns
                       if col not in ['hospitalization_id', 'event_dttm']
                       and not col.startswith('crrt_')]
crrt_therapy_df = crrt_therapy_df.rename(columns={col: f'crrt_{col}' for col in crrt_cols_to_rename})

# Merge with wide_df
wide_df = wide_df.merge(
    crrt_therapy_df,
    on=['hospitalization_id', 'event_dttm'],
    how='outer'
)

# Re-sort by hospitalization and time
wide_df = wide_df.sort_values(['hospitalization_id', 'event_dttm']).reset_index(drop=True)

print(wide_df.columns)

################################################################################
# ADT
################################################################################

clif.load_table(
        'adt',
        columns= clif_config['adt_required_columns'], 
        filters={'hospitalization_id': list(final_hosp_ids)}
    )
adt_df = clif.adt.df

# =============================================================================
# Add ADT rows to wide_df (location transitions as events)
# =============================================================================
adt_df['in_dttm'] = pd.to_datetime(adt_df['in_dttm'])
adt_df['out_dttm'] = pd.to_datetime(adt_df['out_dttm'])

# Sort ADT by hospitalization and in_dttm
adt_sorted = adt_df.sort_values(['hospitalization_id', 'in_dttm'])

# Create ADT wide format with in_dttm as event_dttm
adt_wide = adt_df[['hospitalization_id', 'in_dttm', 'location_category', 'location_name', 'location_type']].copy()
adt_wide = adt_wide.rename(columns={'in_dttm': 'event_dttm'})

# Add prefix to ADT columns
adt_wide = adt_wide.rename(columns={
    'location_category': 'adt_location_category',
    'location_name': 'adt_location_name',
    'location_type': 'adt_location_type'
})

# Merge with wide_df
wide_df = wide_df.merge(
    adt_wide,
    on=['hospitalization_id', 'event_dttm'],
    how='outer'
)

# Sort by hospitalization and time
wide_df = wide_df.sort_values(['hospitalization_id', 'event_dttm']).reset_index(drop=True)

# Forward fill location within each hospitalization so every row knows the current location
wide_df['adt_location_category'] = wide_df.groupby('hospitalization_id')['adt_location_category'].ffill()
wide_df['adt_location_name'] = wide_df.groupby('hospitalization_id')['adt_location_name'].ffill()
wide_df['adt_location_type'] = wide_df.groupby('hospitalization_id')['adt_location_type'].ffill()

print(f"wide_df shape: {wide_df.shape}")

output_path = os.path.join(OUTPUT_INTERMEDIATE_DIR, f"wide_df.parquet")
wide_df.to_parquet(output_path, index=False)

################################################################################
# Save Final Wide Dataset
################################################################################
