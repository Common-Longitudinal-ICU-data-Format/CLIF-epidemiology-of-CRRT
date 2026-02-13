#!/usr/bin/env python
# coding: utf-8

# # Epidemiology of CRRT
# 
# Author: Kaveri Chhikara
# 
# This script identifies the cohort using CLIF 2.1 tables
# **Requirements**
# 
# * Required table filenames should be `clif_patient`, `clif_hospitalization`, `clif_adt`, `clif_vitals`, `clif_labs`, `clif_medication_admin_continuous`, `clif_respiratory_support` ,`crrt_therapy`, `clif_hospital_diagnosis`
# * Within each table, the following variables and categories are required.
# 
# | Table Name | Required Variables | Required Categories |
# | --- | --- | --- |
# | **clif_patient** | `patient_id`, `race_category`, `ethnicity_category`, `sex_category`, `death_dttm` | - |
# | **clif_hospitalization** | `patient_id`, `hospitalization_id`, `admission_dttm`, `discharge_dttm`, `age_at_admission`, `discharge_category` | - |
# | **clif_adt** |  `hospitalization_id`, `hospital_id`,`in_dttm`, `out_dttm`, `location_category`, `location_type` | - |
# | **clif_vitals** | `hospitalization_id`, `recorded_dttm`, `vital_category`, `vital_value` | heart_rate, resp_rate, sbp, dbp, map, spo2, weight_kg, height_cm |
# | **clif_labs** | `hospitalization_id`, `lab_result_dttm`, `lab_category`, `lab_value` | sodium, potassium, chloride, bicarbonate, bun, creatinine, glucose_serum, calcium_total, lactate, magnesium, ph_arterial, ph_venous, po2_arterial |
# | **clif_medication_admin_continuous** | `hospitalization_id`, `admin_dttm`, `med_name`, `med_category`, `med_dose`, `med_dose_unit` | norepinephrine, epinephrine, phenylephrine, vasopressin, dopamine, angiotensin, dobutamine, milrinone, isoproterenol |
# | **clif_respiratory_support** | `hospitalization_id`, `recorded_dttm`, `device_category`, `mode_category`, `tracheostomy`, `fio2_set`, `lpm_set`, `resp_rate_set`, `peep_set`, `resp_rate_obs`, `tidal_volume_set`, `pressure_control_set`, `pressure_support_set`, `peak_inspiratory_pressure_set`, `tidal_volume_obs` | - |
# | **clif_crrt_therapy** | `hospitalization_id`, `recorded_dttm`, `crrt_mode_name`, `crrt_mode_category`, `device_id`, `blood_flow_rate`, `dialysate_flow_rate`, `pre_filter_replacement_fluid_rate`,`post_filter_replacement_fluid_rate`, `ultrafilteration_out` | - |
# | **clif_hospital_diagnosis** | `hospitalization_id`, `diagnosis_code`, `present_on_admission` | - |

# # Setup

# In[1]:


import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for script mode
import matplotlib.pyplot as plt
import gc
from pathlib import Path
import json
import pyarrow
import warnings
import clifpy
from typing import Union
from tqdm import tqdm

import sys
import clifpy
import os

print("=== Environment Verification ===")
print(f"Python executable: {sys.executable}")
print(f"Python version: {sys.version}")
print(f"clifpy version: {clifpy.__version__}")
print(f"clifpy location: {clifpy.__file__}")

print(f"\nclifpy version: {clifpy.__version__}")

print(f"\n=== Working Directory ===")
print(f"Current directory: {os.getcwd()}")


# In[2]:


# Load configuration
config_path = "../config/config.json"
with open(config_path, 'r') as f:
    config = json.load(f)

## import outlier json
# with open('../config/outlier_config.json', 'r', encoding='utf-8') as f:
#     outlier_cfg = json.load(f)

print(f"\n=== Configuration:")
has_crrt_settings = config.get('has_crrt_settings', True)

print(f"   Data directory: {config['tables_path']}")
print(f"   File type: {config['file_type']}")
print(f"   Timezone: {config['timezone']}")
print(f"   Has CRRT settings: {has_crrt_settings}")


# In[3]:


import os
# Create output directories if they do not exist
os.makedirs("../output/final/graphs", exist_ok=True)
os.makedirs("../output/intermediate", exist_ok=True)


# # Required columns and categories

# In[4]:


print("\n" + "=" * 80)
print("Defining Required Data Elements")
print("=" * 80)

# Full patient table 

# Full hospitalization table 

# Full ADT table

# Vitals
vitals_required_columns = [
    'hospitalization_id',
    'recorded_dttm',
    'vital_category',
    'vital_value'
]
vitals_of_interest = ['heart_rate', 'respiratory_rate', 'sbp', 'dbp', 'map', 'spo2', 'weight_kg', 'height_cm']

#Labs
labs_required_columns = [
    'hospitalization_id',
    'lab_result_dttm',
    'lab_category',
    'lab_value',
    'lab_value_numeric'
]
labs_of_interest = ['po2_arterial','pco2_arterial', 'ph_arterial','ph_venous', 'bicarbonate','so2_arterial',
                    'sodium', 'potassium', 'chloride', 'calcium_total', 'magnesium', 'creatinine', 
                    'bun', 'glucose_serum', 'lactate', 'hemoglobin' ]

# Continuous administered meds
meds_required_columns = [
    'hospitalization_id',
    'admin_dttm',
    'med_name',
    'med_category',
    'med_dose',
    'med_dose_unit'
]
meds_of_interest = [
    'norepinephrine', 'epinephrine', 'phenylephrine', 'vasopressin',
    'dopamine', 'angiotensin', 'dobutamine', 'milrinone', 'isoproterenol',
    'propofol', 'midazolam', 'lorazepam', 'dexmedetomidine', 
    'vecuronium', 'rocuronium', 'cisatracurium', 'pancuronium'
]

# Respiratory Support 
rst_required_columns = [
    'hospitalization_id',
    'recorded_dttm',
    'device_name',
    'device_category',
    'mode_name', 
    'mode_category',
    'tracheostomy',
    'fio2_set',
    'lpm_set',
    'resp_rate_set',
    'peep_set',
    'resp_rate_obs',
    'tidal_volume_set',
    'pressure_control_set',
    'pressure_support_set',
    'peak_inspiratory_pressure_set',
    'peak_inspiratory_pressure_obs',
    'plateau_pressure_obs',
    'minute_vent_obs',
]

# CRRT table — columns depend on has_crrt_settings flag
if has_crrt_settings:
    crrt_required_columns = [
        'hospitalization_id',
        'recorded_dttm',
        'crrt_mode_category',
        'blood_flow_rate',
        'pre_filter_replacement_fluid_rate',
        'post_filter_replacement_fluid_rate',
        'dialysate_flow_rate',
        'ultrafiltration_out'
    ]
else:
    crrt_required_columns = [
        'hospitalization_id',
        'recorded_dttm'
    ]


# # Cohort Identification
# 
# 
# **Inclusion**
# 1. Adults
# 2. Admitted between January 1, 2018 to December, 31, 2024
# 3. Receiving CRRT- must have DFR or UF documented at any point in the hospitalization
# 4. Data completeness- Must have weight & CRRT settings  documented
# 
# **Exclusion**
# 1. Prior to admission ICD codes for ESRD

# In[5]:


strobe_counts = {}


# ## Step0: Load Core Tables

# In[6]:


print("\n" + "=" * 80)
print("Loading CLIF Tables")
print("=" * 80)

from clifpy.clif_orchestrator import ClifOrchestrator

# Initialize ClifOrchestrator
clif = ClifOrchestrator(
    data_directory=config['tables_path'],
    filetype=config['file_type'],
    timezone=config['timezone']
)


# In[7]:


# ============================================================================
# STEP 0: Load Core Tables (Patient, Hospitalization, ADT)
# ============================================================================
print("\n" + "=" * 80)
print("Step 0: Load Core Tables (Patient, Hospitalization, ADT)")
print("=" * 80)
core_tables = ['patient', 'hospitalization', 'adt']

print(f"\nLoading {len(core_tables)} core tables...")
for table_name in core_tables:
    print(f"   Loading {table_name}...", end=" ")
    try:
        clif.load_table(table_name)
        table = getattr(clif, table_name)
        print(f"✓ ({len(table.df):,} rows)")
    except Exception as e:
        print(f"✗ Error: {e}")
        raise

print("\nCore tables loaded successfully!")


# In[8]:


hosp_df = clif.hospitalization.df
adt_df = clif.adt.df

# Merge to get age information
all_encounters = pd.merge(
    hosp_df[["patient_id", "hospitalization_id", "admission_dttm", "discharge_dttm", 
             "age_at_admission", "discharge_category"]],
    adt_df[["hospitalization_id", "hospital_id", "in_dttm", "out_dttm", 
            "location_category", "location_type"]],
    on='hospitalization_id',
    how='inner'
)


# In[9]:


# Check for duplicates by ['hospitalization_id', 'in_dttm', 'out_dttm']
dup_counts = all_encounters.duplicated(subset=['hospitalization_id', 'in_dttm', 'out_dttm']).sum()
if dup_counts > 0:
    print(f"Warning: {dup_counts} duplicate (hospitalization_id, in_dttm, out_dttm) entries found in all_encounters.")
else:
    print("No duplicate (hospitalization_id, in_dttm, out_dttm) entries found in all_encounters.")


# ## Step1: Date & Age filter

# In[10]:


# ============================================================================
# STEP 1: Identify Adult Patients (Age >= 18) and Admissions 2018-2024
# ============================================================================
print("\n" + "=" * 80)
print("Step 1: Identifying Adult Patients (Age >= 18) and Admissions 2018-2024")
print("=" * 80)

print("Applying initial cohort filters...")

# Use only the relevant columns from all_encounters
adult_encounters = all_encounters[
    [
        'patient_id', 'hospitalization_id', 'admission_dttm', 'discharge_dttm',
        'age_at_admission', 'discharge_category', 'hospital_id',
        'in_dttm', 'out_dttm', 'location_category', 'location_type'
    ]
].copy()

# Filter for adult patients (age >= 18) and valid age
adult_encounters = adult_encounters[
    (adult_encounters['age_at_admission'] >= 18) & (adult_encounters['age_at_admission'].notna())
]

# Filter for admission years 2018-2024
adult_encounters = adult_encounters[
    (adult_encounters['admission_dttm'].dt.year >= 2018) & (adult_encounters['admission_dttm'].dt.year <= 2024)
]

print(f"\nFiltering Results:")
print(f"   Total hospitalizations: {len(all_encounters['hospitalization_id'].unique()):,}")
print(f"   Adult hospitalizations (age >= 18, 2018-2024): {len(adult_encounters['hospitalization_id'].unique()):,}")
print(f"   Excluded (age < 18 or outside 2018-2024): {len(all_encounters['hospitalization_id'].unique()) - len(adult_encounters['hospitalization_id'].unique()):,}")


strobe_counts["0_total_hospitalizations"] = len(all_encounters['hospitalization_id'].unique())
strobe_counts["1_adult_hospitalizations"] = len(adult_encounters['hospitalization_id'].unique())
# Get list of adult hospitalization IDs for filtering
adult_hosp_ids = set(adult_encounters['hospitalization_id'].unique())
print(f"\n   Unique adult hospitalization IDs: {len(adult_hosp_ids):,}")


# ## Step1B: Stitch hospitalizations

# In[11]:


from clifpy.utils.stitching_encounters import stitch_encounters

# Instead of multiple copies, work with references and clean up
hosp_filtered = clif.hospitalization.df[clif.hospitalization.df['hospitalization_id'].isin(adult_hosp_ids)]
adt_filtered = clif.adt.df[clif.adt.df['hospitalization_id'].isin(adult_hosp_ids)]

hosp_stitched, adt_stitched, encounter_mapping = stitch_encounters(
    hospitalization=hosp_filtered,
    adt=adt_filtered,
    time_interval=6  
)

# Direct assignment without additional copies
clif.hospitalization.df = hosp_stitched
clif.adt.df = adt_stitched

# Store the encounter mapping in the orchestrator for later use
clif.encounter_mapping = encounter_mapping

# Clean up intermediate variables
del hosp_filtered, adt_filtered
gc.collect()


# In[12]:


# After your stitching code, add these calculations:

# Calculate stitching statistics
strobe_counts['1b_before_stitching'] = len(adult_hosp_ids)  # Original adult hospitalizations
strobe_counts['1b_after_stitching'] = len(hosp_stitched['encounter_block'].unique())  # Unique encounter blocks after stitching
strobe_counts['1b_stitched_hosp_ids'] = strobe_counts['1b_before_stitching'] - strobe_counts['1b_after_stitching']  # Number of hospitalizations that were linked

print(f"\nEncounter Stitching Results:")
print(f"   Number of unique hospitalizations before stitching: {strobe_counts['1b_before_stitching']:,}")
print(f"   Number of unique encounter blocks after stitching: {strobe_counts['1b_after_stitching']:,}")
print(f"   Number of linked hospitalization ids: {strobe_counts['1b_stitched_hosp_ids']:,}")

# Optional: Show the encounter mapping details
print(f"\nEncounter Mapping Details:")
print(f"   Total encounter mappings created: {len(encounter_mapping):,}")
if len(encounter_mapping) > 0:
    # Show some examples of how many original hospitalizations were combined
    mapping_counts = encounter_mapping.groupby('encounter_block').size()
    print(f"   Encounter blocks with multiple hospitalizations: {(mapping_counts > 1).sum():,}")
    print(f"   Maximum hospitalizations combined into one block: {mapping_counts.max()}")


# In[13]:


cohort_df = encounter_mapping.copy()


# ## Step2: Identify CRRT Encounters

# In[14]:


print(f"\nLoading crrt_therapy table...")
try:
    clif.load_table(
        'crrt_therapy',
        filters={'hospitalization_id': list(adult_hosp_ids)}
    )
    print(f"   CRRT therapy loaded: {len(clif.crrt_therapy.df):,} rows")
    print(f"   Unique CRRT therapy hospitalizations: {clif.crrt_therapy.df['hospitalization_id'].nunique()}")
except Exception as e:
    print(f"   CRRT therapy not available or error: {e}")
    raise


# In[15]:


# Update CRRT therapy DataFrame with encounter blocks
clif.crrt_therapy.df = clif.crrt_therapy.df.merge(
    clif.encounter_mapping[['hospitalization_id', 'encounter_block']],
    on='hospitalization_id',
    how='left'
)

n_crrt_hosp = clif.crrt_therapy.df['hospitalization_id'].nunique()
n_crrt_blocks = clif.crrt_therapy.df['encounter_block'].nunique()
crrt_hosp_ids = set(clif.crrt_therapy.df['hospitalization_id'].unique())

print(f"Updated CRRT therapy DataFrame:")
print(f"   Total CRRT records: {len(clif.crrt_therapy.df):,}")
print(f"   Records with encounter blocks: {clif.crrt_therapy.df['encounter_block'].notna().sum():,}")
print(f"   Unique encounter blocks in CRRT data: {n_crrt_blocks}")
print(f"   Unique hospitalizations  in CRRT data: {n_crrt_hosp}")

strobe_counts["2_crrt_hospitalizations"] = n_crrt_hosp
strobe_counts["2_crrt_blocks"] = n_crrt_blocks

# Filter cohort_df to only hospitalizations present in CRRT data
cohort_df = cohort_df[cohort_df['hospitalization_id'].isin(crrt_hosp_ids)].copy()
crrt_df = clif.crrt_therapy.df

# --- Filter to valid CRRT modes (exclude IHD and other non-continuous modalities) ---
VALID_CRRT_MODES = {'scuf', 'cvvh', 'cvvhd', 'cvvhdf', 'avvh'}
crrt_df['crrt_mode_category'] = crrt_df['crrt_mode_category'].str.lower().str.strip()

invalid_mask = ~crrt_df['crrt_mode_category'].isin(VALID_CRRT_MODES) & crrt_df['crrt_mode_category'].notna()
n_invalid = invalid_mask.sum()
if n_invalid > 0:
    invalid_modes = crrt_df.loc[invalid_mask, 'crrt_mode_category'].value_counts()
    print(f"   Excluded {n_invalid:,} rows with non-CRRT modes:")
    print(invalid_modes.to_string())
    strobe_counts['excluded_non_crrt_modes'] = int(n_invalid)
    crrt_df = crrt_df[~invalid_mask].copy()

# Recount after filtering
n_crrt_blocks = crrt_df['encounter_block'].nunique()
strobe_counts["2_crrt_blocks"] = n_crrt_blocks
print(f"   CRRT encounter blocks after mode filter: {n_crrt_blocks:,}")


# In[16]:


if has_crrt_settings:
    from utils import handle_crrt_outliers

    # Apply outlier removal
    crrt_df, outlier_summary = handle_crrt_outliers(
        crrt_df,
        config_path='../config/outlier_config.json'
    )

if has_crrt_settings:
    # ============================================================================
    # Validate Outlier Ranges
    # ============================================================================
    import matplotlib.pyplot as plt
    import seaborn as sns

    print("\n" + "=" * 80)
    print("CRRT Parameter Distribution Validation")
    print("=" * 80)

    # CRRT parameters to validate and plot
    crrt_params = {
        'blood_flow_rate': {'expected_range': [150, 5000], 'unit': 'mL/min'},
        'dialysate_flow_rate': {'expected_range': [0, 20000], 'unit': 'mL/hr'},
        'pre_filter_replacement_fluid_rate': {'expected_range': [0, 20000], 'unit': 'mL/hr'},
        'post_filter_replacement_fluid_rate': {'expected_range': [0, 20000], 'unit': 'mL/hr'},
        'ultrafiltration_out': {'expected_range': [0, 20000], 'unit': 'mL/hr'}
    }

    # Check if means are within expected ranges
    print("\n1. Validating parameter distributions...")
    unit_warnings = []

    for param, info in crrt_params.items():
        if param not in crrt_df.columns:
            continue

        values = crrt_df[param].dropna()
        if len(values) == 0:
            continue

        mean_val = values.mean()
        min_range, max_range = info['expected_range']
        unit = info['unit']

        print(f"\n   {param}:")
        print(f"     Mean: {mean_val:.0f} {unit}")
        print(f"     Expected range: {min_range}-{max_range} {unit}")

        # Check if mean is within range
        if mean_val < min_range or mean_val > max_range:
            warning_msg = f"⚠️  WARNING: Mean ({mean_val:.0f}) is OUTSIDE expected range [{min_range}-{max_range}]"
            print(f"     {warning_msg}")
            print(f"     → Check if units are correct (expected: {unit})")
            unit_warnings.append({
                'parameter': param,
                'mean': mean_val,
                'expected_range': [min_range, max_range],
                'unit': unit
            })
        else:
            print(f"     ✓ Mean is within expected range")

    if unit_warnings:
        print("\n" + "!" * 80)
        print("UNIT MISMATCH WARNINGS:")
        for w in unit_warnings:
            print(f"\n   {w['parameter']}:")
            print(f"     Mean: {w['mean']:.0f}")
            print(f"     Expected range: {w['expected_range']} {w['unit']}")
            print(f"     → Data may be in different units than expected!")
        print("!" * 80)
    else:
        print("\n✓ All parameter means are within expected ranges")

    print("=" * 80)

    # ============================================================================
    # Generate GRID Histograms: Parameters (rows) × Modes (columns)
    # ============================================================================
    print("\n" + "=" * 80)
    print("Generating CRRT Parameter Histograms Grid by Mode")
    print("=" * 80)

    # Get unique modes (sorted)
    modes = sorted(crrt_df['crrt_mode_category'].dropna().unique())
    n_modes = len(modes)
    n_params = len(crrt_params)

    # Create grid: rows = parameters, columns = modes
    fig, axes = plt.subplots(n_params, n_modes, figsize=(5 * n_modes, 4 * n_params))

    # Handle single row/column cases
    if n_params == 1 and n_modes == 1:
        axes = np.array([[axes]])
    elif n_params == 1:
        axes = axes.reshape(1, -1)
    elif n_modes == 1:
        axes = axes.reshape(-1, 1)

    fig.suptitle('CRRT Parameter Distributions by Mode (After Outlier Removal)',
                fontsize=16, fontweight='bold', y=0.995)

    # Plot grid
    for row_idx, (param, info) in enumerate(crrt_params.items()):
        for col_idx, mode in enumerate(modes):
            ax = axes[row_idx, col_idx]

            # Check if parameter exists in data
            if param not in crrt_df.columns:
                ax.text(0.5, 0.5, f'{param}\nNot Available',
                        ha='center', va='center', fontsize=10)
                ax.set_title(f'{mode.upper()}')
                continue

            # Filter data for this mode and parameter
            mode_data = crrt_df[
                (crrt_df['crrt_mode_category'] == mode) &
                (crrt_df[param].notna())
            ][param]

            if len(mode_data) == 0:
                ax.text(0.5, 0.5, 'No Data', ha='center', va='center', fontsize=10)
                ax.set_title(f'{mode.upper()}')
                continue

            # Create histogram
            ax.hist(mode_data, bins=20, alpha=0.7, color='steelblue', edgecolor='black')

            # Calculate statistics
            n = len(mode_data)
            mean_val = mode_data.mean()
            median_val = mode_data.median()

            # Add vertical lines for mean and median
            ax.axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_val:.0f}')
            ax.axvline(median_val, color='green', linestyle=':', linewidth=2, label=f'Median: {median_val:.0f}')

            # Add statistics text box
            stats_text = f'N = {n:,}\nMean = {mean_val:.0f}\nMedian = {median_val:.0f}'
            ax.text(0.98, 0.97, stats_text, transform=ax.transAxes,
                    verticalalignment='top', horizontalalignment='right',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                    fontsize=9)

            # Set title (mode name) for top row only
            if row_idx == 0:
                ax.set_title(f'{mode.upper()}', fontsize=12, fontweight='bold')

            # Set ylabel (parameter name) for first column only
            if col_idx == 0:
                ax.set_ylabel(f'{param.replace("_", " ").title()}\n({info["unit"]})',
                            fontsize=10, fontweight='bold')

            # Set xlabel for bottom row only
            if row_idx == n_params - 1:
                ax.set_xlabel(f'{info["unit"]}', fontsize=9)

            ax.grid(True, alpha=0.3)
            ax.legend(loc='upper left', fontsize=8)

    plt.tight_layout()

    # Save figure
    Path("../output/final").mkdir(parents=True, exist_ok=True)
    plt.savefig('../output/final/graphs/crrt_parameter_histograms_grid.png', dpi=300, bbox_inches='tight')
    plt.close()

    print("\n✓ Grid histograms saved to: output/final/crrt_parameter_histograms_grid.png")
    print("=" * 80)

    # ============================================================================
    # Summary Statistics: Distribution of Settings by Mode Category
    # ============================================================================
    print("\n" + "=" * 80)
    print("CRRT Settings Distribution by Mode Category")
    print("=" * 80)

    summary_data = []

    for mode in modes:
        mode_df = crrt_df[crrt_df['crrt_mode_category'] == mode]

        row = {
            'Mode': mode.upper(),
            'N_Total': len(mode_df)
        }

        for param, info in crrt_params.items():
            if param in mode_df.columns:
                values = mode_df[param].dropna()
                n = len(values)

                if n > 0:
                    row[f'{param}_N'] = n
                    row[f'{param}_Mean'] = round(values.mean(), 1)
                    row[f'{param}_SD'] = round(values.std(), 1)
                    row[f'{param}_Median'] = round(values.median(), 1)
                    row[f'{param}_Q25'] = round(values.quantile(0.25), 1)
                    row[f'{param}_Q75'] = round(values.quantile(0.75), 1)
                    row[f'{param}_Min'] = round(values.min(), 1)
                    row[f'{param}_Max'] = round(values.max(), 1)
                else:
                    row[f'{param}_N'] = 0
                    row[f'{param}_Mean'] = np.nan
                    row[f'{param}_SD'] = np.nan
                    row[f'{param}_Median'] = np.nan
                    row[f'{param}_Q25'] = np.nan
                    row[f'{param}_Q75'] = np.nan
                    row[f'{param}_Min'] = np.nan
                    row[f'{param}_Max'] = np.nan

        summary_data.append(row)

    summary_df = pd.DataFrame(summary_data)

    # Display summary
    print("\nSummary by Mode (first few columns):")
    display_cols = ['Mode', 'N_Total'] + [col for col in summary_df.columns if '_Mean' in col or '_Median' in col]
    print(summary_df[display_cols].to_string(index=False))

    # Save detailed summary
    summary_df.to_csv('../output/final/crrt_settings_distribution_by_mode.csv', index=False)
    print(f"\n✓ Detailed summary saved to: output/final/crrt_settings_distribution_by_mode.csv")

    # Also create a simplified summary for quick reference
    simple_summary = []
    for mode in modes:
        mode_df = crrt_df[crrt_df['crrt_mode_category'] == mode]

        row = {'Mode': mode.upper(), 'N': len(mode_df)}

        for param, info in crrt_params.items():
            if param in mode_df.columns:
                values = mode_df[param].dropna()
                if len(values) > 0:
                    # Format as "Median [Q25-Q75]"
                    row[param] = f"{values.median():.0f} [{values.quantile(0.25):.0f}-{values.quantile(0.75):.0f}]"
                else:
                    row[param] = "No data"

        simple_summary.append(row)

    simple_summary_df = pd.DataFrame(simple_summary)

    # Display simple summary
    print("\n" + "=" * 80)
    print("Quick Reference: Median [IQR] by Mode")
    print("=" * 80)
    print(simple_summary_df.to_string(index=False))

    # Save simple summary
    simple_summary_df.to_csv('../output/final/crrt_settings_summary_simple.csv', index=False)
    print(f"\n✓ Simple summary saved to: output/final/crrt_settings_summary_simple.csv")
    print("=" * 80)


else:
    print("CRRT settings not available — skipping outlier validation and parameter histograms")

# In[17]:


if has_crrt_settings:
    # Drop rows where crrt_mode_category == "scuf"
    n_rows_before = len(crrt_df)
    n_encounter_blocks_before = crrt_df['encounter_block'].nunique()

    # Identify encounter_blocks where all rows are scuf (i.e., only scuf)
    only_scuf_blocks = (
        crrt_df.groupby('encounter_block')['crrt_mode_category']
        .apply(lambda x: set(x.dropna()) == set(['scuf']))
        .loc[lambda x: x].index
    )
    n_only_scuf_blocks = len(only_scuf_blocks)

    # Drop rows where crrt_mode_category == "scuf"
    crrt_df = crrt_df[crrt_df['crrt_mode_category'].str.lower() != "scuf"]
    n_rows_after = len(crrt_df)
    n_encounter_blocks_after = crrt_df['encounter_block'].nunique()
    n_rows_dropped = n_rows_before - n_rows_after

    print(f"Dropped {n_rows_dropped} rows where crrt_mode_category == 'scuf'")
    print(f"Encounter_blocks with only SCUF (all records scuf): {n_only_scuf_blocks}")
    strobe_counts["blocks_with_only_scuf"] = n_only_scuf_blocks


# In[18]:


print("\n" + "=" * 80)
print("Processing CRRT Data")
print("=" * 80)
print(f"   CRRT therapy loaded: {len(clif.crrt_therapy.df):,} rows")
print(f"   Unique CRRT therapy hospitalizations: {clif.crrt_therapy.df['hospitalization_id'].nunique()}")
# ============================================================================
#  Define CRRT Initiation Time
# ============================================================================
print("\n1. Defining CRRT initiation time...")

# Filter crrt_df to only include hospitalization_ids present in the cohort
crrt_cohort = crrt_df[crrt_df['hospitalization_id'].isin(cohort_df['hospitalization_id'])].copy()

# Sort by encounter_block and time
crrt_cohort = crrt_cohort.sort_values(['encounter_block', 'recorded_dttm'])

if has_crrt_settings:
    # Create indicator for any CRRT activity (any non-null flow rate)
    crrt_cohort['has_crrt_activity'] = (
        (
            crrt_cohort['dialysate_flow_rate'].notna() & (crrt_cohort['dialysate_flow_rate'] > 0)
        ) |
        (
            crrt_cohort['pre_filter_replacement_fluid_rate'].notna() & (crrt_cohort['pre_filter_replacement_fluid_rate'] > 0)
        ) |
        (
            crrt_cohort['post_filter_replacement_fluid_rate'].notna() & (crrt_cohort['post_filter_replacement_fluid_rate'] > 0)
        )
    )

    # Get CRRT initiation time (first non-null flow rate per encounter_block)
    crrt_initiation = (crrt_cohort[crrt_cohort['has_crrt_activity']]
                        .groupby('encounter_block')
                        .agg({'recorded_dttm': 'min'})
                        .reset_index())
    crrt_initiation.rename(columns={'recorded_dttm': 'crrt_initiation_time'}, inplace=True)
else:
    # Without CRRT settings, use the first crrt_therapy record as initiation
    crrt_initiation = (crrt_cohort
                        .groupby('encounter_block')
                        .agg({'recorded_dttm': 'min'})
                        .reset_index())
    crrt_initiation.rename(columns={'recorded_dttm': 'crrt_initiation_time'}, inplace=True)

print(f"   CRRT initiation times identified for: {len(crrt_initiation):,} encounter_blocks")
print(f"   Date range: {crrt_initiation['crrt_initiation_time'].min()} to {crrt_initiation['crrt_initiation_time'].max()}")


# In[19]:


# Show number of unique encounter_blocks in crrt_cohort and crrt_initiation
num_blocks_cohort = crrt_cohort['encounter_block'].nunique()
num_blocks_initiation = crrt_initiation['encounter_block'].nunique()

print(f"Unique encounter_blocks in crrt_cohort: {num_blocks_cohort}")
print(f"Unique encounter_blocks in crrt_initiation: {num_blocks_initiation}")

# Analyze and explain the apparent paradox
blocks_cohort_set = set(crrt_cohort['encounter_block'])
blocks_initiation_set = set(crrt_initiation['encounter_block'])

missing_initiation = blocks_cohort_set - blocks_initiation_set
missing_cohort = blocks_initiation_set - blocks_cohort_set

print(f"Encounter_blocks in crrt_cohort not in crrt_initiation: {len(missing_initiation)}")
print(f"Encounter_blocks in crrt_initiation not in crrt_cohort: {len(missing_cohort)}")

# Filter final_crrt to only include encounter_blocks present in crrt_initiation['encounter_block']
final_crrt = clif.crrt_therapy.df[
    clif.crrt_therapy.df['encounter_block'].isin(crrt_initiation['encounter_block'])
].copy()


# In[20]:


# ============================================================================
# Get CRRT Parameters at Initiation
# ============================================================================
print("\n4. Getting CRRT parameters at initiation time...")
# Merge initiation times back to crrt_cohort
crrt_cohort = crrt_cohort.merge(crrt_initiation, on='encounter_block', how='left')
# Filter to records at exactly the initiation time
crrt_at_initiation = crrt_cohort[
    crrt_cohort['recorded_dttm'] == crrt_cohort['crrt_initiation_time']
].copy()

print(f"   CRRT records at initiation: {len(crrt_at_initiation):,}")
print(f"   Unique encounter blocks: {crrt_at_initiation['encounter_block'].nunique():,}")



# ## Step3: Exclude ESRD encounters
# 
# Prior to admission ICD codes for ESRD

# In[21]:


print(f"\nLoading Hospital dx table...")
try:
    clif.load_table(
        'hospital_diagnosis',
        filters={'hospitalization_id': list(crrt_hosp_ids)}
    )
    print(f"   Hospital dx loaded: {len(clif.hospital_diagnosis.df):,} rows")
    print(f"   Unique Hospital dx hospitalizations: {clif.hospital_diagnosis.df['hospitalization_id'].nunique()}")

    print("Merge encounter blocks with diagnosis")
    clif.hospital_diagnosis.df = clif.hospital_diagnosis.df.merge(
                    clif.encounter_mapping[['hospitalization_id', 'encounter_block']],
                    on='hospitalization_id',
                    how='left')

    n_dx_hosp = clif.hospital_diagnosis.df['hospitalization_id'].nunique()
    n_dx_blocks = clif.hospital_diagnosis.df['encounter_block'].nunique()
    cohort_hosp_ids = set(clif.hospital_diagnosis.df['hospitalization_id'].unique())
    cohort_blocks = set(clif.hospital_diagnosis.df['encounter_block'].unique())
    print(f"   Total Hospital dx records: {len(clif.hospital_diagnosis.df):,}")
    print(f"   Records with encounter blocks: {clif.hospital_diagnosis.df['encounter_block'].notna().sum():,}")
    print(f"   Unique encounter blocks in Hospital dx data: {n_dx_blocks}")
    print(f"   Unique hospitalizations  in Hospital dx data: {n_dx_hosp}")
except Exception as e:
    print(f"   Hospital dx not available or error: {e}")
    raise


# In[22]:


hospital_diagnosis_df = clif.hospital_diagnosis.df.copy()

print("Hospital dx column names :", hospital_diagnosis_df.columns)
# Clean and standardize diagnosis codes
hospital_diagnosis_df['diagnosis_code'] = hospital_diagnosis_df['diagnosis_code'].str.replace('.', '').str.lower()

if 'present_on_admission' in hospital_diagnosis_df.columns:
    hospital_diagnosis_df = hospital_diagnosis_df.rename(columns={'present_on_admission': 'poa_present'})

# Check present_on_admission column type and standardize to int8
if 'poa_present' in hospital_diagnosis_df.columns:
    # Only allow 1 (present on admission) or 0 (not present on admission)
    # Any other value (including Exempt, Unknown, Unspecified, NA) is set to 0
    hospital_diagnosis_df['poa_present'] = hospital_diagnosis_df['poa_present'].astype(str).str.lower()
    hospital_diagnosis_df['poa_present'] = hospital_diagnosis_df['poa_present'].map(
        {'yes': 1, 'y': 1, 'true': 1, '1': 1, 'no': 0, 'n': 0, 'false': 0, '0': 0}
    ).fillna(0).astype('int8')


# In[23]:


# Define ESRD diagnosis codes
# Let's debug why we're not finding ESRD codes
esrd_codes = [
    'z992',    # Dependence on renal dialysis
    'z9115',   # Patient's noncompliance with renal dialysis
    'i120',    # Hypertensive chronic kidney disease with stage 5 CKD or ESRD
    'n186',    # End stage renal disease
    'i132',    # Hypertensive heart and chronic kidney disease with heart failure and ESRD
    'z91158',  # Patient's noncompliance with renal dialysis (alternate code)
    'i1311',   # Hypertensive heart and chronic kidney disease with heart failure and stage 5 CKD
    "z4931",   # Encounter for continuous renal replacement therapy (CRRT) for ESRD (CMS/HCC)
    "z4901",   # Encounter regarding vascular access for dialysis for end-stage renal disease (CMS/HCC)
    "i272",    # Pulmonary hypertension associated with ESRD on dialysis (CMS/HCC)
    '5856',     #ICD9 :End stage renal disease
    '40391',    #ICD9: Hypertensive chronic kidney disease, unspecified, with chronic kidney disease stage V or end stage renal disease
    '40311',     #ICD9: Hypertensive chronic kidney disease, benign, with chronic kidney disease stage V or end stage renal disease
    'v4511',     #ICD9: Renal dialysis status
    'v4512'     #ICD9: Noncompliance with renal dialysis
]

# Get hospitalization IDs with ESRD diagnoses and print debug info
print("\nNumber of rows matching ESRD codes:", hospital_diagnosis_df['diagnosis_code'].isin(esrd_codes).sum())


# Count how many ESRD codes have present_on_admission = 1, 0, or NA
esrd_poa_counts = hospital_diagnosis_df[
    hospital_diagnosis_df['diagnosis_code'].isin(esrd_codes)
]['poa_present'].value_counts(dropna=False)
print("Present_on_admission values for ESRD codes:")
print(esrd_poa_counts)

# Use a more inclusive approach for ESRD identification
# Include cases where present_on_admission is 1 OR NA (assuming NA means unknown/possible)
esrd_mask = (
    hospital_diagnosis_df['diagnosis_code'].isin(esrd_codes) & 
    ((hospital_diagnosis_df['poa_present'] == 1) | 
        (hospital_diagnosis_df['poa_present'].isna()))
)
hosp_ids_with_esrd = hospital_diagnosis_df[esrd_mask]['hospitalization_id'].unique()
blocks_with_esrd = hospital_diagnosis_df[esrd_mask]['encounter_block'].unique()

print(f"Hospitalizations with ESRD (including NA present_on_admission): {len(hosp_ids_with_esrd)}")


strobe_counts['3_hospitalizations_with_esrd'] = len(hosp_ids_with_esrd)
strobe_counts['3_encounter_blocks_with_esrd'] = len(blocks_with_esrd)


# Filter out hospitalizations with ESRD
cohort_df = cohort_df[~cohort_df['hospitalization_id'].isin(hosp_ids_with_esrd)].copy()
cohort_hosp_ids = set(cohort_df['hospitalization_id'].unique())
cohort_blocks = set(cohort_df['encounter_block'].unique())
# Create cohort subset excluding hospitalizations with ESRD
strobe_counts['3_encounter_blocks_without_esrd'] = len(cohort_blocks)  # Count blocks without ESRD
strobe_counts['3_hospitalizations_without_esrd'] = len(cohort_hosp_ids)  # Count hospitalizations without ESRD

strobe_counts


# ## Step4: Data availability, and CRRT Settings

# In[24]:


print(f"\nLoading labs table...")
clif.load_table(
    'vitals',
    columns=vitals_required_columns,
    filters={
        'hospitalization_id': list(cohort_hosp_ids)
    }
)
print(f"   Vitals loaded: {len(clif.vitals.df):,} rows")
print(f"   Unique vitals categories: {clif.vitals.df['vital_category'].nunique()}")
print(f"   Unique vitals hospitalizations: {clif.vitals.df['hospitalization_id'].nunique()}")

clif.vitals.df = clif.vitals.df.merge(
    clif.encounter_mapping[['hospitalization_id', 'encounter_block']],
    on='hospitalization_id',
    how='left'
)


# In[25]:


vitals_range = clif.vitals.df.groupby('encounter_block').agg({
    'recorded_dttm': ['min', 'max']
}).reset_index()
vitals_range.columns = ['encounter_block', 'first_vital_dttm', 'last_vital_dttm']


# In[26]:


# Keep only rows where vital_category is 'weight_kg'
weight_df = clif.vitals.df[clif.vitals.df['vital_category'] == 'weight_kg'].copy()
# Identify the number of hospitalizations that do not have weight recorded
hosp_with_weight = set(weight_df['hospitalization_id'].unique())
hosp_without_weight = cohort_hosp_ids - hosp_with_weight
print(f"Number of hospitalizations without recorded weight: {len(hosp_without_weight)}")

clif.vitals.df = None ## clear from memory


# In[27]:


cohort_df = cohort_df[~cohort_df['hospitalization_id'].isin(hosp_without_weight)].copy()
cohort_hosp_ids = set(cohort_df['hospitalization_id'].unique())
cohort_blocks = set(cohort_df['encounter_block'].unique())


# In[28]:


crrt_initiation = crrt_initiation[crrt_initiation['encounter_block'].isin(cohort_blocks)].copy()
crrt_at_initiation = crrt_at_initiation[crrt_at_initiation['encounter_block'].isin(cohort_blocks)].copy()
print(f"   CRRT initiation blocks in cohort: {len(crrt_initiation):,}")
strobe_counts['4_encounter_blocks_with_weight'] = crrt_initiation['encounter_block'].nunique()
strobe_counts['4_encounter_blocks_without_weight'] = len(cohort_blocks - set(crrt_initiation['encounter_block'].unique()))
strobe_counts['4_hospitalizations_with_weight'] = crrt_at_initiation['hospitalization_id'].nunique()  # Count hospitalizations without weight


# In[29]:


strobe_counts


# ## Add weights to crrt initiation

# In[30]:


# ============================================================================
#  Get Closest Weight to CRRT Initiation
# ============================================================================
print("\nFinding closest weights to CRRT initiation time...")
if 'vital_value' in weight_df.columns:
    weight_df = weight_df.rename(columns={'vital_value': 'weight_kg'})
if 'vital_category' in weight_df.columns:
    weight_df = weight_df.drop(columns=['vital_category'])
print(f"   Weight records available: {len(weight_df):,}")

combined = crrt_initiation.merge(weight_df, on='encounter_block', how='inner')

before_mask = combined['recorded_dttm'] <= combined['crrt_initiation_time']
combined_before = combined[before_mask].copy()
combined_before_sorted = combined_before.sort_values(['encounter_block', 'recorded_dttm'])
closest_before = (combined_before_sorted
                  .groupby('encounter_block')
                  .last()
                  .reset_index())

all_blocks = set(combined['encounter_block'])
blocks_with_before = set(closest_before['encounter_block'])
blocks_missing = all_blocks - blocks_with_before
print("Blocks without weight recorded before initiation time:", len(blocks_missing))


# In[31]:


after_mask = (
    combined['encounter_block'].isin(blocks_missing) &
    (combined['recorded_dttm'] > combined['crrt_initiation_time']) &
    (combined['weight_kg'].notnull())
)
combined_after = combined[after_mask].copy()
combined_after_sorted = combined_after.sort_values(['encounter_block', 'recorded_dttm'])
first_after = (combined_after_sorted
               .groupby('encounter_block')
               .first()
               .reset_index())

num_taken_after = len(first_after)
print(f"   Number of weights from after initiation: {num_taken_after}")


# In[32]:


combined_final = pd.concat([closest_before, first_after], axis=0, ignore_index=True)
combined = combined_final

closest_weights = (combined
                .sort_values(['encounter_block', 'recorded_dttm'])
                .groupby('encounter_block')
                .last()
                .reset_index())

closest_weights = closest_weights[['encounter_block', 'weight_kg']]

print(f"   Weights found for: {len(closest_weights):,} encounter_blocks")

# Calculate and print the number of encounter blocks for which weights were not found
all_encounter_blocks_with_crrt = set(crrt_initiation['encounter_block'].unique())
encounter_blocks_with_weights = set(closest_weights['encounter_block'].unique())
encounter_blocks_without_weights = all_encounter_blocks_with_crrt - encounter_blocks_with_weights
print(f"   Weights NOT found for: {len(encounter_blocks_without_weights):,} encounter_blocks")


# In[33]:


# ============================================================================
# STEP 6: Combine CRRT Data with Weights
# ============================================================================
print("\n6. Combining CRRT data with weights...")

index_crrt_df = crrt_at_initiation.merge(
    closest_weights,
    on='encounter_block',
    how='inner'
)

print(f"   Final dataset: {len(index_crrt_df):,} records")
print(f"   Records with weights: {index_crrt_df['weight_kg'].notna().sum():,}")
if has_crrt_settings:
    print(f"   Records with CRRT mode: {index_crrt_df['crrt_mode_category'].notna().sum():,}")


# In[34]:


# Filter crrt_df to only include hospitalization_id present in cohort_df
index_crrt_df = index_crrt_df[index_crrt_df['hospitalization_id'].isin(cohort_df['hospitalization_id'])]

if has_crrt_settings:
    # Identify encounters without any CRRT settings documented
    crrt_settings_cols = [
        'pre_filter_replacement_fluid_rate',
        'post_filter_replacement_fluid_rate',
        'dialysate_flow_rate',
        'ultrafiltration_out'
    ]
    # Find encounter_blocks with ANY crrt settings recorded
    crrt_settings_present = index_crrt_df.groupby('encounter_block')[crrt_settings_cols].apply(
        lambda df: df.notnull().any().any()
    )
    crrt_blocks_with_settings = set(crrt_settings_present[crrt_settings_present].index)
    crrt_blocks_without_settings = set(crrt_df['encounter_block'].unique()) - crrt_blocks_with_settings
    num_encounters_without_crrt_settings = len(crrt_blocks_without_settings)
    print(f"Number of encounter blocks without any recorded CRRT settings: {num_encounters_without_crrt_settings}")

    # Filter cohort_df to only include encounter_blocks with at least one CRRT setting recorded
    cohort_df = cohort_df[cohort_df['encounter_block'].isin(crrt_blocks_with_settings)].copy()
    cohort_hosp_ids = set(cohort_df['hospitalization_id'].unique())
    cohort_blocks = set(cohort_df['encounter_block'].unique())
    strobe_counts['5_encounter_blocks_with_crrt_settings'] = len(cohort_blocks)
    strobe_counts['5_hospitalizations_with_crrt_settings'] = len(cohort_hosp_ids)
else:
    # Without CRRT settings, skip this filter — keep all encounters
    print("CRRT settings not available — skipping CRRT settings filter")
    cohort_hosp_ids = set(cohort_df['hospitalization_id'].unique())
    cohort_blocks = set(cohort_df['encounter_block'].unique())
strobe_counts


# ## Step5: Labs availability

# In[35]:


#Labs
labs_required_columns = [
    'hospitalization_id',
    'lab_result_dttm',
    'lab_category',
    'lab_value',
    'lab_value_numeric'
]
# labs_of_interest = ['po2_arterial','pco2_arterial', 'ph_arterial','ph_venous', 'bicarbonate','so2_arterial',
#                     'sodium', 'potassium', 'chloride', 'calcium_total', 'magnesium', 'creatinine', 
#                     'bun', 'glucose_serum', 'lactate', 'hemoglobin' ]

# labs_of_interest = ['ph_arterial', 'lactate', 'bicarbonate', 'potassium']
labs_of_interest = ['lactate', 'bicarbonate', 'potassium']

print(f"\nLoading labs table...")
clif.load_table(
    'labs',
    columns=labs_required_columns,
    filters={
        'hospitalization_id': cohort_df['hospitalization_id'].unique().tolist(),
        'lab_category': labs_of_interest
    }
)
print(f"   Labs loaded: {len(clif.labs.df):,} rows")
print(f"   Unique lab categories: {clif.labs.df['lab_category'].nunique()}")
print(f"   Unique lab hospitalizations: {clif.labs.df['hospitalization_id'].nunique()}")

clif.labs.df = clif.labs.df.merge(
    clif.encounter_mapping[['hospitalization_id', 'encounter_block']],
    on='hospitalization_id',
    how='left'
)
# Get labs dataframe
labs_df = clif.labs.df.copy()


# In[36]:


# ============================================================================
# Get Most Recent Labs Within -12 hours to +3 hours of CRRT initiation
# ============================================================================
print("\n" + "=" * 80)
print("Processing Labs - Most Recent Within -12 hours to +3 hours of CRRT initiation")
print("=" * 80)


# Count unique encounter blocks with initiation time before merging
n_unique_encounter_blocks_before = labs_df['encounter_block'].nunique() if 'encounter_block' in labs_df.columns else 0
print(f"   Unique encounter blocks in labs_df before merging: {n_unique_encounter_blocks_before:,}")

# Count unique encounter blocks with initiation time available in crrt_initiation
n_unique_encounter_blocks_with_init = index_crrt_df['encounter_block'].nunique()
print(f"   Unique encounter blocks with initiation time in index_crrt_df: {n_unique_encounter_blocks_with_init:,}")

# Merge with CRRT initiation times to get the reference time
labs_with_crrt_time = labs_df.merge(
    index_crrt_df[['encounter_block', 'crrt_initiation_time']],
    on='encounter_block',
    how='inner'
)

# Count unique encounter blocks after merging
n_unique_encounter_blocks_after = labs_with_crrt_time['encounter_block'].nunique()
print(f"   Labs after merging with CRRT cohort: {len(labs_with_crrt_time):,}")
print(f"   Unique encounter blocks after merging: {n_unique_encounter_blocks_after:,}")


# In[37]:


# Filter labs to window: -12 hours to +3 hours of CRRT initiation
time_lower = labs_with_crrt_time['crrt_initiation_time'] -  pd.Timedelta(hours=24)
time_upper = labs_with_crrt_time['crrt_initiation_time'] + pd.Timedelta(hours=3)

labs_in_window = labs_with_crrt_time[
    (labs_with_crrt_time['lab_result_dttm'] >= time_lower) &
    (labs_with_crrt_time['lab_result_dttm'] <= time_upper)
].copy()

print(f"Labs within -24h to +3h window: {len(labs_in_window)} lab records")
print(f"Unique encounter_blocks with labs in window: {labs_in_window['encounter_block'].nunique()}")


# In[38]:


# Get encounter_blocks with each required lab
blocks_with_labs = {}
for lab in labs_of_interest:
    blocks = set(labs_in_window[labs_in_window['lab_category'] ==
lab]['encounter_block'].unique())
    blocks_with_labs[lab] = blocks
    print(f"Encounter blocks with {lab}: {len(blocks)}")

# Find encounter_blocks with ALL required labs (lactate not required for inclusion)
blocks_with_all_labs = blocks_with_labs['bicarbonate'].intersection(
    blocks_with_labs['potassium'])

print(f"\nEncounter blocks with ALL required labs:  {len(blocks_with_all_labs)}")
print(f"Encounter blocks missing at least one required lab:  {len(cohort_blocks) - len(blocks_with_all_labs)}")


# In[39]:


# Filter cohort to only include encounter_blocks with all required labs

cohort_blocks_before = len(cohort_blocks)
cohort_blocks_after = cohort_blocks.intersection(blocks_with_all_labs)

# Identify dropped encounter_blocks
dropped_blocks = cohort_blocks - cohort_blocks_after
print(f"Dropped encounter_blocks for missing labs: {len(dropped_blocks)} encounter blocks")

# Save hospitalization_id and encounter_block for dropped cases
if len(dropped_blocks) > 0:
    dropped_encounters_df = cohort_df[cohort_df['encounter_block'].isin(dropped_blocks)][['hospitalization_id', 'encounter_block']]
    dropped_encounters_df.to_parquet('../output/intermediate/dropped_missing_labs_blocks.parquet', index=False)
    print(f"Saved dropped hospitalization_id and encounter_block to ../output/intermediate/dropped_missing_labs_blocks.parquet")
else:
    print("No encounter_blocks dropped for missing labs.")

# Update the cohort_blocks set
cohort_blocks = cohort_blocks_after

# Update cohort_df
cohort_df = cohort_df[cohort_df['encounter_block'].isin(cohort_blocks)]

print(f"Cohort before lab filter: {cohort_blocks_before} encounter blocks")
print(f"Cohort after lab filter: {len(cohort_blocks)} encounter blocks")
print(f"Excluded for missing labs: {cohort_blocks_before -  len(cohort_blocks)} encounter blocks")

# Update STROBE counts
strobe_counts['6_encounter_blocks_with_required_labs'] = len(cohort_blocks)

crrt_at_initiation = crrt_at_initiation[crrt_at_initiation['encounter_block'].isin(cohort_blocks)]
index_crrt_df = index_crrt_df[index_crrt_df['encounter_block'].isin(cohort_blocks)]
crrt_initiation = crrt_initiation[crrt_initiation['encounter_block'].isin(cohort_blocks)]


# # Cohort Sanity Checks

# ## AKI
# 
# Majority of the cohort should have an ICD code for AKI

# In[40]:


# AKI Codes Sanity check

# Define AKI ICD-10 codes
aki_codes = [
    # ICD-10 codes for acute kidney injury
    'n170', 'n171', 'n172', 'n178', 'n179',  # Acute kidney failure codes
    'r34',   # Anuria and oliguria
    'n990', # Post-procedural kidney failure
    't795',  # Traumatic anuria
    '5845',  # ICD9 Acute kidney failure with lesion of tubular necrosis
    '5849',  # ICD9- Acute kidney failure, unspecified
    "5848"    # ICD9 - Acute kidney failure with other specified pathological lesion in kidney
]

# Filter to non-ESRD encounters first
non_esrd_encounters = hospital_diagnosis_df[hospital_diagnosis_df['encounter_block'].isin(cohort_df['encounter_block'])]

# Create mask for AKI diagnoses on the filtered data
aki_mask = non_esrd_encounters['diagnosis_code'].isin(aki_codes)

# Get encounter blocks with AKI diagnoses
blocks_with_aki = non_esrd_encounters[aki_mask]['encounter_block'].unique()
total_non_esrd_blocks = cohort_df['encounter_block'].nunique()
strobe_counts['6_encounter_blocks_with_AKI_no_esrd'] = len(blocks_with_aki) 

# Calculate percentage
aki_percentage = (len(blocks_with_aki) / total_non_esrd_blocks) * 100

print(f"\nPercentage of non-ESRD encounter blocks with AKI codes: {aki_percentage:.1f}%")
print(f"({len(blocks_with_aki)} out of {total_non_esrd_blocks} blocks)")
strobe_counts['6_Percentage_non_ESRD_encounter_blocks_with_AKI_codes'] = aki_percentage
# Show sample of AKI diagnoses
aki_diagnoses = non_esrd_encounters[aki_mask][['hospitalization_id', 'diagnosis_code','poa_present']].drop_duplicates()
print("\nSample of AKI-related diagnoses found: ")
aki_diagnoses['diagnosis_code'].value_counts()


# ## ICU
# 
# Cohort should ideally be an ICU hospitalization

# In[41]:


# Filter ADT data to only include hospitalizations in all_ids
adt_final_stitched = adt_stitched[adt_stitched['hospitalization_id'].isin(cohort_df['hospitalization_id'])].copy()
adt_final_stitched = adt_final_stitched.sort_values(by=['encounter_block', 'in_dttm'])
desired_order = ['hospitalization_id', 'encounter_block', 'hospital_id', 'in_dttm', 'out_dttm']
remaining_cols = [col for col in adt_final_stitched.columns if col not in desired_order]
adt_final_stitched = adt_final_stitched[desired_order + remaining_cols]

print("\n=== Validating ICU Administration ===")

adt_final_stitched['is_icu'] = adt_final_stitched['location_category'] == 'icu'

# Check if each hospitalization had at least one ICU stay
hosp_icu_status = adt_final_stitched.groupby('encounter_block')['is_icu'].any()
non_icu_hosps = hosp_icu_status[~hosp_icu_status].index.tolist()
strobe_counts["6_number_hosp_without_ICU_stay"] = len(non_icu_hosps)
print(f"\nNumber of CRRT hospitalizations without any ICU stay: {len(non_icu_hosps)}")
if len(non_icu_hosps) > 0:
    print("WARNING: Found CRRT hospitalizations without ICU stays")
    print("Number of hospitalization IDs without ICU stays:", len(non_icu_hosps), "check crrt_non_icu_df df")
else:
    print("All CRRT hospitalizations had at least one ICU stay")

crrt_non_icu_df = crrt_df[crrt_df['encounter_block'].isin(non_icu_hosps)]
crrt_non_icu_df = crrt_non_icu_df.sort_values(by=['hospitalization_id', 'encounter_block', 'recorded_dttm'])
desired_order = ['hospitalization_id', 'encounter_block', 'recorded_dttm', 'crrt_mode_category']
remaining_cols = [col for col in crrt_non_icu_df.columns if col not in desired_order]
crrt_non_icu_df = crrt_non_icu_df[desired_order + remaining_cols]
adt_df_non_icu_hosps = adt_stitched[adt_stitched['encounter_block'].isin(non_icu_hosps)]
adt_df_non_icu_hosps.to_csv('../output/intermediate/adt_df_non_icu_hosps.csv', index=False)


# # Strobe

# In[42]:


import pandas as pd

# Display strobe counts
print(strobe_counts)

# Save strobe counts to CSV in ../output/intermediate
strobe_counts_df = pd.DataFrame(list(strobe_counts.items()), columns=['counter', 'value'])
strobe_counts_df["site"] = config["site_name"]
strobe_counts_df.to_csv('../output/final/strobe_counts.csv', index=False)


# In[43]:


from pathlib import Path
from typing import Dict, Union
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

def create_consort_diagram_straight_flow(
    strobe_counts: Dict,
    output_dir: Union[str, Path] = "../output/final/graphs"
) -> Path:
    """
    Creates a CONSORT flow diagram with a straight vertical main flow 
    and exclusions branching off to the right side, connecting from the 
    vertical line *segment between nodes* (center of vertical arrow connecting top and bottom box).
    Exclusion arrows stop at the edge of exclusion box, not inside/overlap.
    """    
    # ------------------------------------------------------------------
    # 0. Setup and Data Derivation (Dynamic counts)
    # ------------------------------------------------------------------
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Fetch counts dynamically
    get_count = lambda key, default: strobe_counts.get(key, default)
    
    start_n = get_count("1b_after_stitching", get_count("1_adult_hospitalizations", 0))
    
    n_crrt = get_count("2_crrt_blocks", 0)
    n_no_esrd = get_count("3_encounter_blocks_without_esrd", n_crrt)
    n_with_weight = get_count("4_encounter_blocks_with_weight", n_no_esrd)
    n_with_settings = get_count("5_encounter_blocks_with_crrt_settings", n_with_weight)
    n_with_labs = get_count("6_encounter_blocks_with_required_labs", n_with_settings)

    # Build step list
    current_parent_n = start_n
    
    steps = []
    
    # 1. CRRT vs No CRRT
    steps.append({
        'name': 'CRRT', 'parent_n': current_parent_n, 'remaining_n': n_crrt,
        'remaining_label': "CRRT hospitalizations",
        'excluded_n': max(current_parent_n - n_crrt, 0), 'excluded_label': "Excluded: No CRRT"
    })
    current_parent_n = n_crrt

    # 2. ESRD exclusion
    steps.append({
        'name': 'ESRD', 'parent_n': current_parent_n, 'remaining_n': n_no_esrd,
        'remaining_label': "After ESRD exclusion",
        'excluded_n': max(current_parent_n - n_no_esrd, 0), 'excluded_label': "Excluded: ESRD diagnosis"
    })
    current_parent_n = n_no_esrd

    # 3. Missing weight
    steps.append({
        'name': 'Weight', 'parent_n': current_parent_n, 'remaining_n': n_with_weight,
        'remaining_label': "With documented weight",
        'excluded_n': max(current_parent_n - n_with_weight, 0), 'excluded_label': "Excluded: Missing weight"
    })
    current_parent_n = n_with_weight

    # 4. Missing CRRT settings (Omitted if excluded_n is 0, as requested)
    excluded_settings = max(current_parent_n - n_with_settings, 0)
    if excluded_settings > 0:
         steps.append({
            'name': 'Settings', 'parent_n': current_parent_n, 'remaining_n': n_with_settings,
            'remaining_label': "With CRRT settings",
            'excluded_n': excluded_settings, 'excluded_label': "Excluded: Missing CRRT settings"
        })
         current_parent_n = n_with_settings
    else:
        # If skipped, the parent count for labs comes from n_with_weight
        n_with_settings = n_with_weight

    # 5. Missing required labs
    steps.append({
        'name': 'Labs', 'parent_n': n_with_settings, 'remaining_n': n_with_labs,
        'remaining_label': "With required labs",
        'excluded_n': max(n_with_settings - n_with_labs, 0), 'excluded_label': "Excluded: Missing required labs"
    })
    
    # Filter out steps where the remaining count is zero, unless it's the starting count
    steps = [step for step in steps if step['remaining_n'] > 0 or step['parent_n'] > 0]


    # ------------------------------------------------------------------
    # 1. Figure setup and Geometry
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Geometry (X coordinates are the critical part)
    box_h = 0.08
    box_w = 0.40
    
    # Main flow (left column) X-position
    x_main_start = 0.05
    x_main_center = x_main_start + box_w / 2 # 0.25
    
    # Exclusion flow (right column) X-position
    x_excl_start = 0.55
    x_excl_center = x_excl_start + box_w / 2 # 0.75

    # Where to end the exclusion arrows (leave a gap before the exclusion box)
    excl_arrow_gap = 0.015  # in axes units, tweak for clarity
    
    v_spacing = 0.14
    
    # Helper to draw rounded box with centered text
    def draw_box(x, y, w, h, text, fontsize=11, weight="normal"):
        rect = FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.01",
            linewidth=2,
            edgecolor="black",
            facecolor="white"
        )
        ax.add_patch(rect)
        ax.text(
            x + w/2, y + h/2,
            text,
            ha="center", va="center",
            fontsize=fontsize,
            fontweight=weight,
            wrap=True
        )
        return x + w/2, y # return box center x and bottom y

    # Title
    ax.text(0.5, 0.98, "CRRT Cohort Selection", ha="center", va="center",
            fontsize=16, fontweight="bold")

    arrow_main = dict(arrowstyle="->", lw=2, color="black")

    # ------------------------------------------------------------------
    # 2. Top box (Starting cohort)
    # ------------------------------------------------------------------
    top_y = 0.90 - box_h
    # Top box is placed in the main column (left side)
    top_center_x, top_bottom_y = draw_box(
        x_main_start,
        top_y,
        box_w,
        box_h,
        "All adult hospitalizations\n"
        "(2018-2024)\n"
        f"n = {start_n:,}",
        fontsize=11
    )
    
    # For caching main column box vertical coordinates for exclusion elbows
    main_box_ys = [top_y]
    for i in range(len(steps)):
        main_box_ys.append(top_y - ((i + 1) * v_spacing))
    
    # ------------------------------------------------------------------
    # 3. Draw Rows (Remaining and Excluded)
    # ------------------------------------------------------------------
    
    # Collect the centers and box edges for exclusion elbows
    for i, step in enumerate(steps):
        
        # --- A. Calculate Y coordinates ---
        current_y = top_y - ((i + 1) * v_spacing)
        
        if i == 0:
            y_parent = top_y
        else:
            y_parent = top_y - (i * v_spacing)
        
        # Y coordinates for bottom and top of the two main boxes
        y_top_box = y_parent
        y_bottom_box = current_y
        # Vertical center (for arrow elbow): halfway between two box centers
        box_center_y_top = y_top_box + (box_h / 2)
        box_center_y_bottom = y_bottom_box + (box_h / 2)
        arrow_vertical_center = (box_center_y_top + box_center_y_bottom) / 2

        # --- B. Draw Remaining Box (Main Column) ---
        remain_text = (
            f"Remaining hospitalizations\n"
            f"{step['remaining_label']}\n"
            f"n = {step['remaining_n']:,}"
        )
        
        remain_center_x, remain_bottom_y = draw_box(
            x_main_start,
            current_y,
            box_w,
            box_h,
            remain_text,
            fontsize=11
        )
        
        # --- C. Draw Excluded Box (Right Column) ---
        if step['excluded_n'] > 0:
            excl_text = (
                f"{step['excluded_label']}\n"
                f"n = {step['excluded_n']:,}"
            )
            draw_box(
                x_excl_start,
                arrow_vertical_center - box_h / 2, # Center box vertically at elbow
                box_w,
                box_h,
                excl_text,
                fontsize=11
            )
        
        # --- D. Draw Arrows ---
        
        # 1. Main vertical flow line
        ax.annotate(
            "",
            xy=(remain_center_x, current_y + box_h),
            xytext=(remain_center_x, y_parent),
            arrowprops=arrow_main,
        )
        
        # 2. Exclusion flow line (branching from vertical line to excluded box)
        if step['excluded_n'] > 0:
            # P1: Center of the vertical line segment (center between box centers)
            p1 = (remain_center_x, arrow_vertical_center) 
            
            # The arrow now ends just before the exclusion box starts (leaving a gap)
            p2 = (x_excl_start - excl_arrow_gap, arrow_vertical_center)
            
            # Draw horizontal line segment from P1 to P2 with an arrowhead
            ax.annotate(
                "",
                xy=p2,
                xytext=p1,
                arrowprops=dict(arrowstyle="->", lw=2, color="black"),
                annotation_clip=False
            )

    # ------------------------------------------------------------------
    # 4. Add analytical note box to the PNG itself
    # ------------------------------------------------------------------
    # Fetch AKI and ICU info from strobe_counts
    n_non_esrd = strobe_counts.get('3_encounter_blocks_without_esrd', np.nan)
    n_aki_no_esrd = strobe_counts.get('6_encounter_blocks_with_AKI_no_esrd', np.nan)
    perc_aki_no_esrd = strobe_counts.get('6_Percentage_non_ESRD_encounter_blocks_with_AKI_codes', np.nan)
    n_total_crrt_hosp = strobe_counts.get('2_crrt_hospitalizations', np.nan)
    n_hosp_without_icu = strobe_counts.get('6_number_hosp_without_ICU_stay', np.nan)
    n_with_labs = get_count("6_encounter_blocks_with_required_labs", np.nan)

    print("n_with_labs", n_with_labs, n_hosp_without_icu)
    print("np.isnan(n_with_labs)",np.isnan(n_with_labs), "np.isnan(n_hosp_without_icu)", np.isnan(n_hosp_without_icu))
    if not (np.isnan(n_with_labs) or np.isnan(n_hosp_without_icu)):
        n_with_icu = int(n_with_labs) - int(n_hosp_without_icu)
        perc_with_icu = 100.0 * n_with_icu / int(n_with_labs) if n_with_labs else np.nan
        print("icu block", n_with_icu, perc_with_icu)
    else:
        n_with_icu = None
        perc_with_icu = None

    # Compose note text
    note_lines = []
    if not (np.isnan(n_non_esrd) or np.isnan(n_aki_no_esrd) or np.isnan(perc_aki_no_esrd)):
        note_lines.append(
            f"AKI codes present (non-ESRD): {int(n_aki_no_esrd):,} / {int(n_with_labs):,} ({perc_aki_no_esrd:.1f}%)"
        )
    else:
        note_lines.append("AKI data not available")
    if not (n_with_icu is None or perc_with_icu is None):
        note_lines.append(
            f"CRRT hospitalizations with ICU admission: {n_with_icu:,} / {int(n_with_labs):,} ({perc_with_icu:.1f}%)"
        )
    else:
        note_lines.append("ICU admission data not available")

    note_text = '\n'.join(note_lines)

    # Place the note as a small text box at the bottom left, under the boxes, in the figure
    # (y position: low enough to not collide; x position: left side)
    fig_height = 1.0
    y_note = 0.035  # Try bottom margin; adjust if needed
    x_note = 0.05

    # Add a subdued rectangle as background
    bbox_width = 0.89
    bbox_height = 0.10 if len(note_lines) > 2 else 0.085
    fancy_rect = FancyBboxPatch(
        (x_note, y_note - 0.03), bbox_width, bbox_height,
        boxstyle="round,pad=0.01",
        linewidth=1,
        edgecolor="#bbbbbb",
        facecolor="#f6f6ee",
        zorder=0
    )
    ax.add_patch(fancy_rect)

    # Put note text on top of the box
    ax.text(
        x_note + bbox_width/2,
        y_note + bbox_height/2 - 0.01,
        note_text,
        ha="center", va="center",
        fontsize=11, fontweight="normal",
        color="black",
        wrap=True
    )

    # ------------------------------------------------------------------
    # 5. Save and return path
    # ------------------------------------------------------------------
    consort_file = output_dir / "consort_diagram_straight_flow_right_excl.png"
    plt.savefig(consort_file, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    return consort_file

create_consort_diagram_straight_flow(strobe_counts)


# # Outcomes

# In[44]:


# ============================================================================
# OUTCOMES CALCULATION 
# ============================================================================

# 1. ICU LENGTH OF STAY 
print("\n1. Processing ICU segments...")
icu_segs = adt_final_stitched.copy()
icu_segs = icu_segs[
    (icu_segs['location_category'] == 'icu') &
    (icu_segs['in_dttm'].notna()) &
    (icu_segs['out_dttm'].notna()) &
    (icu_segs['out_dttm'] > icu_segs['in_dttm'])
]

print(f"   ICU segments identified: {len(icu_segs):,}")

# Calculate ICU LOS as sum of all ICU segment durations
icu_los = icu_segs[icu_segs['encounter_block'].isin(cohort_df['encounter_block'])].copy()
icu_los['seg_days'] = (icu_los['out_dttm'] - icu_los['in_dttm']).dt.total_seconds() / (24 * 3600)

icu_los_summary = icu_los.groupby('encounter_block').agg({
    'seg_days': 'sum'
}).reset_index()
icu_los_summary.rename(columns={'seg_days': 'icu_los_days'}, inplace=True)

print(f"   Median ICU LOS: {icu_los_summary['icu_los_days'].median():.2f} days")

# ============================================================================
# 2. HOSPITAL LENGTH OF STAY (discharge_dttm - admission_dttm)
# ============================================================================
print("\n3. Calculating Hospital Length of Stay...")
# Use admission/discharge dates from hospitalization table
# Join with encounter_mapping to get encounter_block, then aggregate for stitched encounters
hosp_dates = hosp_df[['hospitalization_id', 'admission_dttm', 'discharge_dttm']].merge(
    clif.encounter_mapping[['hospitalization_id', 'encounter_block']],
    on='hospitalization_id',
    how='inner',
)
hosp_dates_agg = hosp_dates.groupby('encounter_block').agg(
    admission_dttm=('admission_dttm', 'min'),
    discharge_dttm=('discharge_dttm', 'max'),
).reset_index()

hosp_los = cohort_df[['encounter_block']].merge(
    hosp_dates_agg,
    on='encounter_block',
    how='left'
)

# Hospital LOS = discharge_dttm - admission_dttm
hosp_los['hosp_los_days'] = (
    hosp_los['discharge_dttm'] - hosp_los['admission_dttm']
).dt.total_seconds() / (24 * 3600)

# Ensure non-negative values
hosp_los['hosp_los_days'] = hosp_los['hosp_los_days'].apply(
    lambda x: max(x, 0) if pd.notna(x) and np.isfinite(x) else np.nan
)

print(f"   Median Hospital LOS: {hosp_los['hosp_los_days'].median():.2f} days")

# ============================================================================
# 4. DEATH STATUS AND FINAL OUTCOME DATETIME
# ============================================================================
print("\n4. Determining death status and final outcome datetime...")

# Get discharge category and death_dttm from hospitalization and patient tables
patient_df = clif.patient.df[['patient_id', 'death_dttm', 'race_category', 'sex_category', 'ethnicity_category']]

death_info = cohort_df.merge(
    hosp_df[['hospitalization_id', 'patient_id', 'discharge_category', 'age_at_admission', 'admission_type_category']],
    on='hospitalization_id',
    how='left'
).merge(
    patient_df,
    on='patient_id',
    how='left'
).merge(
    vitals_range,
    on='encounter_block',
    how='left'
)

# Drop 'hospitalization_id' from death_info
if 'hospitalization_id' in death_info.columns:
    death_info = death_info.drop(columns=['hospitalization_id'])

# Collapse to unique encounter_block, aggregating required columns
death_info = death_info.sort_values('encounter_block')  

agg_dict = {
    'admission_type_category': 'last',
    'discharge_category': 'last',
    'race_category': 'last',
    'sex_category': 'last',
    'ethnicity_category': 'last',
    'death_dttm': 'last',
    'first_vital_dttm': 'min',
    'last_vital_dttm': 'max'
}

# Include all other columns not being aggregated with "first" to keep at least one value per group, unless they are non-aggregatable
for col in death_info.columns:
    if col not in agg_dict and col not in ['encounter_block']:
        agg_dict[col] = 'first'

death_info = death_info.groupby('encounter_block', as_index=False).agg(agg_dict)

# Standardize discharge category
death_info['discharge_category'] = death_info['discharge_category'].str.lower()

# Step 1: Determine if patient died (based on discharge_category)
death_info['died'] = death_info['discharge_category'].isin(['expired', 'hospice']).astype(int)

# Step 2: Determine final_outcome_dttm
# If died: use death_dttm if available, otherwise use last_vital_dttm
# If not died: use last_vital_dttm
death_info['final_outcome_dttm'] = (
    death_info['death_dttm']
    .fillna(death_info['last_vital_dttm'])  # Fallback to last_vital
    .where(death_info['died'] == 1, pd.NaT)  # Only keep for died==1, else NaT
)

print(f"   Patients identified as died (expired/hospice): {death_info['died'].sum():,}")

num_with_death_dttm = ((death_info['died'] == 1) & (death_info['death_dttm'].notna())).sum()
num_using_last_vital = ((death_info['died'] == 1) & (death_info['death_dttm'].isna())).sum()

print(f"   - With death_dttm: {num_with_death_dttm:,}")
print(f"   - Using last_vital_dttm: {num_using_last_vital:,}")

# ============================================================================
# 5. MORTALITY CALCULATIONS
# ============================================================================
print("\n5. Calculating mortality outcomes...")

# In-hospital death: died AND final_outcome_dttm is between first and last vital
death_info['in_hosp_death'] = (
    (death_info['died'] == 1) &
    (death_info['final_outcome_dttm'].notna()) &
    (death_info['final_outcome_dttm'] >= death_info['first_vital_dttm']) &
    (death_info['final_outcome_dttm'] <= death_info['last_vital_dttm'])
).astype(int)

# 30-day mortality: died AND final_outcome_dttm within 30 days of first vital
death_info['death_30d'] = (
    (death_info['died'] == 1) &
    (death_info['final_outcome_dttm'].notna()) &
    (death_info['final_outcome_dttm'] <= (death_info['first_vital_dttm'] + pd.Timedelta(days=30)))
).astype(int)


print(f"   In-hospital deaths: {death_info['in_hosp_death'].sum():,} ({death_info['in_hosp_death'].mean()*100:.1f}%)")
print(f"   30-day deaths: {death_info['death_30d'].sum():,} ({death_info['death_30d'].mean()*100:.1f}%)")

# ============================================================================
# 6. COMBINE ALL OUTCOMES
# ============================================================================
print("\n6. Combining all outcomes...")
outcomes_df = cohort_df[['hospitalization_id', 'encounter_block']].merge(
    icu_los_summary, on='encounter_block', how='left'
).merge(
    hosp_los[['encounter_block', 'hosp_los_days', 'discharge_dttm']], on='encounter_block', how='left'
).merge(
    death_info, on='encounter_block', how='left'
)

print(f"\nFinal outcomes dataset:")
print(f"   Total records: {len(outcomes_df):,}")
print(f"   Records with ICU LOS: {outcomes_df['icu_los_days'].notna().sum():,}")
print(f"   Records with Hospital LOS: {outcomes_df['hosp_los_days'].notna().sum():,}")
print(f"   In-hospital mortality rate: {outcomes_df['in_hosp_death'].mean()*100:.1f}%")
print(f"   30-day mortality rate: {outcomes_df['death_30d'].mean()*100:.1f}%")

# Display summary statistics
print("\n" + "="*60)
print("OUTCOMES SUMMARY STATISTICS")
print("="*60)
print(f"ICU LOS (days):")
print(f"  Median [IQR]: {outcomes_df['icu_los_days'].median():.1f} [{outcomes_df['icu_los_days'].quantile(0.25):.1f}-{outcomes_df['icu_los_days'].quantile(0.75):.1f}]")
print(f"\nHospital LOS (days):")
print(f"  Median [IQR]: {outcomes_df['hosp_los_days'].median():.1f} [{outcomes_df['hosp_los_days'].quantile(0.25):.1f}-{outcomes_df['hosp_los_days'].quantile(0.75):.1f}]")
print(f"\nMortality:")
print(f"  In-hospital: {outcomes_df['in_hosp_death'].sum():,}/{len(outcomes_df):,} ({outcomes_df['in_hosp_death'].mean()*100:.1f}%)")
print(f"  30-day: {outcomes_df['death_30d'].sum():,}/{len(outcomes_df):,} ({outcomes_df['death_30d'].mean()*100:.1f}%)")
print("="*60)

# Convert specified columns to lowercase (if they exist)
category_cols = [
    'admission_type_category', 'discharge_category',
    'race_category', 'sex_category', 'ethnicity_category'
]
for col in category_cols:
    if col in outcomes_df.columns:
        outcomes_df[col] = outcomes_df[col].str.lower()

# Arrange columns: patient_id, hospitalization_id, encounter_block, then everything else
front_cols = [col for col in ['patient_id', 'hospitalization_id', 'encounter_block'] if col in outcomes_df.columns]
other_cols = [col for col in outcomes_df.columns if col not in front_cols]
outcomes_df = outcomes_df[front_cols + other_cols]


# # CRRT Dose
#
# Calculate the dose for each time point , and then take the median of first 3 hours for the dose and the initiation time.

# In[45]:

crrt_cohort = crrt_cohort[crrt_cohort['encounter_block'].isin(cohort_blocks)]

if not has_crrt_settings:
    # Without CRRT settings, create index_crrt_df without dose columns
    print("CRRT settings not available — skipping dose calculation")
    index_crrt_df = crrt_initiation.merge(
        closest_weights[['encounter_block', 'weight_kg']].drop_duplicates(),
        on='encounter_block',
        how='inner'
    )
    # Filter to cohort
    index_crrt_df = index_crrt_df[index_crrt_df['encounter_block'].isin(cohort_blocks)]
    # Add hospitalization_id
    index_crrt_df = index_crrt_df.merge(
        cohort_df[['encounter_block', 'hospitalization_id']].drop_duplicates(),
        on='encounter_block',
        how='left'
    )

# --- BEGIN CRRT DOSE CALCULATION (only when has_crrt_settings=True) ---
if has_crrt_settings:
    # Define desired column order
    desired_order = [
        'hospitalization_id', 'encounter_block', 'recorded_dttm', 
        'crrt_mode_category', 'dialysate_flow_rate', 'pre_filter_replacement_fluid_rate', 
        'post_filter_replacement_fluid_rate', 'ultrafiltration_out', 
        'blood_flow_rate', 'crrt_initiation_time'
    ]
    # Only keep columns that exist in crrt_cohort
    front_cols = [col for col in desired_order if col in crrt_cohort.columns]
    other_cols = [col for col in crrt_cohort.columns if col not in front_cols]
    crrt_cohort = crrt_cohort[front_cols + other_cols]

    # Sort as specified
    crrt_cohort = crrt_cohort.sort_values(['encounter_block', 'recorded_dttm'])


    # In[46]:


    print("\n Calculating CRRT dose at initiation...")

    # Filter to dose-eligible modes (SCUF/AVVH excluded — no standard dose formula)
    DOSE_MODES = {'cvvh', 'cvvhd', 'cvvhdf'}

    print(f"   Total CRRT records before dose filtering: {len(crrt_cohort):,}")
    crrt_df_filtered = crrt_cohort[crrt_cohort['crrt_mode_category'].isin(DOSE_MODES)].copy()
    print(f"   Records after filtering for dose modes (cvvh, cvvhd, cvvhdf): {len(crrt_df_filtered):,}")
    print(f"   Excluded from dose calc: {len(crrt_cohort) - len(crrt_df_filtered):,}")

    # Fill NaN values with 0 for flow rate columns
    flow_cols = ['dialysate_flow_rate', 'pre_filter_replacement_fluid_rate',
                'post_filter_replacement_fluid_rate']

    # Drop rows where all 3 variables are missing
    crrt_df_filtered = crrt_df_filtered.dropna(subset=flow_cols, how='all')

    # Then fill remaining NaNs in those columns with 0
    crrt_df_filtered[flow_cols] = crrt_df_filtered[flow_cols].fillna(0)

    print("\n   Mode distribution across all time points:")
    print(crrt_df_filtered['crrt_mode_category'].value_counts())

    # Calculate mode-specific dose at EACH time point
    conditions = [
        crrt_df_filtered['crrt_mode_category'] == 'cvvhd',
        crrt_df_filtered['crrt_mode_category'] == 'cvvh',
        crrt_df_filtered['crrt_mode_category'] == 'cvvhdf'
    ]

    choices = [
        crrt_df_filtered['dialysate_flow_rate'],
        crrt_df_filtered['pre_filter_replacement_fluid_rate'] + crrt_df_filtered['post_filter_replacement_fluid_rate'],
        crrt_df_filtered['dialysate_flow_rate'] + crrt_df_filtered['pre_filter_replacement_fluid_rate'] +
        crrt_df_filtered['post_filter_replacement_fluid_rate']
    ]

    # Mode-specific total flow rate at each time point
    crrt_df_filtered['total_flow_rate'] = np.select(conditions, choices, default=np.nan)

    # Also calculate full flow rate (all components) at each time point
    crrt_df_filtered['total_flow_rate_full'] = (
        crrt_df_filtered['dialysate_flow_rate'] +
        crrt_df_filtered['pre_filter_replacement_fluid_rate'] +
        crrt_df_filtered['post_filter_replacement_fluid_rate']
    )

    # Merge weight data (assuming weight_df has encounter_block and weight_kg)
    crrt_df_filtered = crrt_df_filtered.merge(
        closest_weights[['encounter_block', 'weight_kg']].drop_duplicates(),
        on='encounter_block',
        how='left'
    )

    # Calculate dose at each time point
    crrt_df_filtered['crrt_dose_ml_kg_hr'] = np.where(
        (crrt_df_filtered['weight_kg'] > 0) & (crrt_df_filtered['total_flow_rate'] > 0),
        crrt_df_filtered['total_flow_rate'] / crrt_df_filtered['weight_kg'],
        np.nan
    )

    crrt_df_filtered['crrt_dose_ml_kg_hr_full'] = np.where(
        (crrt_df_filtered['weight_kg'] > 0) & (crrt_df_filtered['total_flow_rate_full'] > 0),
        crrt_df_filtered['total_flow_rate_full'] / crrt_df_filtered['weight_kg'],
        np.nan
    )

    print(f"\n   Dose calculations at individual time points:")
    print(f"     Mode-specific doses calculated: {crrt_df_filtered['crrt_dose_ml_kg_hr'].notna().sum():,}")
    print(f"     Full doses calculated: {crrt_df_filtered['crrt_dose_ml_kg_hr_full'].notna().sum():,}")


    # In[47]:


    # ============================================================================
    # Calculate Median Dose for First 3 Hours + Get Initiation Values
    # ============================================================================

    # Filter to first 3 hours after CRRT initiation
    crrt_first_3hrs = crrt_df_filtered[
        crrt_df_filtered['recorded_dttm'] <= (crrt_df_filtered['crrt_initiation_time'] + pd.Timedelta(hours=3))
    ].copy()

    print(f"   Records within first 3 hours: {len(crrt_first_3hrs):,}")

    # Define columns to aggregate
    dose_columns = ['dialysate_flow_rate', 'blood_flow_rate', 
                    'pre_filter_replacement_fluid_rate', 
                    'post_filter_replacement_fluid_rate', 'ultrafiltration_out',
                    'total_flow_rate', 'total_flow_rate_full',
                    'crrt_dose_ml_kg_hr', 'crrt_dose_ml_kg_hr_full']

    # Calculate medians for first 3 hours
    median_3hr = crrt_first_3hrs.groupby('encounter_block').agg({
        'hospitalization_id': 'first',
        'crrt_initiation_time': 'first',
        'weight_kg': 'first',
        'crrt_mode_category': lambda x: x.mode()[0] if not x.empty else np.nan,  # Most frequent mode
        **{col: 'median' for col in dose_columns}
    }).reset_index()

    print(f"   Encounters with median values calculated: {len(median_3hr):,}")

    # Get values at initiation time (original values)
    print("\n   Getting original values at initiation time...")

    crrt_at_init = crrt_df_filtered[
        crrt_df_filtered['recorded_dttm'] == crrt_df_filtered['crrt_initiation_time']
    ].copy()

    # Group by encounter_block and take first (should be unique at initiation time)
    init_values = crrt_at_init.groupby('encounter_block').agg({
        col: 'first' for col in dose_columns
    }).reset_index()

    # Rename init columns to add _not_avged suffix
    init_values = init_values.rename(columns={
        col: f'{col}_not_avged' for col in dose_columns
    })

    print(f"   Encounters with initiation values: {len(init_values):,}")

    # Merge median and initiation values
    final_df = median_3hr.merge(init_values, on='encounter_block', how='left')

    # Now arrange columns in the requested order
    # First the main columns (with median values)
    main_columns = [
        'encounter_block', 'hospitalization_id', 'crrt_initiation_time',
        'weight_kg', 'crrt_mode_category',
        'dialysate_flow_rate', 'blood_flow_rate',
        'pre_filter_replacement_fluid_rate', 'post_filter_replacement_fluid_rate',
        'ultrafiltration_out', 'total_flow_rate', 'total_flow_rate_full',
        'crrt_dose_ml_kg_hr', 'crrt_dose_ml_kg_hr_full'
    ]

    # Then the initiation columns (with _not_avged suffix)
    not_avged_columns = [f'{col}_not_avged' for col in dose_columns]

    # Combine all columns
    all_columns = main_columns + not_avged_columns

    # Select and reorder columns
    final_df = final_df[all_columns]

    print(f"\n   Final dataframe created:")
    print(f"     Total rows: {len(final_df):,} (one per encounter)")
    print(f"     Total columns: {len(final_df.columns)}")

    # Assign to your desired variable name
    index_crrt_df = final_df.copy()

    print("\n✅ Final dataframe created with one row per encounter!")
    print(f"   Stored as 'index_crrt_df' with {len(index_crrt_df)} encounters")


    # In[48]:


    # ============================================================================
    # Visualization: Overlaid Histograms of Median vs Initiation Doses
    # ============================================================================
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np

    print("\n" + "=" * 80)
    print("Creating Overlaid Histogram Visualization")
    print("=" * 80)

    # Set style for better visualization
    plt.style.use('seaborn-v0_8-darkgrid')
    sns.set_palette("husl")

    # Create a single figure
    fig, ax = plt.subplots(figsize=(10, 6))

    # Remove NaN values for cleaner plotting
    dose_full_median = final_df['crrt_dose_ml_kg_hr_full'].dropna()
    dose_full_init = final_df['crrt_dose_ml_kg_hr_full_not_avged'].dropna()

    # Create overlaid histograms
    ax.hist(dose_full_median, bins=30, alpha=0.5, label='Median (First 3hr)',
            color='blue', edgecolor='darkblue', density=True)
    ax.hist(dose_full_init, bins=30, alpha=0.5, label='At Initiation',
            color='red', edgecolor='darkred', density=True)

    # Add vertical lines for means
    ax.axvline(dose_full_median.mean(), color='blue', linestyle='--',
                linewidth=2, label=f'Mean (3hr): {dose_full_median.mean():.1f}')
    ax.axvline(dose_full_init.mean(), color='red', linestyle='--',
                linewidth=2, label=f'Mean (Init): {dose_full_init.mean():.1f}')

    # Labels and title
    ax.set_xlabel('CRRT Dose (mL/kg/hr)', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title('Full CRRT Dose: Median (3hr) vs Initiation', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # Save figure
    plt.tight_layout()
    fig.savefig("../output/final/dose_comparison.png")
    plt.close(fig)


    # In[49]:


    # ============================================================================
    # Plot 2: Mode-Specific CRRT Dose Comparison by Mode Category
    # ============================================================================
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np

    # Get unique mode categories (excluding NaN)
    mode_categories = final_df['crrt_mode_category'].dropna().unique()
    mode_categories = sorted(mode_categories)  # Sort for consistent ordering

    # Create figure with subplots for each mode
    n_modes = len(mode_categories)
    fig, axes = plt.subplots(1, n_modes, figsize=(6*n_modes, 6))

    # If only one mode, make axes iterable
    if n_modes == 1:
        axes = [axes]

    # Color palette
    colors_median = ['steelblue', 'royalblue', 'dodgerblue']
    colors_init = ['coral', 'salmon', 'lightsalmon']

    # Process each mode
    for idx, mode in enumerate(mode_categories):
        ax = axes[idx]

        # Filter data for this mode
        mode_data = final_df[final_df['crrt_mode_category'] == mode]

        # Get dose values for this mode
        dose_median = mode_data['crrt_dose_ml_kg_hr'].dropna()
        dose_init = mode_data['crrt_dose_ml_kg_hr_not_avged'].dropna()

        # Skip if no data
        if len(dose_median) == 0 and len(dose_init) == 0:
            ax.text(0.5, 0.5, f'No data for {mode.upper()}',
                    ha='center', va='center', fontsize=12)
            ax.set_title(f'{mode.upper()}', fontsize=14, fontweight='bold')
            continue

        # Create overlaid histograms
        if len(dose_median) > 0:
            ax.hist(dose_median, bins=20, alpha=0.5, label=f'Median 3hr (n={len(dose_median)})',
                    color=colors_median[idx % len(colors_median)], edgecolor='darkblue', density=True)
            ax.axvline(dose_median.mean(), color=colors_median[idx % len(colors_median)],
                        linestyle='--', linewidth=2, label=f'Mean 3hr: {dose_median.mean():.1f}')

        if len(dose_init) > 0:
            ax.hist(dose_init, bins=20, alpha=0.5, label=f'At Init (n={len(dose_init)})',
                    color=colors_init[idx % len(colors_init)], edgecolor='darkred', density=True)
            ax.axvline(dose_init.mean(), color=colors_init[idx % len(colors_init)],
                        linestyle='--', linewidth=2, label=f'Mean Init: {dose_init.mean():.1f}')

        # Labels and title
        ax.set_xlabel('CRRT Dose (mL/kg/hr)', fontsize=11)
        ax.set_ylabel('Density' if idx == 0 else '', fontsize=11)
        ax.set_title(f'{mode.upper()}\n(n={len(mode_data)} encounters)', fontsize=13, fontweight='bold')
        ax.legend(loc='upper right', fontsize=9)
        ax.grid(True, alpha=0.3)

    # Overall title
    fig.suptitle('Mode-Specific CRRT Dose Comparison by Mode Category\nMedian (First 3hr) vs At Initiation',
                fontsize=15, fontweight='bold', y=1.02)

    plt.tight_layout()

    # Save figure
    output_path = '../output/final/graphs/crrt_dose_comparison_by_mode.png'
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"\n✓ Saved mode-specific comparison to: {output_path}")


    # --- END CRRT DOSE CALCULATION ---

# # CRRT duration

# In[50]:


import pandas as pd
import numpy as np

def calculate_crrt_duration(crrt_cohort):
    """
    Calculate CRRT duration for each encounter.
    Duration is defined as the time from crrt_initiation_time to the last recorded setting,
    considering CRRT ended when there's a 24-hour gap in recordings.
    
    Parameters:
    -----------
    crrt_cohort : pd.DataFrame
        DataFrame with CRRT time series data including columns:
        - encounter_block: patient identifier
        - crrt_initiation_time: start of CRRT
        - recorded_dttm: timestamp column for each setting recording
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with encounter_block and calculated duration metrics
    """

    # Make a copy to avoid modifying original
    df = crrt_cohort.copy()

    # Ensure datetime columns are properly formatted
    time_column = 'recorded_dttm'
    df[time_column] = pd.to_datetime(df[time_column])
    df['crrt_initiation_time'] = pd.to_datetime(df['crrt_initiation_time'])

    # CRITICAL: Drop all rows where recorded_dttm is before crrt_initiation_time
    initial_rows = len(df)
    df = df[df[time_column] >= df['crrt_initiation_time']]
    rows_dropped = initial_rows - len(df)

    if rows_dropped > 0:
        print(f"Dropped {rows_dropped} recordings that occurred before CRRT initiation time")
        print(f"Remaining recordings: {len(df)}")

    # Sort by encounter and time
    df = df.sort_values(['encounter_block', time_column])

    # Function to calculate duration for each encounter
    def get_encounter_duration(group):
        """Calculate CRRT duration for a single encounter with 24-hour gap detection"""

        # Get initiation time
        initiation_time = group['crrt_initiation_time'].iloc[0]

        # Get all recorded times (already filtered to be >= initiation_time)
        recorded_times = group[time_column].dropna().sort_values()

        if len(recorded_times) == 0:
            # No recordings after initiation
            return pd.Series({
                'crrt_initiation_time': initiation_time,
                'crrt_end_time': initiation_time,
                'duration_hours': 0,
                'duration_days': 0,
                'num_recordings': 0,
                'had_24hr_gap': False
            })

        # Check for 24-hour gaps
        time_diffs = recorded_times.diff()

        # Find if there's any gap >= 24 hours
        gaps_24hr = time_diffs > pd.Timedelta(hours=24)

        if gaps_24hr.any():
            # CRRT ended at the last recording before the first 24-hour gap
            first_gap_idx = gaps_24hr.idxmax()
            # Get the index position of the gap
            gap_position = recorded_times.index.get_loc(first_gap_idx)
            # The end time is the recording just before the gap
            end_time = recorded_times.iloc[gap_position - 1]
            had_gap = True
        else:
            # No 24-hour gap, use the last recording
            end_time = recorded_times.iloc[-1]
            had_gap = False

        # Calculate duration
        duration = end_time - initiation_time
        duration_hours = duration.total_seconds() / 3600
        duration_days = duration_hours / 24

        # Count recordings
        num_recordings = len(recorded_times)

        return pd.Series({
            'crrt_initiation_time': initiation_time,
            'crrt_end_time': end_time,
            'duration_hours': duration_hours,
            'duration_days': duration_days,
            'num_recordings': num_recordings,
            'had_24hr_gap': had_gap
        })

    # Apply to each encounter
    duration_df = df.groupby('encounter_block').apply(get_encounter_duration).reset_index()

    # Add summary statistics
    print("\n=== CRRT Duration Summary ===")
    print(f"Total encounters: {len(duration_df)}")
    print(f"Encounters with recordings: {len(duration_df[duration_df['num_recordings'] > 0])}")
    print(f"Encounters without recordings: {len(duration_df[duration_df['num_recordings'] == 0])}")

    # Stats for encounters with recordings
    valid_durations = duration_df[duration_df['duration_hours'] > 0]

    if len(valid_durations) > 0:
        print(f"\nDuration Statistics (for {len(valid_durations)} encounters with valid recordings):")
        print(f"\nDuration (hours):")
        print(f"  Mean: {valid_durations['duration_hours'].mean():.1f}")
        print(f"  Median: {valid_durations['duration_hours'].median():.1f}")
        print(f"  Q25: {valid_durations['duration_hours'].quantile(0.25):.1f}")
        print(f"  Q75: {valid_durations['duration_hours'].quantile(0.75):.1f}")
        print(f"  Min: {valid_durations['duration_hours'].min():.1f}")
        print(f"  Max: {valid_durations['duration_hours'].max():.1f}")

        print(f"\nDuration (days):")
        print(f"  Mean: {valid_durations['duration_days'].mean():.1f}")
        print(f"  Median: {valid_durations['duration_days'].median():.1f}")
        print(f"  Q25: {valid_durations['duration_days'].quantile(0.25):.1f}")
        print(f"  Q75: {valid_durations['duration_days'].quantile(0.75):.1f}")

    print(f"\nEncounters with 24-hour gap: {duration_df['had_24hr_gap'].sum()} ({duration_df['had_24hr_gap'].mean()*100:.1f}%)")

    # Add duration categories
    duration_df['duration_category'] = pd.cut(
        duration_df['duration_days'],
        bins=[-0.001, 0, 1, 3, 7, 14, float('inf')],
        labels=['No duration', '<1 day', '1-3 days', '3-7 days', '7-14 days', '>14 days'],
        include_lowest=True
    )

    print(f"\nDuration categories:")
    print(duration_df['duration_category'].value_counts().sort_index())

    return duration_df

# Pre-processing function to check for pre-initiation recordings
def check_pre_initiation_recordings(crrt_cohort):
    """
    Check how many recordings occur before CRRT initiation time
    """
    df = crrt_cohort.copy()
    df['recorded_dttm'] = pd.to_datetime(df['recorded_dttm'])
    df['crrt_initiation_time'] = pd.to_datetime(df['crrt_initiation_time'])

    # Find pre-initiation recordings
    pre_init = df[df['recorded_dttm'] < df['crrt_initiation_time']]

    if len(pre_init) > 0:
        print("=== Pre-Initiation Recordings Found ===")
        print(f"Total pre-initiation recordings: {len(pre_init)} ({len(pre_init)/len(df)*100:.1f}% of all recordings)")
        print(f"Affected encounters: {pre_init['encounter_block'].nunique()}")

        # Calculate how early these recordings are
        pre_init['hours_before'] = (pre_init['crrt_initiation_time'] - pre_init['recorded_dttm']).dt.total_seconds() / 3600

        print(f"\nTiming statistics (hours before initiation):")
        print(f"  Mean: {pre_init['hours_before'].mean():.1f} hours")
        print(f"  Median: {pre_init['hours_before'].median():.1f} hours")
        print(f"  Max: {pre_init['hours_before'].max():.1f} hours")

        # Show a few examples
        print("\nExample pre-initiation recordings:")
        sample = pre_init.nlargest(5, 'hours_before')[['encounter_block', 'recorded_dttm', 'crrt_initiation_time', 'hours_before']]
        print(sample)
    else:
        print("No pre-initiation recordings found - data is clean!")

    return pre_init

# First check for pre-initiation recordings (optional)
pre_init_check = check_pre_initiation_recordings(crrt_cohort)

# Calculate duration (automatically drops pre-initiation recordings)
duration_df = calculate_crrt_duration(crrt_cohort)

# Merge back with main cohort
index_crrt_df = index_crrt_df.merge(
    duration_df[['encounter_block', 'duration_hours', 'duration_days', 'duration_category', 'had_24hr_gap']],
    on='encounter_block',
    how='left'
)


# # Respiratory support
# 
# For duration on IMV, start time as the first IMV row, and end time as not on IMV for 24 hours

# In[51]:


print(f"\nLoading respiratory support table...")
try:
    clif.load_table(
        'respiratory_support',
        filters={'hospitalization_id': list(cohort_df["hospitalization_id"].unique())}
    )
    print(f"   respiratory_support loaded: {len(clif.respiratory_support.df):,} rows")
    print(f"   Unique respiratory_support hospitalizations: {clif.respiratory_support.df['hospitalization_id'].nunique()}")
except Exception as e:
    print(f"   respiratory_support not available or error: {e}")
    raise


# In[52]:


clif.respiratory_support.df = clif.respiratory_support.df.merge(
    clif.encounter_mapping[['hospitalization_id', 'encounter_block']],
    on='hospitalization_id',
    how='left'
)
clif.respiratory_support = clif.respiratory_support.waterfall()


# In[53]:


resp_support_df = clif.respiratory_support.df
# del clif


# In[54]:


resp_support_df.columns


# In[55]:


def calculate_imv_duration(resp_support_df):
    """
    Calculate IMV duration for each encounter block.
    Uses device_category transitions (after waterfall) to identify continuous
    IMV stretches. A stretch ends when a non-IMV device is recorded.
    Total IMV duration is the sum of all IMV stretches.

    Parameters:
    -----------
    resp_support_df : pd.DataFrame
        DataFrame with columns: encounter_block, recorded_dttm, device_category
        (should have waterfall already applied)

    Returns:
    --------
    pd.DataFrame
        DataFrame with encounter_block and IMV duration metrics
    """
    df = resp_support_df.copy()
    df['recorded_dttm'] = pd.to_datetime(df['recorded_dttm'])
    df = df.sort_values(['encounter_block', 'recorded_dttm'])

    # Mark IMV rows using device_category
    df['is_imv'] = df['device_category'].str.lower() == 'imv'

    # Assign stretch IDs: increment whenever is_imv status changes
    df['device_change'] = df.groupby('encounter_block')['is_imv'].transform(
        lambda x: (x != x.shift()).cumsum()
    )

    # Filter to only IMV stretches
    imv_stretches = df[df['is_imv']].copy()

    if len(imv_stretches) == 0:
        return pd.DataFrame({
            'encounter_block': pd.Series(dtype='int'),
            'imv_start_time': pd.Series(dtype='datetime64[ns]'),
            'imv_end_time': pd.Series(dtype='datetime64[ns]'),
            'imv_duration_hours': pd.Series(dtype='float'),
            'imv_duration_days': pd.Series(dtype='float'),
        })

    # Per stretch: duration = last recording - first recording
    stretch_stats = imv_stretches.groupby(
        ['encounter_block', 'device_change']
    )['recorded_dttm'].agg(['first', 'last']).reset_index()
    stretch_stats['stretch_hours'] = (
        (stretch_stats['last'] - stretch_stats['first']).dt.total_seconds() / 3600
    )

    # Per encounter: sum all stretch durations, keep overall first/last times
    imv_duration_df = stretch_stats.groupby('encounter_block').agg(
        imv_start_time=('first', 'min'),
        imv_end_time=('last', 'max'),
        imv_duration_hours=('stretch_hours', 'sum'),
    ).reset_index()
    imv_duration_df['imv_duration_days'] = imv_duration_df['imv_duration_hours'] / 24

    return imv_duration_df

# Usage:
imv_duration_df = calculate_imv_duration(resp_support_df)


# In[56]:


imv_duration_df.columns


# In[57]:


# Merge back with main cohort
index_crrt_df = index_crrt_df.merge(
    imv_duration_df[['encounter_block', 'imv_duration_hours', 'imv_duration_days']],
    on='encounter_block',
    how='left'
)


# # Save Intermediate data

# In[58]:


cohort_df.to_parquet("../output/intermediate/cohort_df.parquet", index=False)
outcomes_df.to_parquet("../output/intermediate/outcomes_df.parquet", index=False)
# Filter weight_df to hospitalization_ids present in cohort_df before saving
weight_df_filtered = weight_df[weight_df["hospitalization_id"].isin(cohort_df["hospitalization_id"])]
weight_df_filtered.to_parquet("../output/intermediate/weight_df.parquet", index=False)
# save or use crrt_initiation df 
crrt_initiation.to_parquet("../output/intermediate/crrt_initiation.parquet", index=False)
crrt_at_initiation.to_parquet("../output/intermediate/crrt_at_initiation.parquet", index=False)
index_crrt_df.to_parquet("../output/intermediate/index_crrt_df.parquet", index=False)
crrt_cohort.to_parquet("../output/intermediate/crrt_cohort.parquet", index=False)


# In[ ]:




