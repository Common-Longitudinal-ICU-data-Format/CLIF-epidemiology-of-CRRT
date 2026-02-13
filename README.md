# CLIF Epidemiology of CRRT

Multi-site analysis of Continuous Renal Replacement Therapy (CRRT) for Acute Kidney Injury (AKI) across the CLIF consortium.

**CLIF Version:** 2.1.0

## Objective

Describe the epidemiology of CRRT for AKI across CLIF consortium hospitals:
- Frequency of CRRT (denominator: all ICU admissions without ESRD)
- Pre-CRRT clinical state (BMP, lactate, vasopressors, IMV, location)
- Treatment details (CRRT mode, settings, dose, duration)
- Outcomes (in-hospital mortality)

## Required CLIF Tables

| Table | Required Columns | Required Categories |
|-------|-----------------|---------------------|
| **clif_patient** | `patient_id`, `race_category`, `ethnicity_category`, `sex_category`, `death_dttm` | - |
| **clif_hospitalization** | `patient_id`, `hospitalization_id`, `admission_dttm`, `discharge_dttm`, `age_at_admission`, `discharge_category` | - |
| **clif_adt** | `hospitalization_id`, `hospital_id`, `in_dttm`, `out_dttm`, `location_category`, `location_type` | - |
| **clif_vitals** | `hospitalization_id`, `recorded_dttm`, `vital_category`, `vital_value` | heart_rate, resp_rate, sbp, dbp, map, spo2, weight_kg, height_cm |
| **clif_labs** | `hospitalization_id`, `lab_result_dttm`, `lab_category`, `lab_value`, `lab_value_numeric` | sodium, potassium, chloride, bicarbonate, bun, creatinine, glucose_serum, calcium_total, lactate, magnesium, ph_arterial, ph_venous, po2_arterial |
| **clif_medication_admin_continuous** | `hospitalization_id`, `admin_dttm`, `med_name`, `med_category`, `med_dose`, `med_dose_unit` | norepinephrine, epinephrine, phenylephrine, vasopressin, dopamine, angiotensin, dobutamine, milrinone, isoproterenol |
| **clif_respiratory_support** | `hospitalization_id`, `recorded_dttm`, `device_category`, `mode_category`, `tracheostomy`, `fio2_set`, `peep_set`, `resp_rate_set`, `tidal_volume_set` | - |
| **clif_crrt_therapy** | `hospitalization_id`, `recorded_dttm` (+ settings columns if `has_crrt_settings=true`) | - |
| **clif_hospital_diagnosis** | `hospitalization_id`, `diagnosis_code`, `present_on_admission` | - |



## Prerequisites

- **Python 3.11+**
- **UV package manager** ([install here](https://docs.astral.sh/uv/))
- Access to CLIF 2.1.0 data tables at your site

## Quick Start

### macOS / Linux

```bash
# 1. Clone the repository
git clone <repo-url>
cd CLIF-epidemiology-of-CRRT

# 2. Create your site config
cp config/config_template.json config/config.json
# Edit config/config.json with your site settings (see Configuration below)

# 3. Install dependencies
uv sync

# 4. Run the pipeline (all 4 steps in order)
bash run_pipeline.sh
# Or run individually:
# uv run python code/00_cohort.py
# uv run python code/01_create_wide_df.py
# uv run python code/02_construct_crrt_tableone.py
# uv run python code/03_crrt_visualizations.py
```

### Windows

```cmd
:: 1. Clone the repository
git clone <repo-url>
cd CLIF-epidemiology-of-CRRT

:: 2. Create your site config
copy config\config_template.json config\config.json
:: Edit config\config.json with your site settings (see Configuration below)

:: 3. Install dependencies
uv sync

:: 4. Run the pipeline (all 4 steps in order)
run_pipeline.bat
:: Or run individually:
:: set PYTHONIOENCODING=utf-8
:: uv run python code\00_cohort.py
:: uv run python code\01_create_wide_df.py
:: uv run python code\02_construct_crrt_tableone.py
:: uv run python code\03_crrt_visualizations.py
```

> **Note:** The batch script sets `PYTHONIOENCODING=utf-8` automatically. If running scripts individually on Windows, set this variable first to avoid Unicode encoding errors.

## Configuration

Edit `config/config.json` (copied from `config/config_template.json`):

```json
{
    "site_name": "Your_Site_Name",
    "tables_path": "/path/to/clif/tables/",
    "file_type": "parquet",
    "timezone": "US/Central",
    "project_root": "/path/to/CLIF-epidemiology-of-CRRT",
    "has_crrt_settings": true
}
```

| Field | Description |
|-------|-------------|
| `site_name` | Your hospital/site identifier (e.g., "UCMC", "MIMIC") |
| `tables_path` | Absolute path to directory containing CLIF 2.1.0 tables (use forward slashes on all platforms, e.g., `Z:/data/clif`) |
| `file_type` | File format of CLIF tables: `"parquet"`, `"csv"`, or `"fst"` |
| `timezone` | Timezone for your data (e.g., `"US/Central"`, `"US/Eastern"`) |
| `project_root` | Absolute path to this repository's root directory |
| `has_crrt_settings` | Set to `true` if your `clif_crrt_therapy` table includes flow rate and mode columns. Set to `false` if it only has `hospitalization_id` and `recorded_dttm`. |

### `has_crrt_settings` Details

- **`true`** (default): Your CRRT table has `crrt_mode_category`, `blood_flow_rate`, `dialysate_flow_rate`, `pre_filter_replacement_fluid_rate`, `post_filter_replacement_fluid_rate`, `ultrafiltration_out`. The pipeline will compute CRRT dose, generate dose visualizations, and include per-mode statistics in Table 1.
- **`false`**: Your CRRT table only has `hospitalization_id` and `recorded_dttm`. The pipeline will skip dose calculations, dose visualizations, and CRRT settings rows in Table 1. CRRT initiation is defined as the first `crrt_therapy` record.

## Pipeline Steps

### Step 0: Cohort Identification (`00_cohort.py`)
Identifies the CRRT cohort with inclusion/exclusion criteria, computes CRRT initiation times, outcomes (mortality, LOS), and IMV duration.

**Outputs:** `output/intermediate/outcomes_df.parquet`, `output/intermediate/index_crrt_df.parquet`, `output/intermediate/crrt_initiation.parquet`, `output/final/strobe_counts.csv`

### Step 1: Wide Dataset (`01_create_wide_df.py`)
Builds a wide longitudinal dataset merging labs, vitals, medications, respiratory support, CRRT therapy, and ADT into a single time-indexed table. Includes vasopressor unit conversion and NEE computation.

**Outputs:** `output/intermediate/wide_df.parquet`

### Step 2: Table 1 (`02_construct_crrt_tableone.py`)
Generates Table 1 with demographics, SOFA scores, labs, vasopressors, respiratory settings, and CRRT details across time windows (baseline, 72h, discharge) stratified by survival.

**Outputs:** `output/final/table1_crrt.csv`, `output/final/table1_crrt.html`

### Step 3: Visualizations (`03_crrt_visualizations.py`)
Generates clinical trajectory figures: CRRT dose over time, lab distributions, MAP, respiratory support, and patient state (stacked area plot).

**Outputs:** `output/final/graphs/` (PNG figures + CSV aggregate data)

