# CLIF Epidemiology of CRRT

Multi-site analysis of Continuous Renal Replacement Therapy (CRRT) for Acute Kidney Injury (AKI) across the CLIF consortium.

**CLIF Version:** 2.1.0

## Objectives

1. **Descriptive epidemiology** of CRRT for AKI across CLIF consortium hospitals:
   - Frequency of CRRT (denominator: all ICU admissions without ESRD)
   - Pre-CRRT clinical state (BMP, lactate, vasopressors, IMV, location)
   - Treatment details (CRRT mode, settings, dose, duration)
   - Outcomes (in-hospital mortality, LOS, IMV duration)

2. **Causal inference** — effect of CRRT dose on 30-day mortality:
   - Point-treatment analysis (PSM, IPTW via Super Learner)
   - Continuous dose-response (natural spline Cox, generalized propensity score)
   - Time-varying marginal structural models (12h and 24h intervals)
   - Hypothesis-generating subgroup analysis (9 clinically motivated subgroups)

## Required CLIF Tables

| Table | Required Columns | Required Categories |
|-------|-----------------|---------------------|
| **clif_patient** | `patient_id`, `race_category`, `ethnicity_category`, `sex_category`, `death_dttm` | - |
| **clif_hospitalization** | `patient_id`, `hospitalization_id`, `admission_dttm`, `discharge_dttm`, `age_at_admission`, `discharge_category` | - |
| **clif_adt** | `hospitalization_id`, `hospital_id`, `in_dttm`, `out_dttm`, `location_category`, `location_type` | - |
| **clif_vitals** | `hospitalization_id`, `recorded_dttm`, `vital_category`, `vital_value` | heart_rate, resp_rate, sbp, dbp, map, spo2, weight_kg, height_cm |
| **clif_labs** | `hospitalization_id`, `lab_result_dttm`, `lab_category`, `lab_value`, `lab_value_numeric` | sodium, potassium, chloride, bicarbonate, bun, creatinine, glucose_serum, calcium_total, lactate, magnesium, ph_arterial, ph_venous, po2_arterial |
| **clif_medication_admin_continuous** | `hospitalization_id`, `admin_dttm`, `med_name`, `med_category`, `med_dose`, `med_dose_unit` | norepinephrine, epinephrine, phenylephrine, vasopressin, dopamine, angiotensin, dobutamine, milrinone, isoproterenol |
| **clif_medication_admin_intermittent** | `hospitalization_id`, `admin_dttm`, `med_category` | - |
| **clif_respiratory_support** | `hospitalization_id`, `recorded_dttm`, `device_category`, `mode_category`, `tracheostomy`, `fio2_set`, `peep_set`, `resp_rate_set`, `tidal_volume_set` | - |
| **clif_crrt_therapy** | `hospitalization_id`, `recorded_dttm` (+ settings columns if `has_crrt_settings=true`) | - |
| **clif_hospital_diagnosis** | `hospitalization_id`, `diagnosis_code`, `present_on_admission` | - |
| **clif_microbiology_culture** | `hospitalization_id`, `collected_dttm`, `result_category` | - |

See `config/clif_data_requirements.yaml` for the full column and category specification.

## Prerequisites

- **Python 3.11+**
- **R 4.x** (for causal inference scripts 05-06b)
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

# 4. Run the full pipeline (descriptive 00-03, causal 04-06b)
bash run_pipeline.sh
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

:: 4. Run the full pipeline (descriptive 00-03, causal 04-06b)
run_pipeline.bat
```

> **Note:** The pipeline scripts run all steps sequentially: Python scripts 00-04 followed by R scripts 05-06b. R scripts are invoked with `--no-init-file` to avoid conflicts with user `.Rprofile` settings. CRAN mirror defaults to `https://cloud.r-project.org` if not configured.

### R Packages

The R scripts auto-install missing CRAN packages on first run. Key dependencies include:

`tidyverse`, `arrow`, `survival`, `cmprsk`, `MatchIt`, `WeightIt`, `cobalt`, `SuperLearner`, `randomForest`, `xgboost`, `gam`, `EValue`, `survminer`, `survey`, `mice`, `gtsummary`, `patchwork`, `splines`

## Configuration

Edit `config/config.json` (copied from `config/config_template.json`):

```json
{
    "site_name": "Your_Site_Name",
    "tables_path": "/path/to/clif/tables/",
    "file_type": "parquet",
    "timezone": "US/Central",
    "project_root": "/path/to/CLIF-epidemiology-of-CRRT",
    "has_crrt_settings": true,
    "admission_year_start": 2018,
    "admission_year_end": null
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
| `admission_year_start` | (Optional) Filter to admissions on or after this year |
| `admission_year_end` | (Optional) Filter to admissions before this year. `null` for no upper bound. |

### `has_crrt_settings` Details

- **`true`** (default): Your CRRT table has `crrt_mode_category`, `blood_flow_rate`, `dialysate_flow_rate`, `pre_filter_replacement_fluid_rate`, `post_filter_replacement_fluid_rate`, `ultrafiltration_out`. The pipeline will compute CRRT dose, generate dose visualizations, and include per-mode statistics in Table 1.
- **`false`**: Your CRRT table only has `hospitalization_id` and `recorded_dttm`. The pipeline will skip dose calculations, dose visualizations, and CRRT settings rows in Table 1. CRRT initiation is defined as the first `crrt_therapy` record.

## Pipeline Steps

### Descriptive Epidemiology (Steps 00-03)

#### Step 00: Cohort Identification (`00_cohort.py`)
Identifies the CRRT cohort with inclusion/exclusion criteria (excludes ESRD via ICD codes), computes CRRT initiation times, outcomes (mortality, LOS), and IMV duration. Generates a CONSORT flow diagram.

**Outputs:**
- `output/intermediate/outcomes_df.parquet`
- `output/intermediate/index_crrt_df.parquet`
- `output/intermediate/crrt_initiation.parquet`
- `output/final/crrt_epi/{site}_strobe_counts.csv`
- `output/final/crrt_epi/graphs/{site}_consort_diagram_straight_flow_right_excl.png`

#### Step 01: Wide Dataset (`01_create_wide_df.py`)
Builds a wide longitudinal dataset merging labs, vitals, medications, respiratory support, CRRT therapy, and ADT into a single time-indexed table. Includes vasopressor unit conversion and norepinephrine-equivalent (NEE) computation.

**Outputs:** `output/intermediate/wide_df.parquet`

#### Step 02: Table 1 (`02_construct_crrt_tableone.py`)
Generates Table 1 with demographics, SOFA scores, labs, vasopressors, respiratory settings, and CRRT details across time windows (baseline -12h to +3h, 72h post-CRRT, discharge) stratified by survival.

**Outputs:**
- `output/final/crrt_epi/{site}_table1_crrt.csv`
- `output/final/crrt_epi/{site}_table1_crrt.html`
- `output/final/crrt_epi/{site}_table1_crrt_long.csv`

#### Step 03: Visualizations (`03_crrt_visualizations.py`)
Generates clinical trajectory figures for the first 7 days post-CRRT: dose over time, lab distributions, MAP, respiratory support, and patient state (stacked area plot).

**Outputs:** `output/final/crrt_epi/graphs/{site}_*.png` (figures) + `{site}_*.csv` (aggregate data)

### Causal Inference (Steps 04-06b)

> These steps require the descriptive pipeline (steps 00-03) to have completed first.

#### Step 04: Competing Risk DataFrame & Causal CONSORT (`04_build_msm_competing_risk_df.py`)
Builds a wide competing-risk DataFrame (58 columns) with labs, SOFA scores, oxygenation, vasopressors, IMV, and Charlson Comorbidity Index at t=0/12/24h windows. Includes sensitivity columns for 24h/48h interval analyses. Generates a missingness heatmap and a CONSORT flow diagram showing cohort narrowing from descriptive analysis through causal eligibility criteria.

**Outputs:**
- `output/intermediate/msm_competing_risk_df.parquet`
- `output/final/crrt_epi/graphs/{site}_missingness_heatmap.png`
- `output/final/psm_iptw/{site}_causal_consort_diagram.{png,pdf}`

#### Step 05: PSM & IPTW (`05_PSM_IPTW_CRRT_dose.R`)
Point-treatment causal analysis of CRRT dose (high >=30 vs low <30 mL/kg/hr) on 30-day mortality. Three branches:
- **PSM**: Nearest-neighbor matching + Fine-Gray competing risk models
- **IPTW**: Super Learner propensity scores + cause-specific Cox models
- **Subgroup analysis**: 9 clinically motivated subgroups with interaction tests

Uses MICE imputation (5 datasets), Rubin's rules pooling, bootstrap CIF curves (500 reps), and E-value sensitivity analysis.

**Outputs:** `output/final/psm_iptw/` — balance plots, CIF curves, model comparison CSV, subgroup forest plot, E-value table, pooled results, Table 1/S1/S2 CSVs, CIF data CSVs. See `psm_iptw_output_guide.md` in that directory for a full listing.

#### Step 05b: Dose-Response Analysis (`05b_dose_response_analysis.R`)
Continuous dose-response analysis with three approaches:
- **Dose decile**: Unadjusted mortality by dose decile
- **Natural spline Cox**: Covariate-adjusted nonlinear dose-response (nonlinearity test)
- **Generalized propensity score**: Doubly-robust causal analysis via CBPS (continuous treatment)

Also generates a target trial emulation specification table and a 24h exclusion sensitivity analysis.

**Outputs:** `output/final/psm_iptw/` — dose-response plots, combined 3-panel figure, GPS balance diagnostics, target trial CSV, dose decile CSV

#### Step 06: Time-Varying MSM (`06_time_varying_MSM.R`)
Time-varying marginal structural model for CRRT dose using 12h intervals. Models dose trajectories as treatment strategies with inverse-probability-of-treatment weighting.

**Outputs:** `output/final/msm/time_varying/` — dose histograms, balance tables, CIF curves, cause-specific Cox results, ESS diagnostics, Table 1/S1 CSVs, CIF data CSV

#### Step 06b: MSM Sensitivity (`06b_time_varying_MSM_sensitivity.R`)
Sensitivity analysis using 24h intervals with 48h eligibility filter (survived >=48h, CRRT >=48h, has dialysate in 0-48h window).

**Outputs:** `output/final/msm/time_varying_sensitivity/` — same structure as primary MSM

## Output Structure

```
output/
├── intermediate/                    # Intermediate parquet files (not committed)
│   ├── outcomes_df.parquet
│   ├── index_crrt_df.parquet
│   ├── crrt_initiation.parquet
│   ├── wide_df.parquet
│   ├── tableone_analysis_df.parquet
│   └── msm_competing_risk_df.parquet
│
└── final/                           # All files are site-prefixed (e.g., UCMC_*)
    ├── crrt_epi/                    # Descriptive epidemiology
    │   ├── {site}_strobe_counts.csv
    │   ├── {site}_table1_crrt.csv
    │   ├── {site}_table1_crrt.html
    │   ├── {site}_table1_crrt_long.csv
    │   ├── {site}_missingness_summary.csv
    │   ├── {site}_crrt_settings_*.csv
    │   └── graphs/                  # Site-prefixed PNGs + CSVs
    │
    ├── psm_iptw/                    # Point-treatment causal analysis
    │   ├── {site}_Table1_unadjusted.{html,csv}
    │   ├── {site}_TableS1_matched.{html,csv}
    │   ├── {site}_TableS2_IPTW.{html,csv}
    │   ├── {site}_IPTW_pooled_results.csv
    │   ├── {site}_ModelComparison_PSMvsIPTW.csv
    │   ├── {site}_subgroup_analysis_results.csv
    │   ├── {site}_PSM_CIF_data.csv  # CIF curve data for multi-site combining
    │   ├── {site}_IPTW_CIF_data.csv # CIF curve data for multi-site combining
    │   ├── {site}_dose_response_*.csv
    │   ├── {site}_evalue_sensitivity.csv
    │   └── {site}_*.png             # PS overlap, love plots, CIF curves, forest plots
    │
    └── msm/                         # Time-varying MSM analysis
        ├── msm_output_guide.md
        ├── msm_sensitivity_findings.md
        ├── time_varying/            # Primary (12h intervals)
        │   ├── {site}_CRRT_30cutoff_Table1_unadjusted.{html,csv}
        │   ├── {site}_CRRT_30cutoff_TableS1_MSM_IPTW.{html,csv}
        │   ├── {site}_CRRT_30cutoff_MSM_CIF_data.csv
        │   ├── {site}_CRRT_30cutoff_MSM_*.csv
        │   └── {site}_CRRT_30cutoff_*.png
        └── time_varying_sensitivity/ # Sensitivity (24h intervals)
            ├── (same structure as time_varying/)
            └── ...
```

## Key Modules

| Module | Description |
|--------|-------------|
| `code/sofa_calculator.py` | Polars-based SOFA score computation with timezone normalization |
| `code/utils.py` | CRRT outlier handling using `config/outlier_config.json` thresholds |
| `code/pipeline_helpers.py` | Config validation, intermediate file loading, safe CLIF table loading |
| `utils/config.py` | Loads `config/config.json` and exposes as module-level `config` dict |
