# CLIF Epidemiology of CRRT

Multi-site analysis of Continuous Renal Replacement Therapy (CRRT) for Acute Kidney Injury (AKI) across the CLIF consortium.

**CLIF Version:** 2.1.0

## Objectives

This is a **descriptive-led CLIF-value study**: the headline is CRRT descriptive epidemiology and risk-standardized mortality, with a point-treatment causal analysis of CRRT dose as one rigorous section.

1. **Descriptive epidemiology** of CRRT for AKI across CLIF consortium hospitals:
   - Incidence/frequency of CRRT (denominator: all ICU admissions without ESRD — a population denominator registries cannot compute)
   - Pre-CRRT clinical state (BMP, lactate, vasopressors, IMV, location) and 30-day longitudinal trajectories
   - Treatment details (CRRT mode, settings, delivered dose, net ultrafiltration, duration)
   - Outcomes (30-day and in-hospital mortality, LOS, IMV duration)

2. **Risk-standardized mortality (SMR)** — observed vs case-mix-expected 30-day mortality per site, against a reference model developed on an external CLIF dataset (MIMIC-IV); the CRRTnet-style "practice varies but risk-adjusted mortality does not" benchmark.

3. **Causal inference** — effect of CRRT dose on 30-day mortality:
   - Point-treatment analysis (PSM with Fine-Gray competing risk; IPTW via Super Learner + cause-specific Cox)
   - Continuous dose-response (natural spline Cox, generalized propensity score) — sensitivity
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
| **clif_crrt_therapy** | `hospitalization_id`, `recorded_dttm`, `crrt_mode_category`, `blood_flow_rate`, `dialysate_flow_rate`, `pre_/post_filter_replacement_fluid_rate`, `ultrafiltration_out` | - |
| **clif_hospital_diagnosis** | `hospitalization_id`, `diagnosis_code`, `present_on_admission` | - |
| **clif_microbiology_culture** | `hospitalization_id`, `collected_dttm`, `result_category` | *(optional)* used for the sepsis (ASE) flag in Table 1; if absent the pipeline sets the sepsis flag to NA and continues |

See `config/clif_data_requirements.yaml` for the full column and category specification.

## Prerequisites

- **Python 3.11+**
- **R 4.x** (for causal inference scripts 05 / 05b)
- **UV package manager** ([install here](https://docs.astral.sh/uv/))
- Access to CLIF 2.1.0 data tables at your site

## Configuration

Edit `config/config.json` (copied from `config/config_template.json`):

```json
{
    "site_name": "Your_Site_Name",
    "tables_path": "/path/to/clif/tables/",
    "file_type": "parquet",
    "timezone": "US/Central",
    "project_root": "/path/to/CLIF-epidemiology-of-CRRT",
    "output_dir": "output",
    "external_dev_site": false
}
```

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

# 4. Run the full pipeline (descriptive + SMR 00-03b + 06, causal 04-05b)
bash run_pipeline.sh
# If it prints "R script(s) failed" / "NOT RUN", run code/05_PSM_IPTW_CRRT_dose.R
# then code/05b_dose_response_analysis.R manually in RStudio.
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
Rscript -e "renv::restore(prompt = FALSE)"

:: 4. Run the full pipeline (descriptive + SMR 00-03b + 06, causal 04-05b)
.\run_pipeline.bat
:: If it prints "R did not complete", run code\05_PSM_IPTW_CRRT_dose.R
:: then code\05b_dose_response_analysis.R manually in RStudio.
```

> **Note:** `run_pipeline.sh` runs all steps sequentially: Python scripts `00 01 02 03 03b 06 04` followed by R scripts `05 05b`. R scripts are invoked with `--no-init-file` to avoid conflicts with user `.Rprofile` settings; CRAN mirror defaults to `https://cloud.r-project.org` if not configured. Use `bash run_pipeline.sh --descriptive-only` to run just the descriptive + SMR deliverable (`00 01 02 03 03b 06`) and skip the causal R stack. The multi-site dashboard/manuscript scripts (`07`, `08`) are **coordinator-only** — run at the pooling site after collecting every site's `output/final_no_phi/`, not by individual sites.

### R Packages (`renv`)

Package versions are pinned in `renv.lock` it works inside RStudio:

1. Open the project folder in **RStudio** (it auto-activates renv via `.Rprofile`).
2. Run once in the R console: `renv::restore()`  — installs the pinned packages into a project-local library.
3. Then run the causal scripts (`code/05_PSM_IPTW_CRRT_dose.R`, `code/05b_dose_response_analysis.R`).

Pinned dependencies include `tidyverse`, `arrow`, `survival`, `cmprsk`, `MatchIt`, `WeightIt`, `cobalt`, `SuperLearner`, `randomForest`, `xgboost`, `gam`, `EValue`, `survminer`, `survey`, `mice`, `gtsummary`, `patchwork` (full list in `renv.lock`). If a site doesn't use renv, the scripts still auto-install missing CRAN packages as a fallback.

## Pipeline Steps

Per-step descriptions live in [`code/README.md`](code/README.md). Output files are documented in [`output/README.md`](output/README.md).

## Output Structure

Results split into **`output/final_no_phi/`** (aggregate, shareable) and **`output/intermediate_phi/`** (patient-level, never share). See [`output/README.md`](output/README.md) for the full tree and data-security rules.

## Key Modules

| Module | Description |
|--------|-------------|
| `code/sofa_calculator.py` | Polars-based SOFA score computation with timezone normalization |
| `code/utils.py` | CRRT outlier handling using `config/outlier_config.json` thresholds |
| `code/pipeline_helpers.py` | Config validation, intermediate file loading, safe CLIF table loading |
| `utils/config.py` | Loads `config/config.json` and exposes as module-level `config` dict |
