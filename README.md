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
| **clif_crrt_therapy** | `hospitalization_id`, `recorded_dttm` (+ settings columns if `has_crrt_settings=true`) | - |
| **clif_hospital_diagnosis** | `hospitalization_id`, `diagnosis_code`, `present_on_admission` | - |
| **clif_microbiology_culture** | `hospitalization_id`, `collected_dttm`, `result_category` | *(optional)* used for the sepsis (ASE) flag in Table 1; if absent the pipeline sets the sepsis flag to NA and continues |

See `config/clif_data_requirements.yaml` for the full column and category specification.

## Prerequisites

- **Python 3.11+**
- **R 4.x** (for causal inference scripts 05 / 05b)
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

# 4. Run the full pipeline (descriptive + SMR 00-03b + 06, causal 04-05b)
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

:: 4. Run the full pipeline (descriptive + SMR 00-03b + 06, causal 04-05b)
.\run_pipeline.bat
```

> **Note:** `run_pipeline.sh` runs all steps sequentially: Python scripts `00 01 02 03 03b 06 04` followed by R scripts `05 05b`. R scripts are invoked with `--no-init-file` to avoid conflicts with user `.Rprofile` settings; CRAN mirror defaults to `https://cloud.r-project.org` if not configured. Use `bash run_pipeline.sh --descriptive-only` to run just the descriptive + SMR deliverable (`00 01 02 03 03b 06`) and skip the causal R stack. The multi-site dashboard/manuscript scripts (`07`, `08`) are **coordinator-only** — run at the pooling site after collecting every site's `output/final/`, not by individual sites.

### R Packages

The R scripts auto-install missing CRAN packages on first run. Key dependencies include:

`tidyverse`, `arrow`, `survival`, `cmprsk`, `MatchIt`, `WeightIt`, `cobalt`, `SuperLearner`, `randomForest`, `xgboost`, `gam`, `EValue`, `survminer`, `survey`, `mice`, `gtsummary`, `patchwork`, `splines`

### Windows workstations (managed / locked-down)

On managed Windows machines, R packages auto-installed mid-run can be blocked by **Windows Application Control / Smart App Control** while *loading* their DLLs, e.g.:

```
unable to load shared object '...\R\win-library\4.5\utf8\libs\x64\utf8.dll':
  LoadLibrary failure: An Application Control policy has blocked this file.
```

Fix it with a **one-time setup** before running the pipeline (installs packages, then strips the "Mark of the Web" that triggers the block):

```powershell
powershell -ExecutionPolicy Bypass -File .\setup_r.ps1
```

Or do the two steps by hand once: install the packages in R/RStudio, then
`Get-ChildItem "$env:LOCALAPPDATA\R\win-library" -Recurse -File | Unblock-File`.

After that, `.\run_pipeline.bat` runs end to end.

### Running the R (causal) scripts manually

The causal stack is just two scripts — you can always run them by hand from the project root (this is the most transparent way to see any error, and mirrors exactly what the pipeline does):

```
Rscript --no-init-file code/05_PSM_IPTW_CRRT_dose.R
Rscript --no-init-file code/05b_dose_response_analysis.R
```

They read `output/intermediate/causal_df.parquet` (produced by `04_build_causal_df.py`) and write to `output/final/psm_iptw/`. So the minimal site workflow is: `run_pipeline.(sh|bat) --descriptive-only` for the descriptive + SMR deliverable, then run these two R scripts for the causal results. On a machine where R packages are already installed and loadable, `run_pipeline` full mode does all of this in one command.

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
Builds the per-encounter analytic cohort (`tableone_analysis_df.parquet`) and generates Table 1 of baseline characteristics (demographics, SOFA total and **non-renal SOFA**, nearest-measured baseline labs, vasopressors, IMV, net UF, outcomes) at CRRT initiation (−12h to +3h), **stratified by initial delivered-dose band** (<20 / 20–30 / >30 mL/kg/hr). Computes SOFA via `sofa_calculator.py` and a sepsis (ASE) flag.

**Outputs:**
- `output/intermediate/tableone_analysis_df.parquet`
- `output/final/crrt_epi/{site}_table1_crrt.csv`
- `output/final/crrt_epi/{site}_table1_crrt.html`

#### Step 03: CRRT Epidemiology (`03_crrt_epidemiology.py`)
Comprehensive per-site descriptive epidemiology: (A) CRRT incidence/utilization by ICU population and co-indication, (B) practice variation + dose-band quality, (C) distribution figures (incidence, delivered dose, net-UF intensity Murugan bands, mortality-vs-UF), and (D) 30-day longitudinal trajectories (delivered dose, **net ultrafiltration: rate over time + cumulative volume + course-average intensity**, labs, vasopressor/NEE, IMV-state, CRRT-state). Net UF is summarized over the whole course (the first-3h initiation snapshot was retired); MAP and respiratory trajectories were cut.

**Outputs:** `output/final/crrt_epi/{site}_crrt_incidence.csv`, `{site}_crrt_practice_quality.csv`, and `output/final/crrt_epi/graphs/{site}_*.png` (+ aggregate `*.csv`)

#### Step 03b: Risk-Standardized Mortality (`03b_crrt_epi_smr.py`)
Computes the per-site standardized mortality ratio (SMR = observed / case-mix-expected 30-day deaths) on the **same analytic cohort as Table 1**, applying a frozen reference model (`config/smr_reference_model.json`, developed on MIMIC-IV) shipped with the repo. Covariates: age, sex, SOFA, baseline lactate, Charlson index. Reports the SMR with an exact-Poisson (Byar) CI plus transfer calibration. **Sites need only the tracked model JSON — no MIMIC data.**

**Outputs:** `output/final/crrt_epi/{site}_smr.csv`, `{site}_smr_calibration.csv`

#### Step 06: Low-Dose Characterization (`06_low_dose_characterization.py`)
Descriptive characterization of the very-low-dose CRRT subcohort (delivered dose 10–15 mL/kg/hr): dose-band tally + very-low-vs-rest baseline comparison. Requires `has_crrt_settings=true`.

**Outputs:** `output/final/low_dose/{site}_low_dose_{counts,long,table}.csv`

### Causal Inference (Steps 04 / 05 / 05b)

> These steps require the descriptive pipeline (steps 00-03) to have completed first.

#### Step 04: Competing Risk DataFrame & Causal CONSORT (`04_build_causal_df.py`)
Builds a wide competing-risk DataFrame (58 columns) with labs, SOFA scores, oxygenation, vasopressors, IMV, and Charlson Comorbidity Index at t=0/12/24h windows. Includes sensitivity columns for 24h/48h interval analyses. Generates a missingness heatmap and a CONSORT flow diagram showing cohort narrowing from descriptive analysis through causal eligibility criteria.

**Outputs:**
- `output/intermediate/causal_df.parquet`
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

> **Archived:** the time-varying marginal structural models (`06_time_varying_MSM.R` / `06b_time_varying_MSM_sensitivity.R`) and the NCT06021288 trial-emulation script (`05c_low_dose_emulation.R`) were cut from the manuscript and moved to `archive/code/`. They are recoverable via the `msm-v1` / `low-dose-emulation-v1` git tags but are no longer part of the pipeline.

## Output Structure

All of `output/` is gitignored (regenerable). For multi-site pooling, each site delivers its `output/final/` tree, collected under `output/multi_site/<SITE>/final/`.

```
output/
├── intermediate/                    # Intermediate parquet files
│   ├── outcomes_df.parquet
│   ├── index_crrt_df.parquet
│   ├── crrt_initiation.parquet
│   ├── wide_df.parquet
│   ├── tableone_analysis_df.parquet # analytic cohort (Table 1 + SMR source)
│   └── causal_df.parquet
│
├── smr/                             # SMR cohort parquet (03b; row-level)
│
└── final/                           # All files are site-prefixed (e.g., UCMC_*)
    ├── crrt_epi/                    # Descriptive epidemiology + SMR
    │   ├── {site}_strobe_counts.csv
    │   ├── {site}_table1_crrt.{csv,html}
    │   ├── {site}_crrt_incidence.csv
    │   ├── {site}_crrt_practice_quality.csv
    │   ├── {site}_smr.csv                 # 03b
    │   ├── {site}_smr_calibration.csv     # 03b
    │   └── graphs/                  # trajectory + distribution CSVs/PNGs
    │
    ├── low_dose/                    # Very-low-dose subcohort (06)
    │   └── {site}_low_dose_{counts,long,table}.csv
    │
    ├── diagnostics/                 # Internal QC (00/04): missingness, settings
    │
    └── psm_iptw/                    # Point-treatment causal analysis (04/05/05b)
        ├── {site}_table2_unadjusted_balance.{html,csv}
        ├── {site}_TableS1_matched.{html,csv}
        ├── {site}_TableS2_IPTW.{html,csv}
        ├── {site}_IPTW_pooled_results.csv
        ├── {site}_fg_psm_pooled_results.csv
        ├── {site}_subgroup_analysis_results.csv
        ├── {site}_PSM_CIF_data.csv  # CIF curve data for pooling (Fig: PSM CIF death)
        ├── {site}_dose_decile_mortality.csv, {site}_dose_response_rcs.csv, {site}_gps_dose_response.csv
        ├── {site}_evalue_sensitivity.csv
        ├── {site}_causal_consort_diagram.{png,pdf}   # 04
        └── {site}_*.png             # PS overlap, love plots, CIF curves, forest plots
```

> The MIMIC SMR-reference **development** tree lives separately in `output_mimic/` (dev/coordinator-only, gitignored). The pooled multi-site dashboard + manuscript artifacts (`output/multi_site/`) are produced by the coordinator scripts `07`/`08`.

## Key Modules

| Module | Description |
|--------|-------------|
| `code/sofa_calculator.py` | Polars-based SOFA score computation with timezone normalization |
| `code/utils.py` | CRRT outlier handling using `config/outlier_config.json` thresholds |
| `code/pipeline_helpers.py` | Config validation, intermediate file loading, safe CLIF table loading |
| `utils/config.py` | Loads `config/config.json` and exposes as module-level `config` dict |
