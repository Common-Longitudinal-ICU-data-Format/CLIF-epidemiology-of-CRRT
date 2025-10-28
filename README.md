# CLIF Epidemiology of CRRT

This project studies the epidemiology of Continuous Renal Replacement Therapy (CRRT) for Acute Kidney Injury (AKI) in the CLIF consortium using competing risk analysis.

## CLIF VERSION

2.1.0

## Objective

We hypothesize that among patients with acute kidney injury (AKI) who require continuous renal replacement therapy (CRRT), a higher initial CRRT dose is independently associated with increased 28-day mortality after adjustment for potential confounding factors.

**Exposure Definition:**  
The exposure is the initial CRRT dose, measured in mL/kg/hr. This is defined as:  
  [dialysate flow rate (mL/hr) + fluid replacement rate (mL/hr; both pre- and post-filter)] divided by weight (kg),  
where weight is the most recent value recorded prior to CRRT initiation.

- For CVVHD, the CRRT dose equals the dialysate flow rate alone (mL/hr/kg).
- For CVVH, the CRRT dose equals the fluid replacement rate alone (mL/hr/kg).
- For CVVHDF, the CRRT dose is the sum of dialysate and fluid replacement rates divided by weight, as described above.

## Project Structure

```
CLIF-epidemiology-of-CRRT/
├── code/
│   ├── 00_cohort.ipynb               # Step 1: Cohort identification
│   ├── 01_process_crrt.ipynb         # Step 2: CRRT processing & competing risk dataset
│   ├── 02_model.R                    # Step 3: Competing risk models (R)
│   ├── utils.py                      # Python utility functions
│   ├── sofa_calculator.py            # SOFA score calculation
│   └── README_ANALYSIS.md            # Detailed analysis documentation
├── config/
│   ├── config.json                   # Site-specific configuration (YOU MUST CREATE THIS)
│   ├── config_template.json          # Template for config.json
│   └── outlier_config.json           # Outlier thresholds for CRRT parameters
├── output/                           # Generated during analysis
│   ├── final/                        # Final results (aggregate data, safe to share)
│   └── intermediate/                 # Intermediate files (patient-level, DO NOT SHARE)
├── pyproject.toml                    # UV project dependencies
├── uv.lock                           # UV lock file
└── README.md
```

## Required CLIF Tables and Fields

Please refer to the online [CLIF data dictionary](https://clif-consortium.github.io/website/data-dictionary.html), [ETL tools](https://github.com/clif-consortium/CLIF/tree/main/etl-to-clif-resources), and [specific table contacts](https://github.com/clif-consortium/CLIF?tab=readme-ov-file#relational-clif) for more information on constructing the required tables and fields.

| Table Name | Required Variables | Required Categories |
| --- | --- | --- |
| **clif_patient** | `patient_id`, `race_category`, `ethnicity_category`, `sex_category`, `death_dttm` | - |
| **clif_hospitalization** | `patient_id`, `hospitalization_id`, `admission_dttm`, `discharge_dttm`, `age_at_admission` | - |
| **clif_adt** |  `hospitalization_id`, `hospital_id`,`in_dttm`, `out_dttm`, `location_category`, `location_type` | - |
| **clif_vitals** | `hospitalization_id`, `recorded_dttm`, `vital_category`, `vital_value` | heart_rate, resp_rate, sbp, dbp, map, spo2, weight_kg, height_cm |
| **clif_labs** | `hospitalization_id`, `lab_result_dttm`, `lab_category`, `lab_value` | sodium, potassium, chloride, bicarbonate, bun, creatinine, glucose_serum, calcium_total, lactate, magnesium, ph_arterial, ph_venous, po2_arterial |
| **clif_medication_admin_continuous** | `hospitalization_id`, `admin_dttm`, `med_name`, `med_category`, `med_dose`, `med_dose_unit` | norepinephrine, epinephrine, phenylephrine, vasopressin, dopamine, angiotensin, dobutamine, milrinone, isoproterenol |
| **clif_respiratory_support** | `hospitalization_id`, `recorded_dttm`, `device_category`, `mode_category`, `tracheostomy`, `fio2_set`, `lpm_set`, `resp_rate_set`, `peep_set`, `resp_rate_obs`, `tidal_volume_set`, `pressure_control_set`, `pressure_support_set`, `peak_inspiratory_pressure_set`, `tidal_volume_obs` | - |
| **clif_crrt_therapy** | `hospitalization_id`, `recorded_dttm`, `crrt_mode_name`, `crrt_mode_category`, `device_id`, `blood_flow_rate`, `dialysate_flow_rate`, `ultrafilteration_out` | - |
| **clif_hospital_diagnosis** | `hospitalization_id`, `diagnosis_code`, `present_on_admission` | - |

## Cohort Identification

### Inclusion Criteria
- Adults (age ≥ 18 years)
- Admissions between 2018-01-01 and 2024-12-31
- Received CRRT therapy
- Data completeness- Must have weight & CRRT settings  documented

### Exclusion Criteria
- Prior End-Stage Renal Disease (ESRD) diagnosis with any of the following codes:
  ```python
  esrd_codes = [
      'Z992',    # Dependence on renal dialysis
      'Z9115',   # Patient's noncompliance with renal dialysis
      'I120',    # Hypertensive chronic kidney disease with stage 5 CKD or ESRD
      'N186',    # End stage renal disease
      'I132',    # Hypertensive heart and chronic kidney disease with heart failure and ESRD
      'Z992',    # Dependence on renal dialysis (alternate code)
      'N186',    # End stage renal disease (alternate code)
      'I120',    # Hypertensive chronic kidney disease with stage 5 CKD or ESRD (alternate code)
      'Z91158',  # Patient's noncompliance with renal dialysis (alternate code)
      'I1311',   # Hypertensive heart and chronic kidney disease with heart failure and stage 5 CKD
      'I132',    # Hypertensive heart and chronic kidney disease with ESRD (alternate code)
      '5856',     #ICD9 :End stage renal disease
      '40391',    #ICD9: Hypertensive chronic kidney disease, unspecified, with chronic kidney disease stage V or end stage renal disease
      '40311',     #ICD9: Hypertensive chronic kidney disease, benign, with chronic kidney disease stage V or end stage renal disease
      'V4511',     #ICD9: Renal dialysis status
      'V4512'     #ICD9: Noncompliance with renal dialysis
  ]
  ```

## Setup Instructions

### Prerequisites

- **Python 3.11+**
- **UV package manager** ([install from here](https://docs.astral.sh/uv/))
- **R 4.0+** with packages: `cmprsk`, `survival`, `arrow`, `jsonlite`
- Access to CLIF 2.1.0 data tables at your site

### Step 1: Configuration

Create your site-specific configuration file:

```bash
cp config/config_template.json config/config.json
```

Edit `config/config.json` with your site settings:

```json
{
  "site_name": "YOUR_SITE_NAME",
  "tables_path": "/path/to/your/clif/tables",
  "file_type": "parquet",
  "timezone": "US/Central"
}
```

### Step 2: Environment Setup

Install dependencies using UV:

```bash
# Sync dependencies from pyproject.toml
uv sync

# Install Jupyter kernel
uv run python -m ipykernel install --user --name=CLIF-CRRT --display-name="CLIF-CRRT"
```

This will:
- Create a virtual environment (`.venv/`)
- Install all Python dependencies
- Set up Jupyter kernel for notebooks

## Usage

Run the analysis in 3 sequential steps:

### Step 1: Cohort Identification (`00_cohort.ipynb`)

**Purpose**: Identify cohort and calculate baseline characteristics
---

### Step 2: CRRT Processing & Competing Risk Dataset (`01_process_crrt.ipynb`)

**Purpose**: Process CRRT data and create competing risk analysis dataset

---

### Step 3: Competing Risk Models (`02_model.R`)

**Purpose**: Fit Fine-Gray competing risk regression models

## Expected Results

The analysis generates comprehensive results saved in the `output/final/` directory. Please upload to box. 




