# CLIF Epidemiology of CRRT

This project studies the epidemiology of Continuous Renal Replacement Therapy (CRRT) for Acute Kidney Injury (AKI) using the Common Longitudinal ICU data Format (CLIF).

## CLIF VERSION 

2.1.0

## Objective

Describe the epidemiology of CRRT for AKI in the CLIF consortium, including:

- **Frequency** of first CRRT mode (denominator: all ICU admissions without ESRD)
- **Pre-CRRT conditions**: Last recorded values of labs, vitals, vasopressors, mechanical ventilation, and location
- **Treatment details**: Settings and duration by CRRT mode
- **Outcomes**: Mode switches and in-hospital mortality

## Project Structure

```
CLIF-epidemiology-of-CRRT/
├── code/
│   ├── 01_cohort_identification.py    # Cohort identification script
│   ├── 02_analysis_summary.py                 # Analysis script
│   ├── pyCLIF.py                      # Utility functions
│   ├── sofa_score.py                  # module for sofa calculation
│   ├── waterfall.py                   # module for waterfall algorithms for respiratory support and crrt therapy tables 
│   └── README.md
├── config/
│   ├── config_template.json
│   └── outlier_config.json           # Outlier thresholds
├── output/                           # Generated during analysis
│   ├── final/                        # Final results (aggregate data only)
│   │   └── graphs/
│   └── intermediate/                 # Intermediate files
├── requirements.txt                  # Python dependencies
├── setup_venv.sh                    # Virtual environment setup script
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
| **clif_labs** | `hospitalization_id`, `lab_result_dttm`, `lab_category`, `lab_value` | sodium, potassium, chloride, bicarbonate, bun, creatinine, glucose_serum, calcium_total, lactate, magnesium, ph_arterial, ph_venous |
| **clif_medication_admin_continuous** | `hospitalization_id`, `admin_dttm`, `med_name`, `med_category`, `med_dose`, `med_dose_unit` | norepinephrine, epinephrine, phenylephrine, vasopressin, dopamine, angiotensin, dobutamine, milrinone, isoproterenol |
| **clif_respiratory_support** | `hospitalization_id`, `recorded_dttm`, `device_category`, `mode_category`, `tracheostomy`, `fio2_set`, `lpm_set`, `resp_rate_set`, `peep_set`, `resp_rate_obs`, `tidal_volume_set`, `pressure_control_set`, `pressure_support_set`, `peak_inspiratory_pressure_set`, `tidal_volume_obs` | - |
| **clif_crrt_therapy** | `hospitalization_id`, `recorded_dttm`, `crrt_mode_name`, `crrt_mode_category`, `device_id`, `blood_flow_rate`, `dialysate_flow_rate`, `ultrafilteration_out` | - |
| **clif_hospital_diagnosis** | `hospitalization_id`, `diagnosis_code` | - |

## Cohort Identification

### Inclusion Criteria
- Adults (age ≥ 18 years)
- Admissions between 2018-01-01 and 2024-12-31
- Received CRRT therapy

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

### 1. Prerequisites

- Python 3.8 or higher
- Access to CLIF data tables at your site or mimic

### 2. Environment Setup

Run the setup script to create a virtual environment and install dependencies:

**For Unix/macOS:**
```bash
./setup_venv.sh
```

**For Windows:**
```cmd
setup_venv.bat
```

This script will:
- Create a Python virtual environment
- Install all required packages
- Set up Jupyter kernel
- Create necessary directories

### 3. Configuration

Update the configuration file for your site:

```bash
cp config/config_template.json config/config.json
```

Edit `config/config.json` with your site-specific settings:
- `site_name`: Your site identifier
- `tables_path`: Path to your CLIF data files
- `file_type`: "parquet" or "csv"
- `timezone`: Your site's timezone

## Usage

### Activate Virtual Environment

**For Unix/macOS:**
```bash
source venv/bin/activate
```

**For Windows:**
```cmd
venv\Scripts\activate.bat
```

### Run Analysis

1. **Cohort Identification** (run first):
```bash
python code/01_cohort_identification.ipynb
```

2. **Analysis** (run after cohort identification):
```bash
python code/02_analysis_summary.pynb
```


## Expected Results

The analysis generates comprehensive results saved in the `output/final/` directory. Upload them to the project box folder. 

## Data Privacy

- **No patient-level data** is saved in `output/final/`
- All final outputs contain only aggregate statistics
- Patient-level data remains in `output/intermediate/` for internal processing only. DO NOT UPLOAD TO BOX. 

## Authors

- Kaveri Chhikara

---


