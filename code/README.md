# Pipeline Steps

Per-step descriptions for the CRRT epidemiology pipeline. Output files are documented in [`../output/README.md`](../output/README.md); setup and run instructions are in the [top-level README](../README.md).

### Descriptive Epidemiology (Steps 00-03)

#### Step 00: Cohort Identification (`00_cohort.py`)
Identifies the CRRT cohort with inclusion/exclusion criteria (excludes ESRD via ICD codes), computes CRRT initiation times, outcomes (mortality, LOS), and IMV duration. Generates a CONSORT flow diagram.

#### Step 01: Wide Dataset (`01_create_wide_df.py`)
Builds a wide longitudinal dataset merging labs, vitals, medications, respiratory support, CRRT therapy, and ADT into a single time-indexed table. Includes vasopressor unit conversion and norepinephrine-equivalent (NEE) computation.

#### Step 02: Table 1 (`02_construct_crrt_tableone.py`)
Builds the per-encounter analytic cohort (`tableone_analysis_df.parquet`) and generates Table 1 of baseline characteristics (demographics, SOFA total and **non-renal SOFA**, nearest-measured baseline labs, vasopressors, IMV, net UF, outcomes) at CRRT initiation (−12h to +3h), **stratified by initial delivered-dose band** (<20 / 20–30 / >30 mL/kg/hr). Computes SOFA via `sofa_calculator.py` and a sepsis (ASE) flag.

#### Step 03: CRRT Epidemiology (`03_crrt_epidemiology.py`)
Comprehensive per-site descriptive epidemiology: (A) CRRT incidence/utilization by ICU population and co-indication, (B) practice variation + dose-band quality, (C) distribution figures (incidence, delivered dose, net-UF intensity Murugan bands, mortality-vs-UF), and (D) 30-day longitudinal trajectories (delivered dose, **net ultrafiltration: rate over time + cumulative volume + course-average intensity**, labs, vasopressor/NEE, IMV-state, CRRT-state). Net UF is summarized over the whole course (the first-3h initiation snapshot was retired); MAP and respiratory trajectories were cut.

#### Step 03b: Risk-Standardized Mortality (`03b_crrt_epi_smr.py`)
Computes the per-site standardized mortality ratio (SMR = observed / case-mix-expected 30-day deaths) on the **same analytic cohort as Table 1**, applying a frozen reference model (`config/smr_reference_model.json`, developed on MIMIC-IV) shipped with the repo. Covariates: age, sex, SOFA, baseline lactate, Charlson index. Reports the SMR with an exact-Poisson (Byar) CI plus transfer calibration. **Sites need only the tracked model JSON — no MIMIC data.**

#### Step 06: Low-Dose Characterization (`06_low_dose_characterization.py`)
Descriptive characterization of the very-low-dose CRRT subcohort (delivered dose 10–15 mL/kg/hr): dose-band tally + very-low-vs-rest baseline comparison.

### Causal Inference (Steps 04 / 05 / 05b)

> These steps require the descriptive pipeline (steps 00-03) to have completed first.

#### Step 04: Competing Risk DataFrame & Causal CONSORT (`04_build_causal_df.py`)
Builds a wide competing-risk DataFrame (58 columns) with labs, SOFA scores, oxygenation, vasopressors, IMV, and Charlson Comorbidity Index at t=0/12/24h windows. Includes sensitivity columns for 24h/48h interval analyses. Generates a missingness heatmap and a CONSORT flow diagram showing cohort narrowing from descriptive analysis through causal eligibility criteria.

#### Step 05: PSM & IPTW (`05_PSM_IPTW_CRRT_dose.R`)
Point-treatment causal analysis of CRRT dose (high >=30 vs low <30 mL/kg/hr) on 30-day mortality. Three branches:
- **PSM**: Nearest-neighbor matching + Fine-Gray competing risk models
- **IPTW**: Super Learner propensity scores + cause-specific Cox models
- **Subgroup analysis**: 9 clinically motivated subgroups with interaction tests

Uses MICE imputation (5 datasets), Rubin's rules pooling, bootstrap CIF curves (500 reps), and E-value sensitivity analysis.

#### Step 05b: Dose-Response Analysis (`05b_dose_response_analysis.R`)
Continuous dose-response analysis with three approaches:
- **Dose decile**: Unadjusted mortality by dose decile
- **Natural spline Cox**: Covariate-adjusted nonlinear dose-response (nonlinearity test)
- **Generalized propensity score**: Doubly-robust causal analysis via CBPS (continuous treatment)

Also generates a target trial emulation specification table and a 24h exclusion sensitivity analysis.

> **Archived:** the time-varying marginal structural models (`06_time_varying_MSM.R` / `06b_time_varying_MSM_sensitivity.R`) and the NCT06021288 trial-emulation script (`05c_low_dose_emulation.R`) were cut from the manuscript and moved to `archive/code/`. They are recoverable via the `msm-v1` / `low-dose-emulation-v1` git tags but are no longer part of the pipeline.
