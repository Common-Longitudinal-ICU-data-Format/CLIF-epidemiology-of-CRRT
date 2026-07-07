## Output directory

Results are split into two folders by **data-sharing safety** — the folder name tells you whether the
contents may leave your site:

* **`final_no_phi/`** — **aggregate, shareable** results only. Everything delivered to the project PI /
  consortium goes here: no `patient_id` or row-level records, every reported statistic at cell size
  **n ≥ 10**, site-prefixed files (e.g. `UCMC_*`). For multi-site pooling each site delivers its
  `final_no_phi/` tree, collected under `output/multi_site/<SITE>/final_no_phi/`.

* **`intermediate_phi/`** — **patient-level working data (NEVER share).** Cohort/analytic parquet files
  and filtered CLIF tables live here. **Never upload to Box or send to the PI.**

All of `output/` is git-ignored (regenerable by re-running the pipeline).

```
output/
├── intermediate_phi/                # patient-level parquet — NEVER SHARE
│   ├── outcomes_df.parquet
│   ├── index_crrt_df.parquet
│   ├── crrt_initiation.parquet
│   ├── wide_df.parquet
│   ├── tableone_analysis_df.parquet # analytic cohort (Table 1 + SMR source)
│   ├── {site}_smr_cohort.parquet    # row-level SMR cohort (03b)
│   ├── causal_df.parquet
│   └── filtered_clif/               # cohort-filtered CLIF tables
│
└── final_no_phi/                    # aggregate, shareable — site-prefixed (e.g. UCMC_*)
    ├── crrt_epi/                    # descriptive epidemiology + SMR
    │   ├── {site}_strobe_counts.csv
    │   ├── {site}_table1_crrt.{csv,html}
    │   ├── {site}_crrt_incidence.csv
    │   ├── {site}_crrt_practice_quality.csv
    │   ├── {site}_smr.csv                 # 03b
    │   ├── {site}_smr_calibration.csv     # 03b
    │   └── graphs/                  # trajectory + distribution CSVs/PNGs
    │
    ├── low_dose/                    # very-low-dose subcohort (06)
    │   └── {site}_low_dose_{counts,long,table}.csv
    │
    ├── diagnostics/                 # internal QC (00/04): missingness, settings, start-location
    │
    └── psm_iptw/                    # point-treatment causal analysis (04/05/05b)
        ├── {site}_table2_unadjusted_balance.{html,csv}
        ├── {site}_TableS1_matched.{html,csv}
        ├── {site}_TableS2_IPTW.{html,csv}
        ├── {site}_IPTW_pooled_results.csv
        ├── {site}_fg_psm_pooled_results.csv
        ├── {site}_subgroup_analysis_results.csv
        ├── {site}_PSM_CIF_data.csv  # CIF curve data for pooling
        ├── {site}_dose_decile_mortality.csv, {site}_dose_response_rcs.csv, {site}_gps_dose_response.csv
        ├── {site}_evalue_sensitivity.csv
        ├── {site}_causal_consort_diagram.{png,pdf}   # 04
        └── {site}_*.png             # PS overlap, love plots, CIF curves, forest plots
```
