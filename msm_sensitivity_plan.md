# MSM Sensitivity Analysis: Wider Time Windows (0-24h / 24-48h)

**Date:** 2026-03-06
**Goal:** Build an alternative MSM dataframe with two 24-hour intervals (0-24h, 24-48h) instead of two 12-hour intervals (0-12h, 12-24h), to test whether wider windows improve or worsen the positivity violation documented in `msm_balance_diagnostics.md`.
**Motivation:** The current MSM uses 12h intervals and suffers from severe positivity violation at t=12h (27% of patients have PS > 0.95) driven by treatment persistence. This sensitivity analysis asks: *do wider windows — which average over more dose variability — reduce the near-deterministic A_lag → A relationship, or does the stickiness of CRRT dose persist regardless of window width?*

---

## 1. Overview of Changes

### Current Pipeline (12h intervals)
| Column | Description | Time Window |
|--------|-------------|-------------|
| `crrt_dose_0_12` | Mean dose 0-12h | 0-12h |
| `crrt_dose_12_24` | Mean dose 12-24h | 12-24h |
| `lactate_0`, `bicarbonate_0`, `potassium_0` | Labs at CRRT start | -12h to +3h |
| `lactate_12`, `bicarbonate_12`, `potassium_12` | Labs at 12h | +3h to +12h |
| `sofa_total_0`, `sofa_total_12` | SOFA scores | 0h, 0-12h |
| `oxygenation_index_0`, `_12` | Oxygenation | -12h/+3h, +3h/+12h |
| `norepinephrine_equivalent_0`, `_12` | Vasopressor burden | -12h/+3h, +3h/+12h |
| `imv_status_0`, `_12` | Mechanical ventilation | -12h/+3h, +3h/+12h |

### Sensitivity Pipeline (24h intervals)
| Column | Description | Time Window |
|--------|-------------|-------------|
| `crrt_dose_0_24` | Mean dose 0-24h | 0-24h |
| `crrt_dose_24_48` | Mean dose 24-48h | 24-48h |
| `lactate_0`, `bicarbonate_0`, `potassium_0` | Labs at CRRT start (unchanged) | -12h to +3h |
| `lactate_24`, `bicarbonate_24`, `potassium_24` | Labs at 24h | +12h to +24h |
| `sofa_total_0` | SOFA at baseline (unchanged) | 0h |
| `sofa_total_24` | SOFA at 24h | 0-24h window |
| `oxygenation_index_0` (unchanged), `_24` | Oxygenation at 24h | +12h to +24h |
| `norepinephrine_equivalent_0` (unchanged), `_24` | NEE at 24h | +12h to +24h |
| `imv_status_0` (unchanged), `_24` | IMV status at 24h | +12h to +24h |

---

## 2. New Columns to Add

All new columns will be appended to the existing dataframe (the original 12h columns remain intact). The sensitivity analysis script will select the appropriate columns.

### 2A. CRRT Dosing Columns

| Column | Type | Derivation | Fallback |
|--------|------|-----------|----------|
| `crrt_dose_0_24` | continuous (mL/kg/hr) | Mean of `crrt_dose_ml_kg_hr` for rows where `hours_from_crrt` in [0, 24) | None — should be available for all patients surviving 24h |
| `crrt_dose_24_48` | continuous (mL/kg/hr) | Mean of `crrt_dose_ml_kg_hr` for rows where `hours_from_crrt` in [24, 48) | Carry forward `crrt_dose_0_24` if no charted values in 24-48h window |

**Implementation note:** `crrt_dose_0_24` can be derived as the weighted or simple average of the existing `crrt_dose_0_12` and `crrt_dose_12_24` columns, but computing fresh from `wide_df` with hours_from_crrt in [0, 24) is more accurate (avoids equal-weighting two intervals with potentially different observation densities). Compute fresh from `crrt_rows`.

### 2B. Lab Columns at t=24h

| Column | Type | Extraction Window | Fallback |
|--------|------|-------------------|----------|
| `lactate_24` | continuous (mmol/L) | Last non-null value in +12h to +24h | `lactate_0` (baseline) |
| `bicarbonate_24` | continuous (mEq/L) | Last non-null value in +12h to +24h | `bicarbonate_0` (baseline) |
| `potassium_24` | continuous (mEq/L) | Last non-null value in +12h to +24h | `potassium_0` (baseline) |

**Window rationale:** The +12h to +24h window is symmetric with the current +3h to +12h window for the 12h analysis — it captures the "state at the interval boundary" using the last observation before that boundary. For the sensitivity analysis, the decision point is at 24h, so we want the most recent lab before 24h that isn't also before 12h (which would overlap the baseline interval).

### 2C. SOFA at t=24h

| Column | Type | Computation Window |
|--------|------|--------------------|
| `sofa_total_24` | continuous | Worst-value SOFA in the 0-24h window via `compute_sofa_polars()` |

**Note:** This uses the same `compute_sofa_polars()` helper as the current t=12h SOFA, just with `end_dttm = crrt_initiation_time + 24h`.

### 2D. Time-Varying Covariates at t=24h

| Column | Type | Extraction Window | Fallback |
|--------|------|-------------------|----------|
| `oxygenation_index_24` | continuous | Last non-null in +12h to +24h | `oxygenation_index_0` |
| `norepinephrine_equivalent_24` | continuous (mcg/kg/min) | Last non-null in +12h to +24h | `norepinephrine_equivalent_0` |
| `imv_status_24` | binary (0/1) | Last non-null in +12h to +24h | `imv_status_0` |

Same `extract_last_in_window()` pattern as the existing code, with window bounds changed to (12, 24).

**NEE imputation:** Same 3-step cascade as current code — patients with no vasopressor records get 0; then widen window; then impute 0 for remaining.

---

## 3. Eligibility for Sensitivity Cohort

The sensitivity analysis uses a **more restrictive** cohort than the primary analysis because patients must survive and remain on CRRT for 48h (vs 24h currently).

### Exclusion Criteria
1. **Died within 48h of CRRT initiation** (`outcome == 2` AND `time_to_event_90d <= 2.0`)
2. **CRRT duration < 48h** (`crrt_duration_days < 2.0`)
3. **SCUF-only** (no dialysate in 0-48h) — extend existing `has_dialysate_24h` check to 48h

### Implementation
- Add a boolean column `eligible_48h_sensitivity` (1/0) to the dataframe so both the primary and sensitivity cohorts coexist in a single parquet file
- The R sensitivity script filters on `eligible_48h_sensitivity == 1`
- Report the sample size reduction (N primary → N sensitivity) in the output

### Expected Sample Reduction
- The current primary analysis excludes patients who died or were off CRRT within 24h
- Extending to 48h will further reduce N; if N drops below ~200, the sensitivity analysis may be underpowered (flag this clearly in output)

---

## 4. Column Schema for Sensitivity Columns

These columns are **appended** to the existing 49-column schema. The output parquet retains all original columns plus these new ones.

| # | Column | Type | Notes |
|---|--------|------|-------|
| 50 | `crrt_dose_0_24` | continuous | Mean CRRT dose (mL/kg/hr) in 0-24h |
| 51 | `crrt_dose_24_48` | continuous | Mean CRRT dose (mL/kg/hr) in 24-48h |
| 52 | `lactate_24` | continuous | mmol/L, last obs in +12h to +24h |
| 53 | `bicarbonate_24` | continuous | mEq/L, last obs in +12h to +24h |
| 54 | `potassium_24` | continuous | mEq/L, last obs in +12h to +24h |
| 55 | `sofa_total_24` | continuous | Worst-value SOFA in 0-24h window |
| 56 | `oxygenation_index_24` | continuous | Last obs in +12h to +24h |
| 57 | `norepinephrine_equivalent_24` | continuous | mcg/kg/min, last obs in +12h to +24h |
| 58 | `imv_status_24` | binary (0/1) | IMV status, last obs in +12h to +24h |
| 59 | `eligible_48h_sensitivity` | binary (0/1) | 1 if patient survived + on CRRT >=48h |

**Total columns after expansion: 59**

---

## 5. Python Script Modifications

### File: `code/04_build_msm_competing_risk_df.py`
Changes are additive — no existing code is removed.

### Step 4 additions (CRRT dose):
```
# After existing dose_12_24 computation:

# 4d: Mean dose 0-24h
dose_0_24 = crrt_rows[(hours >= 0) & (hours < 24)]
    .groupby("encounter_block")["crrt_dose_ml_kg_hr"].mean()
    → rename to "crrt_dose_0_24"

# 4e: Mean dose 24-48h
dose_24_48 = crrt_rows[(hours >= 24) & (hours < 48)]
    .groupby("encounter_block")["crrt_dose_ml_kg_hr"].mean()
    → rename to "crrt_dose_24_48"
    → fallback: fill NaN with crrt_dose_0_24
```

### Step 5 additions (Labs at 24h):
```
# After existing labs_12h extraction:

labs_24h = labs_wide[(hours >= 12) & (hours <= 24)]

For each of lactate, bicarbonate, potassium:
    → last non-null in window → rename to *_24
    → fallback to *_0 (baseline)
```

### Step 6 additions (SOFA at 24h):
```
# After existing sofa_12h computation:

sofa_cohort_24h:
    start_dttm = crrt_initiation_time
    end_dttm   = crrt_initiation_time + 24h
    → compute_sofa_polars(...) → sofa_total_24
```

### Step 7 additions (Oxygenation, NEE, IMV at 24h):
```
# After existing _12 extraction:

For oxygenation_index, norepinephrine_equivalent, imv_status:
    t24 = extract_last_in_window(extra_df, col, 12, 24)
    → rename to *_24
    → fallback to *_0

NEE imputation at 24h: same cascade as current _12 logic
```

### Step 9b additions (Eligibility flag):
```
# After existing 24h exclusion:

eligible_48h = (
    NOT (outcome == 2 AND time_to_event_90d <= 2.0)
    AND crrt_duration_days >= 2.0
    AND has_dialysate_48h == 1  # extend dialysate check to 48h
)
result["eligible_48h_sensitivity"] = eligible_48h.astype(int)
```

### Step 10 (Final column ordering):
- Append 10 new columns to `final_cols` list

---

## 6. R Script: Sensitivity MSM

### File: `code/claude_time_varying_MSM_sensitivity.R`
Copy of `claude_time_varying_MSM.R` with the following changes:

### Data loading:
```r
df <- arrow::read_parquet(data_path)
df <- df %>% filter(eligible_48h_sensitivity == 1)
```

### Treatment naming convention

To avoid confusion with the primary analysis (which uses `A_0` for 0-12h and `A_12` for 12-24h), the sensitivity script uses **interval-encoded names**:

| Variable | Meaning | Primary equivalent |
|----------|---------|-------------------|
| `A_0_24` | Dose group for 0-24h interval | `A_0` (0-12h) |
| `A_24_48` | Dose group for 24-48h interval | `A_12` (12-24h) |

The generic long-format column `A` and `A_lag` remain unchanged (they are positional — `A` is "treatment at this time point", `A_lag` is "treatment at the prior time point").

### Treatment assignment (wider intervals):
```r
df_msm_wide <- df_complete %>%
  mutate(
    A_0_24  = as.integer(crrt_dose_0_24  >= dose_cutoff),   # 0-24h interval
    A_24_48 = as.integer(crrt_dose_24_48 >= dose_cutoff)    # 24-48h interval
  )
```

### Pivot to long (24h intervals):
```r
df_msm_long <- bind_rows(
  # t=0 row: baseline covariates + t=0 labs/severity
  df_msm_wide %>% transmute(
    ..., tpt = 0L, A = A_0_24,
    lactate = lactate_0, bicarbonate = bicarbonate_0, potassium = potassium_0,
    oxygenation_index = oxygenation_index_0,
    norepinephrine_equivalent = norepinephrine_equivalent_0,
    imv_status = imv_status_0,
    sofa_total = sofa_total_0, ...
  ),
  # t=24 row: 24h covariates
  df_msm_wide %>% transmute(
    ..., tpt = 24L, A = A_24_48,
    lactate = lactate_24, bicarbonate = bicarbonate_24, potassium = potassium_24,
    oxygenation_index = oxygenation_index_24,
    norepinephrine_equivalent = norepinephrine_equivalent_24,
    imv_status = imv_status_24,
    sofa_total = sofa_total_24, ...
  )
)
```

### Weighting loop:
```r
for (tt in c(0, 24)) { ... }  # Changed from c(0, 12)
```

### 4-group strategy labels:
```r
# Primary uses: High→High, High→Low, Low→High, Low→Low (for A_0/A_12)
# Sensitivity uses same labels but derived from A_0_24/A_24_48:
df_msm_wide <- df_msm_wide %>%
  mutate(
    strategy_4 = paste0(
      ifelse(A_0_24 == 1, "High", "Low"), " → ",
      ifelse(A_24_48 == 1, "High", "Low")
    )
  )
```

### Output directory:
```r
output_dir <- "output/final/time_varying_sensitivity"
```

---

## 7. Key Comparisons (Sensitivity vs Primary)

The sensitivity analysis report should include a side-by-side comparison table:

| Metric | Primary (12h) | Sensitivity (24h) |
|--------|--------------|-------------------|
| N (analysis sample) | — | — |
| % patients PS > 0.95 at t2 | 27.0% | ? |
| % patients PS > 0.99 at t2 | 5.0% | ? |
| ESS proportion | 90.2% | ? |
| A_lag SMD (weighted) | 0.80 | ? |
| Weight (kg) SMD (weighted) | -0.96 | ? |
| Lactate SMD (weighted) | 0.49 | ? |
| Dose switching rate (t1 → t2) | ~10.4% | ? |

### Hypotheses
- **Optimistic**: Wider windows average over more dose variability, diluting the treatment persistence signal. PS extremes improve, A_lag balance improves.
- **Pessimistic (more likely)**: CRRT dose is fundamentally sticky — the 0-24h mean and 24-48h mean will be nearly identical for most patients, and the switching rate will be similarly low. The positivity violation may even worsen because 48h survivors are more homogeneous (survivorship bias).
- **Either result is informative**: If wider windows don't help, it strengthens the conclusion that the positivity violation is an inherent limitation of MSM for CRRT dosing, not an artifact of interval width.

---

## 8. Implementation Order

1. **Python: Add new columns to `04_build_msm_competing_risk_df.py`**
   - CRRT dose (0-24h, 24-48h)
   - Labs at 24h
   - SOFA at 24h
   - Oxygenation, NEE, IMV at 24h
   - Eligibility flag (`eligible_48h_sensitivity`)
   - Update `final_cols` and save

2. **R: Create sensitivity MSM script**
   - Copy `claude_time_varying_MSM.R` → `claude_time_varying_MSM_sensitivity.R`
   - Swap 12h columns for 24h columns
   - Filter to `eligible_48h_sensitivity == 1`
   - Adjust time points in weighting loop and pivot

3. **Run and compare**
   - Execute Python script to regenerate parquet
   - Execute both R scripts (primary unchanged, sensitivity new)
   - Compare positivity metrics, balance diagnostics, and CIF estimates

---

## 9. Files Referenced

| File | Role |
|------|------|
| `code/intermediate/04_build_msm_competing_risk_df.py` | Python dataframe builder (modify) |
| `code/intermediate/01_create_wide_df.py` | Wide dataset source (read-only) |
| `code/claude_time_varying_MSM.R` | Primary MSM R script (read-only, copy for sensitivity) |
| `plans/msm_base/msm_analysis_df_schema.md` | Current 49-column schema |
| `plans/msm_base/msm_balance_diagnostics.md` | Positivity violation documentation |
| `output/intermediate/msm_competing_risk_df.parquet` | Output dataframe (will gain 10 columns) |
